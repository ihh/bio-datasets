#!/usr/bin/env python3
"""
precompile_panther.py — Pre-compile PANTHER cherry pairs into training shards.

Reads PANTHER .AN.fasta aligned ancestral sequences, selects cherry pairs
(closest by p-distance) plus a few medium-distance pairs, and writes
zstd-compressed JSONL shards compatible with tkf-mixdom's train_pfam.py.

Output format matches precompile_pairs.py: each record has keys
x (ancestor AA string), a (RLE alignment path), y (descendant AA string),
t (estimated evolutionary time), fam (family ID), id (deterministic hash).

Usage:
    # All training families:
    python fetch/panther/precompile_panther.py --split train

    # Quick test:
    python fetch/panther/precompile_panther.py \
        --families PTHR10000,PTHR10003 --out /tmp/test_panther/

    # Full with parallel processing:
    python fetch/panther/precompile_panther.py --split train --workers 4
"""

import argparse
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

_log_start = time.monotonic()


def _log(msg, end='\n'):
    elapsed = time.monotonic() - _log_start
    sys.stderr.write(f"[{elapsed:8.1f}s] {msg}" + end)
    sys.stderr.flush()


# ============================================================
# Constants and helpers (self-contained, no JAX dependency)
# ============================================================

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AA = len(AMINO_ACIDS)
AA_TO_IDX = {a: i for i, a in enumerate(AMINO_ACIDS)}
IDX_TO_AA = {i: a for i, a in enumerate(AMINO_ACIDS)}

# State codes
M_STATE, I_STATE, D_STATE = 1, 2, 3
_STATE_TO_CHAR = {M_STATE: 'M', I_STATE: 'I', D_STATE: 'D'}


def rle_encode(states):
    """Run-length encode M/I/D state list."""
    if not states:
        return ''
    parts = []
    current = states[0]
    count = 1
    for s in states[1:]:
        if s == current:
            count += 1
        else:
            parts.append(f'{_STATE_TO_CHAR[current]}{count}')
            current = s
            count = 1
    parts.append(f'{_STATE_TO_CHAR[current]}{count}')
    return ''.join(parts)


def parse_aligned_fasta(filepath):
    """Parse aligned FASTA. Returns list of (name, aligned_seq)."""
    names = []
    seqs = []
    current_name = None
    current_seq = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name is not None:
                    seqs.append(''.join(current_seq))
                    names.append(current_name)
                current_name = line[1:].split()[0]
                current_seq = []
            elif current_name is not None:
                current_seq.append(line)
        if current_name is not None:
            seqs.append(''.join(current_seq))
            names.append(current_name)
    return names, seqs


def ungap(seq):
    return ''.join(c for c in seq if c in AA_TO_IDX)


def p_distance(seq1, seq2):
    """Compute p-distance from aligned sequences (ignoring gap columns)."""
    matches = mismatches = 0
    for a, b in zip(seq1, seq2):
        if a in AA_TO_IDX and b in AA_TO_IDX:
            if a == b:
                matches += 1
            else:
                mismatches += 1
    total = matches + mismatches
    return mismatches / total if total > 0 else 1.0


def pdist_to_evo_time(pdist):
    if pdist >= 0.95:
        return 5.0
    corrected = 1.0 - pdist * (AA / (AA - 1.0))
    if corrected <= 0.01:
        return 5.0
    return -np.log(corrected)


def _aa_composition(seq):
    comp = np.zeros(AA)
    for c in seq:
        if c in AA_TO_IDX:
            comp[AA_TO_IDX[c]] += 1
    total = comp.sum()
    if total > 0:
        comp /= total
    return comp


def select_cherries(n_seqs, dist_fn, max_pairs=None):
    """Greedy nearest-neighbor cherry pairing."""
    if n_seqs < 2:
        return []
    dists = np.ones((n_seqs, n_seqs)) * 1e10
    for i in range(n_seqs):
        for j in range(i + 1, n_seqs):
            d = dist_fn(i, j)
            dists[i, j] = d
            dists[j, i] = d
    paired = set()
    pairs = []
    while len(paired) < n_seqs - 1:
        best_d = 1e10
        best_i, best_j = -1, -1
        for i in range(n_seqs):
            if i in paired:
                continue
            for j in range(i + 1, n_seqs):
                if j in paired:
                    continue
                if dists[i, j] < best_d:
                    best_d = dists[i, j]
                    best_i, best_j = i, j
        if best_i < 0:
            break
        pairs.append((best_i, best_j))
        paired.add(best_i)
        paired.add(best_j)
        if max_pairs and len(pairs) >= max_pairs:
            break
    return pairs


def select_medium_distance_pairs(n_seqs, dist_fn, rng, n_pairs=3,
                                 min_pdist=0.3, max_pdist=0.6):
    """Select a few pairs at medium evolutionary distance."""
    if n_seqs < 2:
        return []
    candidates = []
    # Sample random pairs and filter by distance
    max_attempts = min(n_seqs * 5, 200)
    for _ in range(max_attempts):
        i, j = rng.randint(0, n_seqs), rng.randint(0, n_seqs)
        if i == j:
            continue
        if i > j:
            i, j = j, i
        d = dist_fn(i, j)
        if min_pdist <= d <= max_pdist:
            candidates.append((i, j, d))
    # Deduplicate
    seen = set()
    unique = []
    for i, j, d in candidates:
        key = (i, j)
        if key not in seen:
            seen.add(key)
            unique.append((i, j))
    # Take up to n_pairs
    if len(unique) > n_pairs:
        indices = rng.choice(len(unique), n_pairs, replace=False)
        unique = [unique[idx] for idx in indices]
    return unique


def _aligned_pair_to_int_arrays(seq1, seq2):
    """Convert aligned string pair to gapped integer arrays (-1 = gap)."""
    a = np.array([AA_TO_IDX.get(c, -1) for c in seq1], dtype=np.int32)
    b = np.array([AA_TO_IDX.get(c, -1) for c in seq2], dtype=np.int32)
    mask = (a >= 0) | (b >= 0)
    return a[mask], b[mask]


def alignment_to_states(ancestor_chars, descendant_chars, gap_token=-1):
    """Convert aligned sequences to state sequence."""
    states = []
    anc_chars = []
    desc_chars = []
    for pos in range(len(ancestor_chars)):
        anc = int(ancestor_chars[pos])
        desc = int(descendant_chars[pos])
        has_anc = (anc != gap_token)
        has_desc = (desc != gap_token)
        if has_anc and has_desc:
            states.append(M_STATE)
            anc_chars.append(anc)
            desc_chars.append(desc)
        elif has_desc:
            states.append(I_STATE)
            desc_chars.append(desc)
        elif has_anc:
            states.append(D_STATE)
            anc_chars.append(anc)
    return states, anc_chars, desc_chars


def encode_pair(anc_seq, desc_seq, states, t_est, family_id,
                row1_name='', row2_name=''):
    """Encode pair as compact X/A/Y record dict."""
    x = ''.join(IDX_TO_AA[int(c)] for c in anc_seq)
    y = ''.join(IDX_TO_AA[int(c)] for c in desc_seq)
    a = rle_encode(states)
    id_str = f'{family_id}:{row1_name}:{row2_name}'
    pair_id = hashlib.sha256(id_str.encode()).hexdigest()[:16]
    return {
        'x': x, 'a': a, 'y': y,
        't': round(float(t_est), 6),
        'fam': family_id,
        'id': pair_id,
    }


# ============================================================
# Per-family processing
# ============================================================

def _process_family(msa_file, max_seqs=500):
    """Process one PANTHER .AN.fasta file.

    Selects cherry pairs + medium-distance pairs, encodes both directions.
    Returns list of record dicts.
    """
    fam = os.path.basename(msa_file).replace('.AN.fasta', '')
    names, seqs = parse_aligned_fasta(msa_file)
    n = len(seqs)
    if n < 2:
        return []

    # Filter: only keep sequences with at least 10 valid AAs
    valid_idx = [i for i in range(n) if len(ungap(seqs[i])) >= 10]
    if len(valid_idx) < 2:
        return []

    # Cap for cherry selection
    max_for_cherries = min(len(valid_idx), max_seqs)
    rng = np.random.RandomState(hash(fam) & 0x7FFFFFFF)
    if len(valid_idx) > max_for_cherries:
        subset = sorted(rng.choice(len(valid_idx), max_for_cherries, replace=False))
        idx_map = [valid_idx[s] for s in subset]
    else:
        idx_map = list(valid_idx)

    n_sub = len(idx_map)
    seqs_sub = [seqs[idx_map[i]] for i in range(n_sub)]

    # Distance function
    if n_sub <= 50:
        def dist_fn(i, j, _seqs=seqs_sub):
            return p_distance(_seqs[i], _seqs[j])
    else:
        comps = np.array([_aa_composition(ungap(s)) for s in seqs_sub])
        def dist_fn(i, j, _c=comps):
            return float(np.sum((_c[i] - _c[j]) ** 2))

    # Select cherry pairs
    cherries = select_cherries(n_sub, dist_fn)

    # Select medium-distance pairs
    # Use actual p_distance for medium pair selection
    def pdist_fn(i, j, _seqs=seqs_sub):
        return p_distance(_seqs[i], _seqs[j])
    medium_pairs = select_medium_distance_pairs(n_sub, pdist_fn, rng, n_pairs=3)

    all_pairs = cherries + medium_pairs

    records = []
    for idx1, idx2 in all_pairs:
        orig1, orig2 = idx_map[idx1], idx_map[idx2]
        aln_i, aln_j = _aligned_pair_to_int_arrays(seqs[orig1], seqs[orig2])
        if len(aln_i) == 0:
            continue

        pd = p_distance(seqs[orig1], seqs[orig2])
        t_est = pdist_to_evo_time(pd)

        states, anc_chars, desc_chars = alignment_to_states(aln_i, aln_j)
        if not states:
            continue

        # Encode both directions
        for ac, dc, r1, r2 in [
            (anc_chars, desc_chars, names[orig1], names[orig2]),
            (desc_chars, anc_chars, names[orig2], names[orig1]),
        ]:
            if ac is desc_chars:
                rev_states, rev_anc, rev_desc = alignment_to_states(aln_j, aln_i)
                rec = encode_pair(rev_anc, rev_desc, rev_states, t_est, fam, r1, r2)
            else:
                rec = encode_pair(ac, dc, states, t_est, fam, r1, r2)
            records.append(rec)

    return records


# ============================================================
# Shard writing
# ============================================================

def _write_one_shard(records, out_dir, shard_idx):
    """Write one shard as zstd-compressed JSONL."""
    import zstandard as zstd

    os.makedirs(out_dir, exist_ok=True)
    shard_name = f'shard_{shard_idx:04d}.jsonl.zst'
    shard_path = os.path.join(out_dir, shard_name)

    rng = np.random.RandomState(shard_idx)
    indices = list(range(len(records)))
    rng.shuffle(indices)

    cctx = zstd.ZstdCompressor(level=3)
    lines = [json.dumps(records[i], separators=(',', ':')) for i in indices]
    data = '\n'.join(lines).encode('utf-8')
    compressed = cctx.compress(data)

    with open(shard_path, 'wb') as f:
        f.write(compressed)

    return len(compressed)


def _write_partial_manifest(out_dir, shard_files, n_pairs, n_families,
                            n_processed):
    """Write manifest (partial or final)."""
    manifest = {
        'n_pairs': n_pairs,
        'n_families': n_families,
        'n_shards': len(shard_files),
        'shard_files': list(shard_files),
        'complete': False,
        'families_processed': n_processed,
    }
    path = os.path.join(out_dir, 'manifest.json')
    with open(path, 'w') as f:
        json.dump(manifest, f, indent=2)


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='Pre-compile PANTHER cherry pairs into training shards')
    parser.add_argument('--msa-dir', default=None,
                        help='Directory with *.AN.fasta files (auto-detected)')
    parser.add_argument('--out', default=None,
                        help='Output directory for shards')
    parser.add_argument('--families', default=None,
                        help='Comma-separated family IDs (default: all)')
    parser.add_argument('--split', default=None, choices=['train', 'val', 'test'],
                        help='Use only families in this split')
    parser.add_argument('--split-file', default=None,
                        help='Path to split JSON file')
    parser.add_argument('--max-seqs-per-family', type=int, default=500,
                        help='Max sequences per family (default: 500)')
    parser.add_argument('--pairs-per-shard', type=int, default=3000,
                        help='Pairs per shard file (default: 3000)')
    parser.add_argument('--workers', type=int, default=1,
                        help='Parallel workers (default: 1)')
    parser.add_argument('--max-families', type=int, default=0,
                        help='Process at most N families (0=all)')
    args = parser.parse_args()

    # Auto-detect directories
    repo_root = Path(__file__).resolve().parent.parent.parent
    data_dir = repo_root / "data"

    if args.msa_dir:
        msa_dir = args.msa_dir
    else:
        msa_dir = str(data_dir / "panther" / "PANTHER19.0_data" / "Tree_MSF")

    if args.out:
        out_dir = args.out
    else:
        out_dir = str(data_dir / "panther" / "precompiled")

    # Find MSA files
    import glob as glob_mod
    if args.families:
        family_list = [f.strip() for f in args.families.split(',')]
        msa_files = []
        for fam in family_list:
            p = os.path.join(msa_dir, f"{fam}.AN.fasta")
            if os.path.exists(p):
                msa_files.append(p)
            else:
                _log(f"Warning: {fam} not found in {msa_dir}")
    else:
        msa_files = sorted(glob_mod.glob(os.path.join(msa_dir, '*.AN.fasta')))

    # Apply split filter
    if args.split:
        split_file = args.split_file
        if not split_file:
            candidates = [
                str(repo_root / "fetch" / "panther" / "splits" / "811.json"),
            ]
            for c in candidates:
                if os.path.exists(c):
                    split_file = c
                    break

        if split_file and os.path.exists(split_file):
            with open(split_file) as f:
                split_data = json.load(f)
            split_fams = set(split_data.get(args.split, []))
            msa_files = [f for f in msa_files
                         if os.path.basename(f).replace('.AN.fasta', '') in split_fams]
            _log(f"Split '{args.split}': {len(msa_files)} families")
        else:
            _log(f"ERROR: --split requires a valid split file")
            sys.exit(1)

    if args.max_families > 0:
        msa_files = msa_files[:args.max_families]

    if not msa_files:
        _log("No MSA files found")
        sys.exit(1)

    _log(f"Processing {len(msa_files)} families -> {out_dir}")

    # Process families
    all_records = []
    families_processed = 0
    families_with_pairs = 0
    shard_files = []
    shard_idx = 0
    total_pairs_written = 0

    if args.workers > 1:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {
                pool.submit(_process_family, f, args.max_seqs_per_family): f
                for f in msa_files
            }
            for future in as_completed(futures):
                msa_file = futures[future]
                fam = os.path.basename(msa_file).replace('.AN.fasta', '')
                try:
                    records = future.result()
                    if records:
                        all_records.extend(records)
                        families_with_pairs += 1
                    families_processed += 1

                    # Flush shard
                    while len(all_records) >= args.pairs_per_shard:
                        shard_batch = all_records[:args.pairs_per_shard]
                        all_records = all_records[args.pairs_per_shard:]
                        _write_one_shard(shard_batch, out_dir, shard_idx)
                        shard_files.append(f'shard_{shard_idx:04d}.jsonl.zst')
                        shard_idx += 1
                        total_pairs_written += len(shard_batch)

                    if families_processed % max(1, len(msa_files) // 20) == 0:
                        _log(f"  {families_processed}/{len(msa_files)} families, "
                             f"{total_pairs_written + len(all_records)} pairs")
                except Exception as e:
                    _log(f"  ERROR {fam}: {e}")
                    families_processed += 1
    else:
        for msa_file in msa_files:
            fam = os.path.basename(msa_file).replace('.AN.fasta', '')
            try:
                records = _process_family(msa_file, args.max_seqs_per_family)
                if records:
                    all_records.extend(records)
                    families_with_pairs += 1
                families_processed += 1

                while len(all_records) >= args.pairs_per_shard:
                    shard_batch = all_records[:args.pairs_per_shard]
                    all_records = all_records[args.pairs_per_shard:]
                    _write_one_shard(shard_batch, out_dir, shard_idx)
                    shard_files.append(f'shard_{shard_idx:04d}.jsonl.zst')
                    shard_idx += 1
                    total_pairs_written += len(shard_batch)
                    _write_partial_manifest(
                        out_dir, shard_files, total_pairs_written,
                        families_with_pairs, families_processed)

                if families_processed % max(1, len(msa_files) // 20) == 0:
                    _log(f"  {families_processed}/{len(msa_files)} families, "
                         f"{total_pairs_written + len(all_records)} pairs, "
                         f"{len(shard_files)} shards")
            except Exception as e:
                _log(f"  ERROR {fam}: {e}")
                families_processed += 1

    # Flush remaining
    if all_records:
        _write_one_shard(all_records, out_dir, shard_idx)
        shard_files.append(f'shard_{shard_idx:04d}.jsonl.zst')
        total_pairs_written += len(all_records)

    _log(f"Total: {total_pairs_written} pairs from "
         f"{families_with_pairs}/{families_processed} families, "
         f"{len(shard_files)} shards")

    if total_pairs_written == 0:
        _log("No pairs found")
        sys.exit(1)

    # Final manifest
    manifest = {
        'n_pairs': total_pairs_written,
        'n_families': families_with_pairs,
        'n_shards': len(shard_files),
        'shard_files': shard_files,
        'pairs_per_shard': args.pairs_per_shard,
        'complete': True,
        'source': {
            'dataset': 'PANTHER 19.0',
            'msa_dir': os.path.abspath(msa_dir),
            'split': args.split,
            'max_seqs_per_family': args.max_seqs_per_family,
        },
    }
    manifest_path = os.path.join(out_dir, 'manifest.json')
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    _log(f"Done: {total_pairs_written} pairs in {len(shard_files)} shards")
    _log(f"Manifest: {manifest_path}")


if __name__ == '__main__':
    main()
