#!/usr/bin/env python3
"""
Prepare Pfam FULL alignments into pairwise training shards.

Streams Pfam-A.full.gz, parses Stockholm blocks, selects diverse pairs
via composition-based clustering, and writes zstd-compressed JSONL shards
compatible with tkf-mixdom's precompile_pairs.py format.

Output structure:
    data/pfam/full/precompiled/
        manifest.json
        shard_0000.jsonl.zst
        shard_0001.jsonl.zst
        ...

Usage:
    python fetch/pfam/prepare_full.py
    python fetch/pfam/prepare_full.py --max-families 100  # test on subset
"""

import argparse
import gzip
import hashlib
import json
import math
import os
import re
import sys
import time
from pathlib import Path

import numpy as np

_log_start = time.monotonic()


def _log(msg, end='\n'):
    elapsed = time.monotonic() - _log_start
    sys.stderr.write(f"[{elapsed:8.1f}s] {msg}" + end)
    sys.stderr.flush()


# Amino acid alphabet (must match tkfmixdom)
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AA = len(AMINO_ACIDS)
AA_TO_IDX = {a: i for i, a in enumerate(AMINO_ACIDS)}
IDX_TO_AA = {i: a for i, a in enumerate(AMINO_ACIDS)}

# State codes
M_STATE, I_STATE, D_STATE = 1, 2, 3
STATE_TO_CHAR = {M_STATE: 'M', I_STATE: 'I', D_STATE: 'D'}


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
            parts.append(f'{STATE_TO_CHAR[current]}{count}')
            current = s
            count = 1
    parts.append(f'{STATE_TO_CHAR[current]}{count}')
    return ''.join(parts)


def aa_composition(seq):
    """20-dim amino acid composition vector (ignoring gaps)."""
    comp = np.zeros(AA, dtype=np.float32)
    for c in seq:
        idx = AA_TO_IDX.get(c)
        if idx is not None:
            comp[idx] += 1
    total = comp.sum()
    if total > 0:
        comp /= total
    return comp


def p_distance(seq1, seq2):
    """Fraction of aligned, non-gap positions that differ."""
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
    """JC-like p-distance to evolutionary time, clamped to [0.01, 10.0]."""
    if pdist >= 0.95:
        return 5.0
    corrected = 1.0 - pdist * (AA / (AA - 1.0))
    if corrected <= 0.01:
        return 5.0
    t = -math.log(corrected)
    return max(0.01, min(10.0, t))


def aligned_pair_to_record(seq1, seq2, family_id, name1, name2):
    """Extract a pairwise alignment record from two aligned MSA rows.

    Returns a dict in X/A/Y format, or None if pair is degenerate.
    """
    states = []
    anc_chars = []
    desc_chars = []

    for c1, c2 in zip(seq1, seq2):
        has1 = c1 in AA_TO_IDX
        has2 = c2 in AA_TO_IDX
        if has1 and has2:
            states.append(M_STATE)
            anc_chars.append(AA_TO_IDX[c1])
            desc_chars.append(AA_TO_IDX[c2])
        elif has2 and not has1:
            states.append(I_STATE)
            desc_chars.append(AA_TO_IDX[c2])
        elif has1 and not has2:
            states.append(D_STATE)
            anc_chars.append(AA_TO_IDX[c1])
        # both gaps: skip (insert columns in MSA)

    if not states:
        return None

    x = ''.join(IDX_TO_AA[c] for c in anc_chars)
    y = ''.join(IDX_TO_AA[c] for c in desc_chars)
    a = rle_encode(states)

    pd = p_distance(seq1, seq2)
    t_est = pdist_to_evo_time(pd)

    id_str = f'{family_id}:{name1}:{name2}'
    pair_id = hashlib.sha256(id_str.encode()).hexdigest()[:16]

    return {
        'x': x,
        'a': a,
        'y': y,
        't': round(t_est, 6),
        'fam': family_id,
        'id': pair_id,
    }


def simple_kmeans(comps, k, rng, max_iter=10):
    """Simple k-means on composition vectors. Returns labels array."""
    n = len(comps)
    if k >= n:
        return np.arange(n)

    # Initialize centroids with k-means++
    idx = [rng.randint(0, n)]
    for _ in range(k - 1):
        dists = np.min([np.sum((comps - comps[c]) ** 2, axis=1) for c in idx], axis=0)
        dists /= dists.sum() + 1e-30
        idx.append(rng.choice(n, p=dists))
    centroids = comps[idx].copy()

    labels = np.zeros(n, dtype=int)
    for _ in range(max_iter):
        # Assign
        dists_to_centroids = np.array([
            np.sum((comps - centroids[c]) ** 2, axis=1) for c in range(k)
        ])  # (k, n)
        labels = np.argmin(dists_to_centroids, axis=0)
        # Update
        new_centroids = np.zeros_like(centroids)
        for c in range(k):
            mask = labels == c
            if mask.any():
                new_centroids[c] = comps[mask].mean(axis=0)
            else:
                new_centroids[c] = centroids[c]
        if np.allclose(centroids, new_centroids, atol=1e-6):
            break
        centroids = new_centroids

    return labels


def select_pairs_clustered(names, seqs, rng, max_seqs=500, max_pairs=200):
    """Select diverse pairs using composition-based clustering.

    Returns list of (idx1, idx2) pairs.
    """
    n = len(seqs)

    # Subsample if needed
    if n > max_seqs:
        subset = sorted(rng.choice(n, max_seqs, replace=False))
        names_sub = [names[i] for i in subset]
        seqs_sub = [seqs[i] for i in subset]
        idx_map = list(subset)
    else:
        names_sub = names
        seqs_sub = seqs
        idx_map = list(range(n))

    n_sub = len(seqs_sub)
    target_pairs = min(n_sub * 3, max_pairs)

    if n_sub < 4:
        # Just pair them all
        pairs = []
        for i in range(n_sub):
            for j in range(i + 1, n_sub):
                pairs.append((idx_map[i], idx_map[j]))
                if len(pairs) >= target_pairs:
                    return pairs
        return pairs

    # Compute compositions
    comps = np.array([aa_composition(s) for s in seqs_sub])

    # K-means clustering
    k = min(max(n_sub // 5, 2), 20)
    labels = simple_kmeans(comps, k, rng)

    # Build cluster index
    clusters = {}
    for i, lab in enumerate(labels):
        clusters.setdefault(int(lab), []).append(i)

    pairs_set = set()
    pairs = []

    def add_pair(i, j):
        key = (min(i, j), max(i, j))
        if key not in pairs_set:
            pairs_set.add(key)
            pairs.append((idx_map[i], idx_map[j]))

    # ~50% within-cluster pairs
    n_within = target_pairs // 2
    within_count = 0
    cluster_ids = list(clusters.keys())
    attempt = 0
    while within_count < n_within and attempt < n_within * 10:
        c = cluster_ids[rng.randint(0, len(cluster_ids))]
        members = clusters[c]
        if len(members) < 2:
            attempt += 1
            continue
        i, j = rng.choice(len(members), 2, replace=False)
        add_pair(members[i], members[j])
        if len(pairs) > within_count:
            within_count = len(pairs)
        attempt += 1

    # ~50% cross-cluster pairs
    n_cross = target_pairs - len(pairs)
    attempt = 0
    while len(pairs) < target_pairs and attempt < n_cross * 10:
        c1, c2 = rng.choice(len(cluster_ids), 2, replace=False) if len(cluster_ids) >= 2 else (0, 0)
        m1 = clusters[cluster_ids[c1]]
        m2 = clusters[cluster_ids[c2]]
        i = m1[rng.randint(0, len(m1))]
        j = m2[rng.randint(0, len(m2))]
        if i != j:
            add_pair(i, j)
        attempt += 1

    return pairs


def parse_stockholm_block(lines):
    """Parse a single Stockholm block from accumulated lines.

    Returns (accession, names, seqs) or (None, None, None) if invalid.
    """
    accession = None
    seq_data = {}
    name_order = []

    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('#=GF AC'):
            parts = line.split()
            if len(parts) >= 3:
                accession = parts[2].split('.')[0]
        elif line.startswith('#') or line.startswith('//') or not line.strip():
            continue
        else:
            parts = line.split()
            if len(parts) >= 2:
                name, seq = parts[0], parts[1]
                if name in seq_data:
                    seq_data[name] += seq
                else:
                    name_order.append(name)
                    seq_data[name] = seq

    if not accession or len(name_order) < 2:
        return None, None, None

    names = name_order
    seqs = [seq_data[n] for n in names]
    return accession, names, seqs


def process_family(accession, names, seqs, rng, max_seqs=500, max_pairs=200):
    """Process one family: select pairs, encode records.

    Returns list of record dicts.
    """
    n = len(seqs)
    if n < 4:
        return []

    pairs = select_pairs_clustered(names, seqs, rng,
                                   max_seqs=max_seqs, max_pairs=max_pairs)
    records = []
    for i, j in pairs:
        # Forward direction
        rec = aligned_pair_to_record(seqs[i], seqs[j], accession, names[i], names[j])
        if rec:
            records.append(rec)
        # Reverse direction (for reversibility)
        rec_rev = aligned_pair_to_record(seqs[j], seqs[i], accession, names[j], names[i])
        if rec_rev:
            records.append(rec_rev)

    return records


def write_shard(records, out_dir, shard_idx):
    """Write one shard as zstd-compressed JSONL. Returns compressed size."""
    import zstandard as zstd

    shard_name = f'shard_{shard_idx:04d}.jsonl.zst'
    shard_path = os.path.join(out_dir, shard_name)

    # Shuffle within shard
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


def write_manifest(out_dir, shard_files, n_pairs, n_families, complete=False,
                   families_processed=0):
    """Write manifest.json."""
    manifest = {
        'n_pairs': n_pairs,
        'n_families': n_families,
        'n_shards': len(shard_files),
        'shard_files': list(shard_files),
        'complete': complete,
        'families_processed': families_processed,
        'source': 'Pfam-A.full.gz',
        'pairs_per_shard': 3000,
    }
    path = os.path.join(out_dir, 'manifest.json')
    with open(path, 'w') as f:
        json.dump(manifest, f, indent=2)


def main():
    parser = argparse.ArgumentParser(
        description='Prepare Pfam full alignments into pairwise training shards')
    parser.add_argument('--data-dir', type=Path, default=None,
                        help='Path to data/ directory (default: auto-detect)')
    parser.add_argument('--max-families', type=int, default=0,
                        help='Process at most N families (0 = all, for testing)')
    parser.add_argument('--max-seqs-per-family', type=int, default=500,
                        help='Max sequences per family (default: 500)')
    parser.add_argument('--max-pairs-per-family', type=int, default=200,
                        help='Max pairs per family (default: 200)')
    parser.add_argument('--pairs-per-shard', type=int, default=3000,
                        help='Pairs per shard file (default: 3000)')
    args = parser.parse_args()

    # Auto-detect data directory
    if args.data_dir:
        data_dir = args.data_dir
    else:
        data_dir = Path(__file__).resolve().parent.parent.parent / "data"

    full_gz = data_dir / "pfam" / "full" / "Pfam-A.full.gz"
    if not full_gz.exists():
        _log(f"ERROR: {full_gz} not found. Run fetch_full.py first.")
        sys.exit(1)

    out_dir = str(data_dir / "pfam" / "full" / "precompiled")
    os.makedirs(out_dir, exist_ok=True)

    # Load split file for train families
    splits_file = Path(__file__).resolve().parent / "splits" / "811-clan-resistant.json"
    if not splits_file.exists():
        _log(f"ERROR: splits file not found: {splits_file}")
        sys.exit(1)

    with open(splits_file) as f:
        split_data = json.load(f)
    train_families = set(split_data['train'])
    _log(f"Train split: {len(train_families)} families")

    gz_size = full_gz.stat().st_size
    _log(f"Input: {full_gz} ({gz_size / (1024**3):.2f} GB)")
    _log(f"Output: {out_dir}")

    # Streaming parse
    rng = np.random.RandomState(42)
    record_buffer = []
    shard_files = []
    shard_idx = 0
    total_pairs = 0
    families_processed = 0
    families_with_pairs = 0
    families_skipped_not_train = 0
    families_skipped_too_small = 0

    block_lines = []
    in_block = False

    _log("Streaming Pfam-A.full.gz...")

    with gzip.open(str(full_gz), 'rt', encoding='latin-1') as f:
        for line in f:
            if line.startswith('# STOCKHOLM 1.0'):
                block_lines = [line]
                in_block = True
                continue

            if in_block:
                block_lines.append(line)

            if line.strip() == '//' and in_block:
                in_block = False

                # Parse the block
                accession, names, seqs = parse_stockholm_block(block_lines)
                block_lines = []

                if accession is None:
                    continue

                # Filter to train split
                if accession not in train_families:
                    families_skipped_not_train += 1
                    continue

                n_seqs = len(names) if names else 0
                if n_seqs < 4:
                    families_skipped_too_small += 1
                    families_processed += 1
                    continue

                # Process family
                try:
                    records = process_family(
                        accession, names, seqs, rng,
                        max_seqs=args.max_seqs_per_family,
                        max_pairs=args.max_pairs_per_family,
                    )
                except Exception as e:
                    _log(f"  ERROR processing {accession}: {e}")
                    families_processed += 1
                    continue

                if records:
                    record_buffer.extend(records)
                    families_with_pairs += 1

                families_processed += 1

                # Flush shard when buffer is full
                while len(record_buffer) >= args.pairs_per_shard:
                    shard_batch = record_buffer[:args.pairs_per_shard]
                    record_buffer = record_buffer[args.pairs_per_shard:]
                    write_shard(shard_batch, out_dir, shard_idx)
                    shard_files.append(f'shard_{shard_idx:04d}.jsonl.zst')
                    shard_idx += 1
                    total_pairs += len(shard_batch)
                    write_manifest(out_dir, shard_files, total_pairs,
                                   families_with_pairs,
                                   families_processed=families_processed)

                if families_processed % 100 == 0:
                    _log(f"  {families_processed} families processed, "
                         f"{families_with_pairs} with pairs, "
                         f"{total_pairs + len(record_buffer)} total pairs, "
                         f"{len(shard_files)} shards written")

                if args.max_families > 0 and families_processed >= args.max_families:
                    _log(f"  Reached --max-families {args.max_families}, stopping")
                    break

    # Flush remaining records
    if record_buffer:
        write_shard(record_buffer, out_dir, shard_idx)
        shard_files.append(f'shard_{shard_idx:04d}.jsonl.zst')
        total_pairs += len(record_buffer)

    # Final manifest
    write_manifest(out_dir, shard_files, total_pairs, families_with_pairs,
                   complete=True, families_processed=families_processed)

    _log(f"Done.")
    _log(f"  Families processed: {families_processed}")
    _log(f"  Families with pairs: {families_with_pairs}")
    _log(f"  Families skipped (not in train): {families_skipped_not_train}")
    _log(f"  Families skipped (< 4 seqs): {families_skipped_too_small}")
    _log(f"  Total pairs: {total_pairs}")
    _log(f"  Shards: {len(shard_files)}")
    _log(f"  Output: {out_dir}")


if __name__ == "__main__":
    main()
