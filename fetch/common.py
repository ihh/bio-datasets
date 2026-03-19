"""
Common utilities for bio-datasets fetch scripts.

All fetch scripts follow this contract:
  - --outdir defaults to data/<dataset>/ relative to repo root
  - Idempotent: skip files that already exist (by size check)
  - Never delete symlinks or existing data
  - Log progress to stderr
"""

import logging
import os
import shutil
import sys
import urllib.request
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("fetch")


def repo_root() -> Path:
    """Return the bio-datasets repo root (parent of fetch/)."""
    return Path(__file__).resolve().parent.parent


def default_outdir(dataset: str) -> Path:
    """Default output directory: <repo_root>/data/<dataset>."""
    return repo_root() / "data" / dataset


def safe_outdir(path: Path) -> Path:
    """Create outdir if needed. Refuse to overwrite symlinks."""
    if path.is_symlink():
        target = path.resolve()
        log.info("Output directory %s is a symlink → %s (using target)", path, target)
        return target
    path.mkdir(parents=True, exist_ok=True)
    return path


def file_exists(path: Path, min_size: int = 1) -> bool:
    """Check if file exists with at least min_size bytes."""
    if path.is_symlink():
        # Resolve symlink and check target
        target = path.resolve()
        return target.exists() and target.stat().st_size >= min_size
    return path.exists() and path.stat().st_size >= min_size


def safe_download(url: str, dest: Path, min_size: int = 1) -> bool:
    """Download url → dest. Skip if exists. Never clobber symlinks.

    Returns True if file is present after call (downloaded or pre-existing).
    """
    if dest.is_symlink():
        log.info("SKIP (symlink): %s → %s", dest.name, dest.resolve())
        return dest.resolve().exists()

    if file_exists(dest, min_size):
        log.info("SKIP (exists): %s (%d bytes)", dest.name, dest.stat().st_size)
        return True

    log.info("Downloading: %s → %s", url, dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    try:
        urllib.request.urlretrieve(url, str(tmp))
        shutil.move(str(tmp), str(dest))
        log.info("OK: %s (%d bytes)", dest.name, dest.stat().st_size)
        return True
    except Exception as e:
        log.error("FAILED: %s — %s", url, e)
        tmp.unlink(missing_ok=True)
        return False


def safe_decompress_gz(gz_path: Path, out_path: Path) -> bool:
    """Decompress .gz file. Skip if output exists. Never clobber symlinks."""
    import gzip

    if out_path.is_symlink():
        log.info("SKIP (symlink): %s", out_path.name)
        return out_path.resolve().exists()

    if file_exists(out_path):
        log.info("SKIP (exists): %s", out_path.name)
        return True

    log.info("Decompressing: %s → %s", gz_path.name, out_path.name)
    with gzip.open(gz_path, "rb") as f_in:
        with open(out_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    return True
