#!/usr/bin/env python
"""
Phase 3 Step 06b: Merge per-chunk motif scan outputs.

Concatenates chunk_NN/motif_hits.tsv.gz files into a single
motif_hits.tsv.gz at the scan root directory.

Usage:
  # Merge one scan:
  python 06b_merge_scan_chunks.py --scan-dir 5_motif_scanning/hits/B-K/B73_coords

  # Merge all 4 scans:
  python 06b_merge_scan_chunks.py --all
"""

import argparse
import gzip
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
HITS_ROOT = PROJECT / "5_motif_scanning" / "hits"

ALL_SCAN_DIRS = [
    HITS_ROOT / "B-K" / "B73_coords",
    HITS_ROOT / "B-K" / "Ki3_coords",
    HITS_ROOT / "B-O" / "B73_coords",
    HITS_ROOT / "B-O" / "Oh43_coords",
]


def merge_one_scan(scan_dir: Path) -> None:
    """Merge all chunk_*/motif_hits.tsv.gz in scan_dir."""
    chunk_files = sorted(scan_dir.glob("chunk_*/motif_hits.tsv.gz"))
    if not chunk_files:
        print(f"[WARN] No chunk files found in {scan_dir}, skipping")
        return

    frames = []
    for cf in chunk_files:
        df = pd.read_csv(cf, sep="\t", compression="gzip")
        chunk_id = cf.parent.name  # e.g. "chunk_03"
        print(f"  {chunk_id}: {len(df):,} hits")
        frames.append(df)

    merged = pd.concat(frames, ignore_index=True)
    # Reassign global hit_id
    merged["hit_id"] = range(len(merged))

    out_path = scan_dir / "motif_hits.tsv.gz"
    with gzip.open(str(out_path), "wt") as f:
        merged.to_csv(f, sep="\t", index=False)

    print(f"[DONE] {out_path} — {len(merged):,} total hits "
          f"from {len(chunk_files)} chunks")


def main():
    p = argparse.ArgumentParser(description="Merge per-chunk motif scan outputs")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--scan-dir", type=Path, help="Single scan directory to merge")
    g.add_argument("--all", action="store_true",
                   help="Merge all 4 scan directories")
    args = p.parse_args()

    if args.all:
        for sd in ALL_SCAN_DIRS:
            print(f"\n[INFO] Merging {sd.relative_to(PROJECT)}")
            merge_one_scan(sd)
    else:
        print(f"[INFO] Merging {args.scan_dir}")
        merge_one_scan(args.scan_dir)


if __name__ == "__main__":
    main()
