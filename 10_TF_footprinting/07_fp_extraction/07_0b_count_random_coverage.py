"""
07_0b_count_random_coverage.py — Count fragments at random non-ACR windows.

Mirrors 04d_count_acr_coverage.py exactly, but counts against the per-parent
random regions written by 07_0a (qc/2_footprinting/random_regions/) instead
of the real ACR set. Used to derive a per-sample "above-random" coverage
threshold in 07_0c / 07_0d.

Reads:
    _logs/04b_jobs.tsv                          — job table (uses frag_out + coord_system)
    qc/2_footprinting/random_regions/{parent}_random_2kb.bed

Writes:
    qc/2_footprinting/random_coverage.tsv.gz    — long format
    qc/2_footprinting/random_coverage_summary.tsv

CLI:
    --task-ids 9,23,51,65    # restrict to a subset of task IDs (testing)

Usage (from project root):
    python 07_0b_count_random_coverage.py
    python 07_0b_count_random_coverage.py --task-ids 9,23,51,65
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--jobs-tsv",
                   default=str(PROJECT_ROOT / "_logs" / "04b_jobs.tsv"))
    p.add_argument("--random-regions-dir",
                   default=str(PROJECT_ROOT / "qc" / "2_footprinting"
                               / "random_regions"))
    p.add_argument("--outdir",
                   default=str(PROJECT_ROOT / "qc" / "2_footprinting"))
    p.add_argument("--task-ids", default=None,
                   help="Comma-separated subset of task IDs (default: all 168)")
    return p.parse_args()


def load_regions(bed_path: str) -> pd.DataFrame:
    df = pd.read_csv(bed_path, sep="\t", header=None,
                     names=["chrom", "start", "end"])
    df["region_idx"] = range(len(df))
    return df


def count_overlaps(frag_path: str, regions: pd.DataFrame) -> np.ndarray:
    n_regions = len(regions)
    counts = np.zeros(n_regions, dtype=np.int64)

    chrom_regions: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for chrom, grp in regions.groupby("chrom"):
        starts = grp["start"].values
        ends = grp["end"].values
        idxs = grp["region_idx"].values
        order = np.argsort(starts)
        chrom_regions[chrom] = (starts[order], ends[order], idxs[order])

    with gzip.open(frag_path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            frag_start = int(parts[1])
            frag_end = int(parts[2])
            midpoint = (frag_start + frag_end) // 2

            cr = chrom_regions.get(chrom)
            if cr is None:
                continue

            starts, ends, idxs = cr
            lo = np.searchsorted(ends, midpoint, side="right")
            hi = np.searchsorted(starts, midpoint, side="right")
            for j in range(lo, hi):
                if starts[j] <= midpoint < ends[j]:
                    counts[idxs[j]] += 1

    return counts


def main() -> int:
    args = parse_args()

    jobs_path = Path(args.jobs_tsv)
    if not jobs_path.exists():
        print(f"[ERROR] {jobs_path} not found", file=sys.stderr)
        return 1
    jobs = pd.read_csv(jobs_path, sep="\t")
    print(f"[INFO] loaded {len(jobs)} jobs from {jobs_path}")

    if args.task_ids:
        keep = sorted({int(x) for x in args.task_ids.split(",")})
        jobs = jobs.iloc[keep].reset_index(drop=True)
        print(f"[INFO] restricted to task-ids {keep} ({len(jobs)} jobs)")

    random_dir = Path(args.random_regions_dir)
    # Map coord_system -> random regions DataFrame (loaded once)
    random_cache: dict[str, pd.DataFrame] = {}

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_rows = []
    summary_rows = []

    for i, job in jobs.iterrows():
        sample_id = job["sample_id"]
        coord_system = job["coord_system"]
        frag_path = job["frag_out"]

        if coord_system not in random_cache:
            bed = random_dir / f"{coord_system}_random_2kb.bed"
            if not bed.exists():
                print(f"[ERROR] missing random BED for {coord_system}: {bed}",
                      file=sys.stderr)
                return 1
            random_cache[coord_system] = load_regions(str(bed))
            print(f"[INFO] loaded random regions for {coord_system}: "
                  f"{len(random_cache[coord_system]):,}")

        regions = random_cache[coord_system]

        if not Path(frag_path).exists():
            print(f"[WARN] missing fragments: {frag_path}")
            summary_rows.append({
                "sample_id": sample_id,
                "cross": job["cross"],
                "coord_system": coord_system,
                "total_fragments": 0,
                "n_regions": 0,
                "median_coverage": 0,
                "frac_zero_coverage": 1.0,
            })
            continue

        print(f"[{i+1}/{len(jobs)}] {sample_id}...", end="", flush=True)
        counts = count_overlaps(frag_path, regions)
        total_frags = counts.sum()
        n_regions = len(regions)
        median_cov = float(np.median(counts))
        frac_zero = float((counts == 0).sum()) / n_regions if n_regions > 0 else 1.0
        print(f" total={total_frags:,} median_cov={median_cov:.1f} zero={frac_zero:.3f}")

        for j, (_, r) in enumerate(regions.iterrows()):
            all_rows.append({
                "cross": job["cross"],
                "coord_system": coord_system,
                "sample_id": sample_id,
                "region_idx": j,
                "chrom": r["chrom"],
                "start": r["start"],
                "end": r["end"],
                "n_fragments": int(counts[j]),
            })

        summary_rows.append({
            "sample_id": sample_id,
            "cross": job["cross"],
            "coord_system": coord_system,
            "total_fragments": int(total_frags),
            "n_regions": n_regions,
            "median_coverage": median_cov,
            "frac_zero_coverage": frac_zero,
        })

    cov_path = out_dir / "random_coverage.tsv.gz"
    print(f"\n[INFO] writing {len(all_rows):,} rows to {cov_path}...")
    pd.DataFrame(all_rows).to_csv(cov_path, sep="\t", index=False,
                                   compression="gzip")
    print(f"[OK] {cov_path}")

    sum_path = out_dir / "random_coverage_per_sample.tsv"
    pd.DataFrame(summary_rows).to_csv(sum_path, sep="\t", index=False)
    print(f"[OK] {sum_path} ({len(summary_rows)} rows)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
