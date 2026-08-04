"""
04d_count_acr_coverage.py — Count fragments per ACR per condition.

For each of 168 (fragment_file, region_bed) pairs, counts the number of
fragments whose midpoint falls within each region. Outputs a long-format
coverage table and a per-condition summary.

These counts are the `n` in the depth correction model:
    δ_hat = z / √n       (depth-corrected effect size)
    SE    = √(1/n_A + 1/n_B)  (standard error of a contrast)

Reads:
    _logs/04b_jobs.tsv                      — job table (fragment paths + region BEDs)
    3_fragments/{...}/*.frags.tsv.gz        — fragment files
    regions/{...}/regions_2000bp.bed        — region BEDs

Writes:
    qc/2_footprinting/acr_coverage.tsv.gz       — long format: (region × condition)
    qc/2_footprinting/acr_coverage_summary.tsv  — per-condition summary

Usage (from project root):
    python 04d_count_acr_coverage.py
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"


def load_regions(bed_path: str) -> pd.DataFrame:
    """Load a 3-column 0-based BED into a DataFrame with an interval index."""
    df = pd.read_csv(bed_path, sep="\t", header=None, names=["chrom", "start", "end"])
    df["region_idx"] = range(len(df))
    return df


def count_overlaps(frag_path: str, regions: pd.DataFrame) -> np.ndarray:
    """Count fragment midpoints falling within each region.

    Uses a simple per-chromosome sorted-array approach — efficient for
    ~100K fragments and ~65K regions on 10 chromosomes.
    """
    n_regions = len(regions)
    counts = np.zeros(n_regions, dtype=np.int64)

    # Build per-chrom region lookup: sorted arrays of (start, end, idx)
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
            # Binary search for candidate regions
            lo = np.searchsorted(ends, midpoint, side="right")
            hi = np.searchsorted(starts, midpoint, side="right")
            for j in range(lo, hi):
                if starts[j] <= midpoint < ends[j]:
                    counts[idxs[j]] += 1

    return counts


def main() -> int:
    with open(CONFIG_PATH) as f:
        config = yaml.safe_load(f)

    jobs_path = PROJECT_ROOT / "_logs" / "04b_jobs.tsv"
    if not jobs_path.exists():
        print(f"[ERROR] {jobs_path} not found", file=sys.stderr)
        return 1

    jobs = pd.read_csv(jobs_path, sep="\t")
    print(f"[INFO] loaded {len(jobs)} jobs from {jobs_path}")

    qc_dir = PROJECT_ROOT / "qc" / "2_footprinting"
    qc_dir.mkdir(parents=True, exist_ok=True)

    # Group jobs by region BED (each coord system shares one BED)
    all_coverage_rows = []
    summary_rows = []

    # Cache loaded regions
    region_cache: dict[str, pd.DataFrame] = {}

    for i, job in jobs.iterrows():
        sample_id = job["sample_id"]
        frag_path = job["frag_out"]
        region_bed = job["region_bed"]

        if not Path(frag_path).exists():
            print(f"[WARN] missing fragments: {frag_path}")
            summary_rows.append({
                "sample_id": sample_id,
                "cross": job["cross"],
                "coord_system": job["coord_system"],
                "total_fragments": 0,
                "n_regions": 0,
                "median_coverage": 0,
                "frac_zero_coverage": 1.0,
            })
            continue

        # Load regions (cached per BED file)
        if region_bed not in region_cache:
            region_cache[region_bed] = load_regions(region_bed)
        regions = region_cache[region_bed]

        print(f"[{i+1}/{len(jobs)}] {sample_id}...", end="", flush=True)
        counts = count_overlaps(frag_path, regions)
        total_frags = counts.sum()
        n_regions = len(regions)
        median_cov = float(np.median(counts))
        frac_zero = float((counts == 0).sum()) / n_regions if n_regions > 0 else 1.0

        print(f" total={total_frags:,} median_cov={median_cov:.1f} zero={frac_zero:.3f}")

        # Long-format coverage rows
        for j, (_, r) in enumerate(regions.iterrows()):
            all_coverage_rows.append({
                "cross": job["cross"],
                "coord_system": job["coord_system"],
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
            "coord_system": job["coord_system"],
            "total_fragments": int(total_frags),
            "n_regions": n_regions,
            "median_coverage": median_cov,
            "frac_zero_coverage": frac_zero,
        })

    # Write coverage matrix (long format, gzipped — can be large)
    cov_path = qc_dir / "acr_coverage.tsv.gz"
    print(f"\n[INFO] writing {len(all_coverage_rows):,} rows to {cov_path}...")
    cov_df = pd.DataFrame(all_coverage_rows)
    cov_df.to_csv(cov_path, sep="\t", index=False, compression="gzip")
    print(f"[OK] {cov_path}")

    # Write summary
    sum_path = qc_dir / "acr_coverage_summary.tsv"
    sum_df = pd.DataFrame(summary_rows)
    sum_df.to_csv(sum_path, sep="\t", index=False)
    print(f"[OK] {sum_path} ({len(sum_df)} rows)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
