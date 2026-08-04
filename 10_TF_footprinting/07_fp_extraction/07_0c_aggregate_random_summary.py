"""
07_0c_aggregate_random_summary.py — Per-sample percentiles of random coverage.

Reads:
    qc/2_footprinting/random_coverage.tsv.gz      (output of 07_0b)
    qc/2_footprinting/acr_coverage_summary.tsv    (output of 04d, optional;
                                                   for ACR vs random comparison)

Writes:
    qc/2_footprinting/random_coverage_summary.tsv.gz
        cols: sample_id, cross, coord_system,
              mean, p50, p75, p90, p95, p99, n_windows

Prints a diagnostic table comparing per-sample P95(random) to median(ACR)
(if 04d summary is available). Ratio < 1 = ACRs not detectably above random.

Usage:
    python 07_0c_aggregate_random_summary.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def main() -> int:
    qc_dir = PROJECT_ROOT / "qc" / "2_footprinting"
    in_path = qc_dir / "random_coverage.tsv.gz"
    if not in_path.exists():
        print(f"[ERROR] {in_path} not found (run 07_0b first)", file=sys.stderr)
        return 1

    print(f"[INFO] loading {in_path}...")
    df = pd.read_csv(in_path, sep="\t", compression="gzip",
                     usecols=["sample_id", "cross", "coord_system",
                              "n_fragments"])
    print(f"[INFO] {len(df):,} rows, {df['sample_id'].nunique()} samples")

    def _pct(q: int):
        def _f(x):
            return float(np.percentile(x, q))
        _f.__name__ = f"p{q}"
        return _f

    summary = (df.groupby(["sample_id", "cross", "coord_system"])
                 ["n_fragments"]
                 .agg([
                     ("mean", "mean"),
                     ("p50", _pct(50)),
                     ("p75", _pct(75)),
                     ("p90", _pct(90)),
                     ("p95", _pct(95)),
                     ("p99", _pct(99)),
                     ("n_windows", "count"),
                 ])
                 .reset_index())

    out_path = qc_dir / "random_coverage_summary.tsv.gz"
    summary.to_csv(out_path, sep="\t", index=False, compression="gzip")
    print(f"[OK] {out_path} ({len(summary)} rows)")

    # Print summary stats
    print("\n[INFO] P95(random) distribution across samples:")
    print(summary["p95"].describe().to_string())
    print("\n[INFO] P95(random) by coord_system:")
    print(summary.groupby("coord_system")["p95"].describe().to_string())

    # Optional: compare to ACR coverage if 04d summary is available
    acr_sum_path = qc_dir / "acr_coverage_summary.tsv"
    if acr_sum_path.exists():
        acr_sum = pd.read_csv(acr_sum_path, sep="\t")
        merged = summary.merge(
            acr_sum[["sample_id", "median_coverage", "frac_zero_coverage"]],
            on="sample_id", how="left",
        )
        merged["acr_median_over_p95"] = (
            merged["median_coverage"] / merged["p95"].replace(0, np.nan)
        )
        print("\n[INFO] Sample-of-samples sanity check (per-coord medians):")
        print(merged.groupby("coord_system").agg(
            p95_random_median=("p95", "median"),
            acr_median_median=("median_coverage", "median"),
            ratio_median=("acr_median_over_p95", "median"),
        ).round(2).to_string())

        # Flag any sample where ACR median <= P95(random)
        flagged = merged[merged["acr_median_over_p95"] <= 1].copy()
        if len(flagged):
            print(f"\n[WARN] {len(flagged)} samples where median ACR coverage "
                  f"≤ P95(random) — these will have very few usable peaks:")
            print(flagged[["sample_id", "p95", "median_coverage",
                           "acr_median_over_p95"]].to_string(index=False))
        else:
            print("\n[OK] All samples: median ACR coverage > P95(random)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
