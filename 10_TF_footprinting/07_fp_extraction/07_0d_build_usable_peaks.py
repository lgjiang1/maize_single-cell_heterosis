"""
07_0d_build_usable_peaks.py — Filter ACRs by coverage per cell type.

A "condition" = one sample (P1 / F1A / F1B / P2) × one peak. A condition
"passes" iff BOTH:
  (a) n_fragments > P95(random) for that sample (detectably above noise)
  (b) n_fragments >= --min-frags-per-condition (default 50 ≈ 100 Tn5 cuts,
      scPrinter's reliable-footprint minimum)

A peak pair is usable iff at least --min-conditions-passing of the 4
forward-direction conditions pass (default 2). Set --require-one-per-side
to additionally require at least one passing condition per coord side.

Set --min-conditions-passing 4 (and no --require-one-per-side) for legacy
"all 4 conditions must pass" strict behaviour.

Reads:
    qc/2_footprinting/random_coverage_summary.tsv.gz   (per-sample P95)
    qc/2_footprinting/acr_coverage.tsv.gz              (per-region n_frags)
    regions/{cross}_peak_pairing.tsv                   (peak_id ↔ regions)
    _logs/04b_jobs.tsv                                 (sample list)

Writes:
    regions/usable_peaks/{cross}/{cell_type}/B73_coord.bed
    regions/usable_peaks/{cross}/{cell_type}/{Ki3,Oh43}_coord.bed
    regions/usable_peaks/{cross}/{cell_type}/pairing.tsv
    regions/usable_peaks/summary.tsv         (per cross × cell type counts)

Usage:
    python 07_0d_build_usable_peaks.py
    python 07_0d_build_usable_peaks.py \
        --only-cell-types C5
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"

CT_PATTERN = re.compile(r"-C(\d+)")


def parse_cell_type(sample_id: str) -> str | None:
    m = CT_PATTERN.search(sample_id)
    return f"C{m.group(1)}" if m else None


def region_key(row: pd.Series) -> str:
    return f"{row['chrom']}:{row['start']}-{row['end']}"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--only-cell-types", default=None,
                   help="Comma-separated cell types to process (default: all). "
                        "Use this when 07_0b was run on a sample subset and "
                        "the full P95 table isn't available yet.")
    p.add_argument("--min-frags-per-condition", type=int, default=50,
                   help="A condition (P1/F1A/F1B/P2 sample × peak) 'passes' "
                        "iff its fragment count meets BOTH this floor "
                        "(default 50 ≈ 100 Tn5 cuts, scPrinter's reliable-"
                        "footprint minimum) AND the per-sample P95(random) "
                        "threshold.")
    p.add_argument("--min-conditions-passing", type=int, default=2,
                   help="A peak pair is usable iff at least this many of "
                        "the 4 conditions (P1, F1A, P2, F1B) pass the "
                        "per-condition floor. Default 2. Use 4 for 'all "
                        "conditions must pass' (legacy strict mode).")
    p.add_argument("--require-one-per-side", action="store_true",
                   help="Additionally require that at least one passing "
                        "condition comes from each coord side (B73-side AND "
                        "Ki3/Oh43-side). Stricter than just N-of-4.")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    only_cts: set[str] | None = None
    if args.only_cell_types:
        only_cts = {x.strip() for x in args.only_cell_types.split(",")}
        print(f"[INFO] restricting to cell types: {sorted(only_cts)}")

    with open(CONFIG_PATH) as f:
        config = yaml.safe_load(f)

    qc_dir = PROJECT_ROOT / "qc" / "2_footprinting"
    regions_dir = PROJECT_ROOT / config["paths"]["regions_dir"]

    rand_sum_path = qc_dir / "random_coverage_summary.tsv.gz"
    acr_cov_path = qc_dir / "acr_coverage.tsv.gz"
    jobs_path = PROJECT_ROOT / "_logs" / "04b_jobs.tsv"
    for p in (rand_sum_path, acr_cov_path, jobs_path):
        if not p.exists():
            print(f"[ERROR] missing: {p}", file=sys.stderr)
            return 1

    print(f"[INFO] loading {rand_sum_path}...")
    rand = pd.read_csv(rand_sum_path, sep="\t", compression="gzip")
    p95_by_sample: dict[str, float] = dict(zip(rand["sample_id"], rand["p95"]))

    print(f"[INFO] loading {acr_cov_path}...")
    cov = pd.read_csv(acr_cov_path, sep="\t", compression="gzip")
    cov["region_str"] = (cov["chrom"].astype(str) + ":"
                          + cov["start"].astype(str) + "-"
                          + cov["end"].astype(str))
    cov["cell_type"] = cov["sample_id"].map(parse_cell_type)
    print(f"[INFO] coverage: {len(cov):,} rows, "
          f"{cov['sample_id'].nunique()} samples, "
          f"{cov['region_str'].nunique():,} unique regions")

    out_root = regions_dir / "usable_peaks"
    out_root.mkdir(parents=True, exist_ok=True)
    summary_rows = []

    for cross_name, cross_cfg in config["crosses"].items():
        p1 = cross_cfg["parents"]["P1"]
        p2 = cross_cfg["parents"]["P2"]

        pairing_path = regions_dir / f"{cross_name}_peak_pairing.tsv"
        if not pairing_path.exists():
            print(f"[ERROR] missing pairing: {pairing_path}", file=sys.stderr)
            return 1
        pairing = pd.read_csv(pairing_path, sep="\t")
        p1_region_col = f"{p1}_region"
        p2_region_col = f"{p2}_region"
        n_pairs_total = len(pairing)
        print(f"\n[INFO] === {cross_name} ({p1}/{p2}): "
              f"{n_pairs_total:,} master pairs ===")

        cross_cov = cov[cov["cross"] == cross_name].copy()

        cts_to_process = config["cell_types"]
        if only_cts is not None:
            cts_to_process = [c for c in cts_to_process if c in only_cts]
            if not cts_to_process:
                print(f"[INFO] {cross_name}: no matching cell types, skipping")
                continue

        for ct in cts_to_process:
            ct_cov = cross_cov[cross_cov["cell_type"] == ct]
            if ct_cov.empty:
                print(f"[WARN] {cross_name} {ct}: no coverage data")
                continue

            # Count passing conditions per region per side (not all-or-none).
            # A peak pair is usable if total passes across both sides
            # >= min_conditions_passing (default 2 of 4).
            n_above_per_side: dict[str, pd.Series] = {}
            sample_counts_per_side: dict[str, int] = {}
            skip_this_ct = False

            for coord in (p1, p2):
                side_cov = ct_cov[ct_cov["coord_system"] == coord]

                # Per (region, sample): is n_frags > p95(sample)?
                side_cov = side_cov.assign(
                    p95=side_cov["sample_id"].map(p95_by_sample)
                )
                missing_mask = side_cov["p95"].isna()
                if missing_mask.any():
                    missing = side_cov.loc[
                        missing_mask, "sample_id"].unique().tolist()
                    print(f"[WARN] {cross_name} {ct} {coord}: missing P95 "
                          f"for {len(missing)} sample(s) ({missing}); "
                          f"skipping them in the filter")
                    side_cov = side_cov.loc[~missing_mask].copy()

                samples_used = side_cov["sample_id"].unique().tolist()
                sample_counts_per_side[coord] = len(samples_used)
                if not samples_used:
                    print(f"[WARN] {cross_name} {ct} {coord}: no samples with "
                          f"P95 available; skipping cell type")
                    skip_this_ct = True
                    break

                # A condition (region × sample) passes iff BOTH:
                #   1. n_fragments > P95(random)  -- detectably above noise
                #   2. n_fragments >= min_frags_per_condition  -- enough cuts
                #      for reliable footprint scoring (~100 Tn5 cuts default)
                side_cov["above_random"] = (
                    side_cov["n_fragments"] > side_cov["p95"])
                side_cov["above_floor"]  = (
                    side_cov["n_fragments"] >= args.min_frags_per_condition)
                side_cov["above"] = side_cov["above_random"] & side_cov["above_floor"]

                # Count passing conditions per region on this side.
                # (Was: groupby.all() -> set of fully-passing regions.)
                n_above_per_side[coord] = (
                    side_cov.groupby("region_str")["above"].sum().astype(int))

            if skip_this_ct:
                continue

            # Pair-level aggregation: total passing conditions across both sides
            pairing_with_use = pairing.copy()
            pairing_with_use[f"n_above_{p1}"] = (
                pairing_with_use[p1_region_col]
                    .map(n_above_per_side[p1]).fillna(0).astype(int))
            pairing_with_use[f"n_above_{p2}"] = (
                pairing_with_use[p2_region_col]
                    .map(n_above_per_side[p2]).fillna(0).astype(int))
            pairing_with_use["n_above_total"] = (
                pairing_with_use[f"n_above_{p1}"]
                + pairing_with_use[f"n_above_{p2}"])

            total_ok = pairing_with_use["n_above_total"] >= args.min_conditions_passing
            if args.require_one_per_side:
                each_side_ok = (
                    (pairing_with_use[f"n_above_{p1}"] >= 1)
                    & (pairing_with_use[f"n_above_{p2}"] >= 1))
                pairing_with_use["pair_ok"] = total_ok & each_side_ok
            else:
                pairing_with_use["pair_ok"] = total_ok

            kept = pairing_with_use[pairing_with_use["pair_ok"]].copy()
            n_kept = len(kept)

            # Per-side usable counts (peaks with >=1 condition passing on that side)
            n_p1_ok = int((pairing_with_use[f"n_above_{p1}"] >= 1).sum())
            n_p2_ok = int((pairing_with_use[f"n_above_{p2}"] >= 1).sum())

            # Write outputs
            ct_dir = out_root / cross_name / ct
            ct_dir.mkdir(parents=True, exist_ok=True)

            for coord, region_col in [(p1, p1_region_col),
                                       (p2, p2_region_col)]:
                bed_rows = kept[region_col].str.extract(
                    r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")
                bed_rows["start"] = bed_rows["start"].astype(int)
                bed_rows["end"] = bed_rows["end"].astype(int)
                bed_path = ct_dir / f"{coord}_coord.bed"
                bed_rows.to_csv(bed_path, sep="\t", header=False, index=False)

            pairing_out = ct_dir / "pairing.tsv"
            kept[["peak_id", p1_region_col, p2_region_col]].to_csv(
                pairing_out, sep="\t", index=False)

            pct = 100.0 * n_kept / n_pairs_total if n_pairs_total else 0.0
            # Distribution of pass-count across all peak pairs (for the log)
            pass_dist = (pairing_with_use["n_above_total"]
                          .value_counts().sort_index())
            n_total_samples = (sample_counts_per_side[p1]
                                + sample_counts_per_side[p2])
            print(f"  {ct}: {n_kept:,}/{n_pairs_total:,} ({pct:.1f}%) usable  "
                  f"(≥{args.min_conditions_passing}/{n_total_samples} conditions"
                  f"{', ≥1 per side' if args.require_one_per_side else ''})")
            print(f"     pass-count distribution: " + ", ".join(
                f"{k}:{v}" for k, v in pass_dist.items()))
            print(f"     {p1}-side ≥1 ok: {n_p1_ok:,}   "
                  f"{p2}-side ≥1 ok: {n_p2_ok:,}   "
                  f"({sample_counts_per_side[p1]}+{sample_counts_per_side[p2]} samples)")

            summary_rows.append({
                "cross": cross_name,
                "cell_type": ct,
                "n_pairs_total": n_pairs_total,
                "n_pairs_usable": n_kept,
                "frac_usable": n_kept / n_pairs_total if n_pairs_total else 0,
                "min_conditions_passing": args.min_conditions_passing,
                "require_one_per_side": args.require_one_per_side,
                f"n_{p1}_side_ge1_ok": n_p1_ok,
                f"n_{p2}_side_ge1_ok": n_p2_ok,
                f"n_samples_{p1}_side": sample_counts_per_side[p1],
                f"n_samples_{p2}_side": sample_counts_per_side[p2],
            })

    sum_path = out_root / "summary.tsv"
    pd.DataFrame(summary_rows).to_csv(sum_path, sep="\t", index=False)
    print(f"\n[OK] {sum_path}")
    print("\n[DONE] usable peak sets written to regions/usable_peaks/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
