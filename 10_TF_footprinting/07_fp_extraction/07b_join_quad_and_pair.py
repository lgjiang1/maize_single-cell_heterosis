#!/usr/bin/env python
"""
Phase 3 Step 07b: Join the 4-condition quad per cell type and pair across
coordinate systems.

For one (cross x cell_type), takes the four per-sample 07a outputs
(P1, F1A in B73 coords; P2, F1B in P2 coords), inner-joins each side on
hit_id, aggregates per (region_str x motif_id), and pairs the two sides
via regions/{cross}_peak_pairing.tsv. Applies the ever-bound filter
(keep peak x signature combinations where fp_present is True in at least
one of the 4 conditions).

Inputs:
  - 4 hit_fp_scores.tsv.gz files (from 07a)
  - regions/{cross}_peak_pairing.tsv (from 04a)
  - config/config.yaml (for parent / f1 names)

Output:
  7_fp_extraction/joined/{cross}/{cell_type}/quad_paired.tsv.gz

Schema (one row per peak x signature, after ever-bound filter):
  peak_id, B73_region, P2_region, motif_id, motif_name,
  n_hits_B73, n_frag_P1, n_frag_F1A,
  max_fp_raw_P1, sum_fp_raw_P1, min_fp_norm_P1, max_z_P1, any_present_P1,
  max_fp_raw_F1A, sum_fp_raw_F1A, min_fp_norm_F1A, max_z_F1A, any_present_F1A,
  n_hits_P2, n_frag_P2, n_frag_F1B,
  max_fp_raw_P2, sum_fp_raw_P2, min_fp_norm_P2, max_z_P2, any_present_P2,
  max_fp_raw_F1B, sum_fp_raw_F1B, min_fp_norm_F1B, max_z_F1B, any_present_F1B,
  ever_bound

Sign convention: 07a now writes -log10(p) z_global (positive = footprint),
so MORE POSITIVE z_global means stronger footprint evidence. We aggregate
with `max` for the "best evidence" stat. We carry both max and sum of
fp_raw because:
  - max captures the strongest single hit's evidence
  - sum captures total binding load (informative when a signature has
    multiple hits per peak)
  - n_hits captures multiplicity (motif gain/loss is itself cis variation)
(min_fp_norm_* is v1-only and applies only to archived z-score inputs.)

Usage:
  python 07b_join_quad_and_pair.py --cross B-K --cell-type C5
"""
from __future__ import annotations

import argparse
import gzip
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cross", required=True, choices=["B-K", "B-O"])
    p.add_argument("--cell-type", required=True,
                   help="C1..C14")
    p.add_argument("--f1-direction", default="forward",
                   choices=["forward", "reverse"],
                   help="forward = BxK/BxO, reverse = KxB/OxB")
    p.add_argument("--fp-dir",
                   default=str(PROJECT_ROOT / "7_fp_extraction"))
    p.add_argument("--peak-pairing-dir",
                   default=str(PROJECT_ROOT / "regions"))
    p.add_argument("--outdir",
                   default=str(PROJECT_ROOT / "7_fp_extraction" / "joined"))
    p.add_argument("--require-ever-bound", action="store_true", default=True,
                   help="Keep only rows with fp_present True in >=1 condition")
    p.add_argument("--keep-all", action="store_true",
                   help="Disable ever-bound filter (overrides --require-ever-bound)")
    p.add_argument("--qc-min-frags", type=int, default=10,
                   help="Per-row qc_pass flag: True if n_frag >= this in "
                        "all 4 conditions. Default 10. Doesn't filter rows; "
                        "08a/08b use it to stratify reports.")
    return p.parse_args()


# ------------------------------------------------------------------
# Loading
# ------------------------------------------------------------------

# Columns we need from each per-sample hit_fp_scores file. Drop hit_center
# (per-hit detail, lost on aggregation). fp_norm is optional (v1 has it,
# v2 drops it because the v2 calibration is coverage-stratified directly).
REQUIRED_COLS = ["hit_id", "region_str", "motif_id", "motif_name",
                 "fp_raw", "z_global", "n_fragments", "fp_present"]
OPTIONAL_COLS = ["fp_norm"]
HIT_DTYPES = {
    "hit_id": np.int64,
    "region_str": "category",
    "motif_id": "category",
    "motif_name": "category",
    "fp_raw": np.float32,
    "fp_norm": np.float32,
    "z_global": np.float32,
    "n_fragments": np.int32,
    "fp_present": "boolean",
}


def detect_input_columns(path: Path) -> tuple[list[str], dict, bool]:
    """Peek at header; return (cols_to_read, dtypes_subset, has_fp_norm)."""
    header = pd.read_csv(path, sep="\t", compression="gzip", nrows=0).columns
    has_fp_norm = "fp_norm" in header
    cols = list(REQUIRED_COLS)
    if has_fp_norm:
        cols.append("fp_norm")
    dtypes = {c: HIT_DTYPES[c] for c in cols if c in HIT_DTYPES}
    return cols, dtypes, has_fp_norm


def load_hits(path: Path, label: str) -> tuple[pd.DataFrame, bool]:
    print(f"[INFO] Loading {label}: {path.name}", flush=True)
    t0 = time.time()
    cols, dtypes, has_fp_norm = detect_input_columns(path)
    df = pd.read_csv(path, sep="\t", compression="gzip",
                     usecols=cols, dtype=dtypes)
    df["region_str"] = df["region_str"].astype(str)
    df["motif_id"] = df["motif_id"].astype(str)
    print(f"[INFO]   loaded {len(df):,} rows in {time.time()-t0:.1f}s "
          f"(fp_norm: {'present' if has_fp_norm else 'absent (v2)'})",
          flush=True)
    return df, has_fp_norm


# ------------------------------------------------------------------
# Within-coord aggregation
# ------------------------------------------------------------------

def aggregate_one_side(df_a: pd.DataFrame, df_b: pd.DataFrame,
                       a_label: str, b_label: str,
                       has_fp_norm: bool) -> pd.DataFrame:
    """Inner-join two same-coord-system hit tables on hit_id, then aggregate
    per (region_str, motif_id).

    a_label / b_label suffix the per-condition columns in the output
    (e.g. 'P1' and 'F1A' on the B73 side).
    has_fp_norm: include min_fp_norm_* aggregates (v1 only).
    """
    print(f"[INFO] Inner-joining {a_label} + {b_label} on hit_id", flush=True)
    t0 = time.time()
    if not df_a["hit_id"].equals(df_b["hit_id"]):
        df_a = df_a.sort_values("hit_id").reset_index(drop=True)
        df_b = df_b.sort_values("hit_id").reset_index(drop=True)
        if not df_a["hit_id"].equals(df_b["hit_id"]):
            raise RuntimeError(
                f"{a_label} and {b_label} hit_id columns differ — "
                "did 07a use different motif hits between the two samples?"
            )

    merged = df_a[["region_str", "motif_id", "motif_name"]].copy()
    for side, df_side in [(a_label, df_a), (b_label, df_b)]:
        merged[f"fp_raw_{side}"]     = df_side["fp_raw"].values
        merged[f"z_{side}"]          = df_side["z_global"].values
        merged[f"n_frag_{side}"]     = df_side["n_fragments"].values
        merged[f"fp_present_{side}"] = df_side["fp_present"].values
        if has_fp_norm:
            merged[f"fp_norm_{side}"] = df_side["fp_norm"].values
    print(f"[INFO]   merged in {time.time()-t0:.1f}s ({len(merged):,} rows)",
          flush=True)

    print(f"[INFO] Aggregating per (region_str x motif_id) for {a_label}/{b_label}",
          flush=True)
    t1 = time.time()
    agg_spec = {
        f"n_frag_{a_label}":     (f"n_frag_{a_label}", "first"),
        f"n_frag_{b_label}":     (f"n_frag_{b_label}", "first"),
        f"max_fp_raw_{a_label}": (f"fp_raw_{a_label}", "max"),
        f"sum_fp_raw_{a_label}": (f"fp_raw_{a_label}", "sum"),
        f"max_z_{a_label}":      (f"z_{a_label}", "max"),
        f"any_present_{a_label}":(f"fp_present_{a_label}", "any"),
        f"max_fp_raw_{b_label}": (f"fp_raw_{b_label}", "max"),
        f"sum_fp_raw_{b_label}": (f"fp_raw_{b_label}", "sum"),
        f"max_z_{b_label}":      (f"z_{b_label}", "max"),
        f"any_present_{b_label}":(f"fp_present_{b_label}", "any"),
    }
    if has_fp_norm:
        agg_spec[f"min_fp_norm_{a_label}"] = (f"fp_norm_{a_label}", "min")
        agg_spec[f"min_fp_norm_{b_label}"] = (f"fp_norm_{b_label}", "min")
    agg = merged.groupby(["region_str", "motif_id"], sort=False).agg(
        motif_name=("motif_name", "first"),
        n_hits=("region_str", "size"),
        **agg_spec,
    ).reset_index()
    print(f"[INFO]   aggregated in {time.time()-t1:.1f}s "
          f"({len(agg):,} (region x motif) rows)", flush=True)
    return agg


# ------------------------------------------------------------------
# Cross-coord pairing
# ------------------------------------------------------------------

def pair_via_peak_map(agg_b73: pd.DataFrame, agg_p2: pd.DataFrame,
                      peak_map: pd.DataFrame, p2_label: str
                      ) -> pd.DataFrame:
    """Pair B73-side and P2-side aggregates via the peak pairing map.

    peak_map columns: peak_id, B73_region, {P2_label}_region.
    We attach peak_id to each side, then outer-join on (peak_id, motif_id).
    """
    print("[INFO] Pairing sides via peak_map", flush=True)
    t0 = time.time()
    # Attach peak_id to each side
    b73 = agg_b73.merge(
        peak_map[["peak_id", "B73_region"]],
        left_on="region_str", right_on="B73_region", how="inner",
    ).drop(columns=["region_str"])

    p2_region_col = f"{p2_label}_region"
    p2 = agg_p2.merge(
        peak_map[["peak_id", p2_region_col]],
        left_on="region_str", right_on=p2_region_col, how="inner",
    ).drop(columns=["region_str"])

    # Rename n_hits per side for clarity, and same with motif_name (will be
    # identical across sides; coalesce later).
    b73 = b73.rename(columns={"n_hits": "n_hits_B73",
                              "motif_name": "motif_name_B73"})
    p2 = p2.rename(columns={"n_hits": f"n_hits_{p2_label}",
                            "motif_name": "motif_name_P2"})

    paired = b73.merge(p2, on=["peak_id", "motif_id"], how="outer")
    # Coalesce motif_name
    paired["motif_name"] = paired["motif_name_B73"].fillna(
        paired["motif_name_P2"])
    paired = paired.drop(columns=["motif_name_B73", "motif_name_P2"])

    # Fill regions from peak_map for rows missing one side
    region_lookup = peak_map.set_index("peak_id")
    paired["B73_region"] = paired["B73_region"].fillna(
        paired["peak_id"].map(region_lookup["B73_region"]))
    paired[p2_region_col] = paired[p2_region_col].fillna(
        paired["peak_id"].map(region_lookup[p2_region_col]))

    # Hit counts: NaN -> 0 (motif absent on that side)
    paired["n_hits_B73"] = paired["n_hits_B73"].fillna(0).astype(int)
    paired[f"n_hits_{p2_label}"] = (
        paired[f"n_hits_{p2_label}"].fillna(0).astype(int)
    )

    # Boolean any_present_* columns: NaN -> False (no hits -> not present)
    for col in [c for c in paired.columns if c.startswith("any_present_")]:
        paired[col] = paired[col].fillna(False).astype(bool)

    print(f"[INFO]   paired in {time.time()-t0:.1f}s "
          f"({len(paired):,} (peak x motif) rows)", flush=True)
    return paired


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

def main() -> int:
    args = parse_args()

    with open(CONFIG_PATH) as f:
        config = yaml.safe_load(f)
    cross_cfg = config["crosses"][args.cross]
    p1_name = cross_cfg["parents"]["P1"]            # B73
    p2_name = cross_cfg["parents"]["P2"]            # Ki3 or Oh43
    f1_name = cross_cfg["f1_directions"][args.f1_direction]  # BxK / KxB / ...

    fp_dir = Path(args.fp_dir)
    sample_dir = lambda coord, sample: (
        fp_dir / args.cross / f"{coord}_coords"
        / f"{args.cross}_{coord}_{sample}" / "hit_fp_scores.tsv.gz"
    )
    paths = {
        "P1":  sample_dir(p1_name, f"{p1_name}-{args.cell_type}"),
        "F1A": sample_dir(p1_name,
                          f"{f1_name}-{args.cell_type}.{p1_name}allele"),
        "P2":  sample_dir(p2_name, f"{p2_name}-{args.cell_type}"),
        "F1B": sample_dir(p2_name,
                          f"{f1_name}-{args.cell_type}.{p2_name}allele"),
    }
    for k, v in paths.items():
        if not v.exists():
            print(f"[ERROR] {k} input missing: {v}", file=sys.stderr)
            return 1
        print(f"[INFO] {k}: {v}")

    peak_pairing_path = (Path(args.peak_pairing_dir)
                         / f"{args.cross}_peak_pairing.tsv")
    if not peak_pairing_path.exists():
        print(f"[ERROR] peak pairing missing: {peak_pairing_path}",
              file=sys.stderr)
        return 1
    print(f"[INFO] peak pairing: {peak_pairing_path}")
    peak_map = pd.read_csv(peak_pairing_path, sep="\t")
    print(f"[INFO]   {len(peak_map):,} paired peaks")

    # The Ki3 / Oh43 column in peak_map is named {P2}_region
    p2_region_col = f"{p2_name}_region"
    if p2_region_col not in peak_map.columns:
        print(f"[ERROR] peak_map missing column {p2_region_col}",
              file=sys.stderr)
        return 1

    # ----- B73 side -----
    df_P1, has_fpn_p1 = load_hits(paths["P1"], "P1")
    df_F1A, has_fpn_f1a = load_hits(paths["F1A"], "F1A")
    has_fp_norm = has_fpn_p1 and has_fpn_f1a   # both sides must agree
    agg_b73 = aggregate_one_side(df_P1, df_F1A, "P1", "F1A", has_fp_norm)
    del df_P1, df_F1A

    # ----- P2 side -----
    df_P2, _ = load_hits(paths["P2"], "P2")
    df_F1B, _ = load_hits(paths["F1B"], "F1B")
    agg_p2 = aggregate_one_side(df_P2, df_F1B, "P2", "F1B", has_fp_norm)
    del df_P2, df_F1B

    # ----- Pair via peak map -----
    paired = pair_via_peak_map(agg_b73, agg_p2, peak_map, p2_name)
    del agg_b73, agg_p2

    # ----- Ever-bound flag + filter -----
    present_cols = ["any_present_P1", "any_present_F1A",
                    "any_present_P2", "any_present_F1B"]
    paired["ever_bound"] = paired[present_cols].any(axis=1)
    n_ever_bound = int(paired["ever_bound"].sum())
    print(f"[INFO] ever_bound = True in {n_ever_bound:,} / {len(paired):,} "
          f"rows ({100.0*n_ever_bound/len(paired):.2f}%)")

    # ----- qc_pass flag: all 4 conditions have at least --qc-min-frags --
    qc_thr = args.qc_min_frags
    n_cols = ["n_frag_P1", "n_frag_P2", "n_frag_F1A", "n_frag_F1B"]
    paired["qc_pass"] = (
        (paired["n_frag_P1"].fillna(0)  >= qc_thr) &
        (paired["n_frag_P2"].fillna(0)  >= qc_thr) &
        (paired["n_frag_F1A"].fillna(0) >= qc_thr) &
        (paired["n_frag_F1B"].fillna(0) >= qc_thr)
    )
    n_qc_pass = int(paired["qc_pass"].sum())
    print(f"[INFO] qc_pass (all 4 conditions >= {qc_thr} frags): "
          f"{n_qc_pass:,} / {len(paired):,} "
          f"({100.0*n_qc_pass/len(paired):.2f}%)")

    keep_all = args.keep_all
    if not keep_all:
        paired = paired[paired["ever_bound"]].copy()
        print(f"[INFO] After ever-bound filter: {len(paired):,} rows")

    # ----- Column ordering -----
    p2 = p2_name  # for column naming
    cols_meta = ["peak_id", "B73_region", p2_region_col,
                 "motif_id", "motif_name"]
    cols_b73 = ["n_hits_B73", "n_frag_P1", "n_frag_F1A",
                "max_fp_raw_P1", "sum_fp_raw_P1",
                "max_z_P1", "any_present_P1",
                "max_fp_raw_F1A", "sum_fp_raw_F1A",
                "max_z_F1A", "any_present_F1A"]
    cols_p2 = [f"n_hits_{p2}", "n_frag_P2", "n_frag_F1B",
               "max_fp_raw_P2", "sum_fp_raw_P2",
               "max_z_P2", "any_present_P2",
               "max_fp_raw_F1B", "sum_fp_raw_F1B",
               "max_z_F1B", "any_present_F1B"]
    # If v1 inputs had fp_norm, insert those cols after the matching fp_raw.
    # NOTE: min_fp_norm_* is v1-only; under v1 the columns were min_-prefixed
    # because z-scores had "more negative = stronger footprint". Keep the v1
    # naming on archived data; for the new pipeline these columns aren't used.
    if has_fp_norm:
        cols_b73 = ["n_hits_B73", "n_frag_P1", "n_frag_F1A",
                    "max_fp_raw_P1", "sum_fp_raw_P1", "min_fp_norm_P1",
                    "max_z_P1", "any_present_P1",
                    "max_fp_raw_F1A", "sum_fp_raw_F1A", "min_fp_norm_F1A",
                    "max_z_F1A", "any_present_F1A"]
        cols_p2 = [f"n_hits_{p2}", "n_frag_P2", "n_frag_F1B",
                   "max_fp_raw_P2", "sum_fp_raw_P2", "min_fp_norm_P2",
                   "max_z_P2", "any_present_P2",
                   "max_fp_raw_F1B", "sum_fp_raw_F1B", "min_fp_norm_F1B",
                   "max_z_F1B", "any_present_F1B"]
    cols_final = cols_meta + cols_b73 + cols_p2 + ["ever_bound", "qc_pass"]
    missing = [c for c in cols_final if c not in paired.columns]
    if missing:
        print(f"[WARN] missing expected columns: {missing}")
        cols_final = [c for c in cols_final if c in paired.columns]
    paired = paired[cols_final]

    # ----- Write output -----
    out_dir = (Path(args.outdir) / args.cross / args.cell_type)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "quad_paired.tsv.gz"
    print(f"[INFO] Writing {out_path}")
    t0 = time.time()
    paired.to_csv(out_path, sep="\t", index=False, compression="gzip",
                  float_format="%.6g")
    print(f"[INFO]   wrote in {time.time()-t0:.1f}s")

    # ----- Summary -----
    print()
    print("=== SUMMARY ===")
    print(f"  cross:       {args.cross}")
    print(f"  cell_type:   {args.cell_type}")
    print(f"  F1:          {f1_name}")
    print(f"  fp version:  {'v1 (has fp_norm)' if has_fp_norm else 'v2 (no fp_norm)'}")
    print(f"  rows out:    {len(paired):,}")
    print(f"  ever_bound:  {n_ever_bound:,} "
          f"({100.0*n_ever_bound/max(1,len(paired)):.1f}% of output)")
    n_qc_in_out = int(paired["qc_pass"].sum())
    print(f"  qc_pass (>={qc_thr} frags all 4): {n_qc_in_out:,} "
          f"({100.0*n_qc_in_out/max(1,len(paired)):.1f}% of output)")

    only_b73 = ((paired["n_hits_B73"] > 0)
                & (paired[f"n_hits_{p2}"] == 0)).sum()
    only_p2 = ((paired["n_hits_B73"] == 0)
               & (paired[f"n_hits_{p2}"] > 0)).sum()
    both = ((paired["n_hits_B73"] > 0)
            & (paired[f"n_hits_{p2}"] > 0)).sum()
    print(f"  motif present only in B73 (cis loss in P2): {only_b73:,}")
    print(f"  motif present only in {p2} (cis loss in B73): {only_p2:,}")
    print(f"  motif present on both sides:                 {both:,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
