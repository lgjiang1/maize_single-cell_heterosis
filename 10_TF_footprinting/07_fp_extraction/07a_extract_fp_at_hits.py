#!/usr/bin/env python
"""
Phase 3 Step 07a: Extract FP scores at motif hit positions with a
coverage-stratified empirical null.

Reads scPrinter's footprint output from FP_neglog10p__ALL.h5ad (when
04c was run with `fp_return_pval=true`) or the legacy FP_zscore__ALL.h5ad
(deprecated). The script doesn't care which one — both are signed
per-position-per-scale scores, and the coverage-stratified null treats
them the same.

Architecture:

  Pass 1 -- build coverage-stratified null:
    For every region with n_fragments >= 1, extract the score at the
    bottom-`bg_percentile`% TFBS tile positions (default 20%). Pool all
    (n_fragments, score) pairs. Bin by n_fragments into --n-bins quantile
    bins. Per bin, compute median(score) and a spread statistic (default
    IQR/1.349; MAD or std selectable via --spread). Save the per-bin table
    to null_model.tsv.

  Spread-estimator rationale (IQR over MAD):
    scPrinter's score is sparse at low/mid coverage: in low-n bins, the
    bottom-TFBS tile distribution has 50-80% of values exactly 0. MAD
    collapses to 0 once >50% of values are 0. IQR/1.349 stays positive
    as long as <75% of values are 0, which rescues many more bins.
    Combined with --bg-percentile 20 (twice as many null tiles per region
    as the original 10), this pushes most coverage bins into a
    non-degenerate spread regime.

  Pass 2 -- score hits with coverage-aware z:
    For each motif hit at region with coverage n:
      look up the bin whose [n_min, n_max] contains n
      mu_null = bin's median fp_raw
      sigma_null = bin's MAD
      z_global = (fp_raw - mu_null) / sigma_null
      fp_present = z_global > --presence-z      # -log10(p): higher = footprint

Outputs (one SLURM task per sample):
  --outdir/{cross}/{coord_system}/{sample_id}/hit_fp_scores.tsv.gz
  --outdir/{cross}/{coord_system}/{sample_id}/null_model.tsv

Output columns:
  hit_id, region_str, motif_id, motif_name, hit_center,
  fp_raw, n_fragments,
  bin_idx, bin_n_med, bin_null_med, bin_null_mad,
  z_global, fp_present

Dropped from v1: fp_norm (no longer used; n captured in the binning).

Usage:
  python 07a_v2_extract_fp_at_hits.py --task-id $SLURM_ARRAY_TASK_ID
"""
from __future__ import annotations

import os
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import argparse
import gzip
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]

# TFBS tile grid (verified from test h5ad: shape (1, 180))
CONTEXT_RADIUS = 100
TILE_SIZE = 10
N_TILES = 180  # (2000 - 2*100) / 10
TILE_BP = np.arange(N_TILES) * TILE_SIZE + CONTEXT_RADIUS + TILE_SIZE // 2


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract FP z-scores at motif hits with coverage-"
                    "stratified empirical null (v2)",
    )
    p.add_argument("--task-id", type=int, required=True,
                   help="SLURM array task ID (0-167)")
    p.add_argument("--jobs-tsv",
                   default=str(PROJECT_ROOT / "_logs" / "04b_jobs.tsv"))
    p.add_argument("--hits-dir",
                   default=str(PROJECT_ROOT / "5_motif_scanning" / "hits"))
    p.add_argument("--coverage",
                   default=str(PROJECT_ROOT / "qc" / "2_footprinting"
                               / "acr_coverage.tsv.gz"))
    p.add_argument("--outdir",
                   default=str(PROJECT_ROOT / "7_fp_extraction"),
                   help="Output root")
    p.add_argument("--scales-min", type=int, default=2)
    p.add_argument("--scales-max", type=int, default=20)
    p.add_argument("--hit-window", type=int, default=5)
    p.add_argument("--bg-percentile", type=float, default=20.0,
                   help="Bottom N%% TFBS tiles as unbound positions "
                        "(default 20; was 10 prior to 2026-05-12)")
    p.add_argument("--n-bins", type=int, default=20,
                   help="Number of coverage bins for the null model")
    p.add_argument("--min-bin-tiles", type=int, default=500,
                   help="Min null tiles per bin (merged if smaller)")
    p.add_argument("--spread", choices=("iqr", "mad", "sd"), default="iqr",
                   help="Per-bin spread estimator. iqr=(P75-P25)/1.349 "
                        "(default, non-degenerate up to 75%% zeros); "
                        "mad=1.4826*MAD (legacy, degenerates at 50%% zeros); "
                        "sd=std. Column written to null_model.tsv is always "
                        "'null_mad' for downstream compatibility.")
    p.add_argument("--use-mean", action="store_true",
                   help="Deprecated alias for --spread sd; if set with the "
                        "default --spread iqr, switches spread to sd and "
                        "centre statistic to mean.")
    p.add_argument("--presence-z", type=float, default=2.0,
                   help="z-threshold for presence. With -log10(p) input "
                        "(current default), z_global > threshold = footprint. "
                        "Was -2.0 (z_global < threshold) in the z-score era.")
    p.add_argument("--usable-peaks", default=None,
                   help="Optional 3-col BED of usable regions (chrom/start/end). "
                        "If set, restrict both null and hit-scoring passes to "
                        "this subset (from 07_0d).")
    p.add_argument("--usable-peaks-dir", default=None,
                   help="Optional dir produced by 07_0d. Per-sample BED is "
                        "auto-derived as "
                        "<dir>/<cross>/<cell_type>/<coord_system>_coord.bed.")
    return p.parse_args()


def load_job_info(jobs_tsv: str, task_id: int) -> pd.Series:
    jobs = pd.read_csv(jobs_tsv, sep="\t")
    if task_id >= len(jobs):
        print(f"[ERROR] task_id {task_id} out of range (max {len(jobs)-1})",
              file=sys.stderr)
        sys.exit(1)
    return jobs.iloc[task_id]


def load_motif_hits(hits_dir: str, cross: str, coord_system: str
                    ) -> dict[str, pd.DataFrame]:
    hits_path = (Path(hits_dir) / cross / f"{coord_system}_coords"
                 / "motif_hits.tsv.gz")
    if not hits_path.exists():
        print(f"[ERROR] motif hits not found: {hits_path}", file=sys.stderr)
        sys.exit(1)
    print(f"[INFO] Loading motif hits from {hits_path}", flush=True)
    t0 = time.time()
    df = pd.read_csv(hits_path, sep="\t", compression="gzip",
                     usecols=["hit_id", "region_str", "hit_center",
                              "motif_id", "motif_name"],
                     dtype={"hit_id": np.int64, "hit_center": np.int64})
    print(f"[INFO]   loaded {len(df):,} hits in {time.time()-t0:.1f}s",
          flush=True)
    grouped = {rstr: grp for rstr, grp in df.groupby("region_str")}
    print(f"[INFO]   {len(grouped):,} regions with hits", flush=True)
    return grouped


def load_acr_coverage(coverage_path: str, sample_id: str) -> dict[str, int]:
    df = pd.read_csv(coverage_path, sep="\t", compression="gzip")
    df = df[df["sample_id"] == sample_id].copy()
    if df.empty:
        print(f"[WARN] No coverage data for {sample_id}", file=sys.stderr)
        return {}
    region_keys = (df["chrom"].astype(str) + ":"
                   + df["start"].astype(str) + "-"
                   + df["end"].astype(str))
    return dict(zip(region_keys, df["n_fragments"].astype(int)))


def extract_fp_at_window(fp_slice: np.ndarray, center_idx: int,
                         window: int, n_pos: int) -> float:
    lo = max(0, center_idx - window)
    hi = min(n_pos, center_idx + window + 1)
    if lo >= hi:
        return np.nan
    chunk = fp_slice[:, lo:hi].copy()
    chunk[~np.isfinite(chunk)] = np.nan
    return float(np.nanmean(chunk))


# ------------------------------------------------------------------
# Pass 1 -- build coverage-stratified null
# ------------------------------------------------------------------

def collect_null_tiles(
    fp_h5: h5py.File,
    tfbs_h5: h5py.File,
    coverage: dict[str, int],
    fp_region_keys: list[str],
    scale_mask: np.ndarray,
    args,
) -> pd.DataFrame:
    """Iterate over all regions with n >= 1; collect (n, fp_raw) pairs at
    TFBS-bottom-bg_percentile% tile positions."""
    n_bg = max(1, int(np.floor(N_TILES * args.bg_percentile / 100.0)))
    print(f"[INFO] Null pass: {n_bg} bottom-TFBS tiles per region",
          flush=True)

    rows = []
    t0 = time.time()
    n_skip_zero = 0
    n_skip_no_tfbs = 0
    n_used = 0
    n_total = len(fp_region_keys)

    for ri, region_str in enumerate(fp_region_keys):
        if (ri + 1) % 5000 == 0 or ri == n_total - 1:
            elapsed = time.time() - t0
            rate = (ri + 1) / elapsed if elapsed > 0 else 0
            print(f"[INFO] Null pass: {ri+1}/{n_total} ({rate:.0f} reg/s, "
                  f"tiles_so_far={len(rows):,})", flush=True)

        n_frags = coverage.get(region_str, 0)
        if n_frags == 0:
            n_skip_zero += 1
            continue
        if region_str not in tfbs_h5["obsm"]:
            n_skip_no_tfbs += 1
            continue

        tfbs_arr = tfbs_h5["obsm"][region_str][:].squeeze()
        bg_tile_indices = np.argsort(tfbs_arr)[:n_bg]
        bg_positions = TILE_BP[bg_tile_indices]

        fp_tensor = fp_h5["obsm"][region_str][:].squeeze()
        fp_slice = fp_tensor[scale_mask, :]
        n_pos = fp_slice.shape[1]

        for bp in bg_positions:
            val = extract_fp_at_window(fp_slice, int(bp), args.hit_window,
                                       n_pos)
            if np.isfinite(val):
                rows.append((n_frags, val))
        n_used += 1

    print(f"[INFO] Null pass: {n_used:,} regions used, "
          f"{n_skip_zero:,} zero-coverage skipped, "
          f"{n_skip_no_tfbs:,} no-TFBS skipped, "
          f"time={time.time()-t0:.1f}s",
          flush=True)

    if len(rows) < 100:
        print(f"[ERROR] Only {len(rows)} null tiles -- cannot proceed",
              file=sys.stderr)
        sys.exit(1)

    return pd.DataFrame(rows, columns=["n", "fp_raw"])


def build_null_model(null_tiles: pd.DataFrame, args) -> pd.DataFrame:
    """Bin null_tiles by n into n_bins quantile bins.  Return per-bin
    table with median, mad, n_med, n_min, n_max, n_tiles."""
    df = null_tiles.copy()
    # Quantile bins.  If many ties (low-n regions), bins may collapse.
    try:
        df["bin"] = pd.qcut(df["n"], q=args.n_bins, duplicates="drop",
                            labels=False)
    except ValueError as e:
        print(f"[ERROR] qcut failed: {e}", file=sys.stderr)
        sys.exit(1)
    df = df[df["bin"].notna()].copy()
    df["bin"] = df["bin"].astype(int)

    # Resolve --use-mean (deprecated) into --spread sd
    spread = args.spread
    if args.use_mean and spread == "iqr":
        spread = "sd"
    print(f"[INFO] Null model: centre={'mean' if spread == 'sd' else 'median'}, "
          f"spread={spread}", flush=True)

    if spread == "sd":
        agg = df.groupby("bin").agg(
            n_med=("n", "median"),
            n_min=("n", "min"),
            n_max=("n", "max"),
            n_tiles=("fp_raw", "count"),
            null_med=("fp_raw", "mean"),
            null_mad=("fp_raw", "std"),
        ).reset_index()
    elif spread == "mad":
        # Median + 1.4826*MAD (legacy robust; degenerate at >=50% zeros)
        def _mad(x):
            return float(1.4826 * np.median(np.abs(x - np.median(x))))
        agg = df.groupby("bin").agg(
            n_med=("n", "median"),
            n_min=("n", "min"),
            n_max=("n", "max"),
            n_tiles=("fp_raw", "count"),
            null_med=("fp_raw", "median"),
            null_mad=("fp_raw", _mad),
        ).reset_index()
    else:  # iqr
        # Median + IQR/1.349 (non-degenerate up to 75% zeros; default)
        def _iqr(x):
            q75, q25 = np.percentile(x, [75, 25])
            return float((q75 - q25) / 1.349)
        agg = df.groupby("bin").agg(
            n_med=("n", "median"),
            n_min=("n", "min"),
            n_max=("n", "max"),
            n_tiles=("fp_raw", "count"),
            null_med=("fp_raw", "median"),
            null_mad=("fp_raw", _iqr),
        ).reset_index()

    # Stash the resolved spread on args for the in-loop merge logic below
    args._resolved_spread = spread

    # Sanity: if any bin has null_mad <= 0 (constant fp_raw), set it to a
    # small floor based on neighbouring bins.
    bad = agg["null_mad"] <= 0
    if bad.any():
        floor = agg.loc[~bad, "null_mad"].median()
        print(f"[WARN] {bad.sum()} bins had null_mad <= 0; "
              f"flooring to {floor:.4f}", flush=True)
        agg.loc[bad, "null_mad"] = floor

    # If small bins, merge upward (rare with 1M+ tiles and 20 bins)
    while (agg["n_tiles"] < args.min_bin_tiles).any() and len(agg) > 1:
        i_small = agg["n_tiles"].idxmin()
        if i_small == len(agg) - 1:
            # last bin -- merge with previous
            i_keep, i_drop = i_small - 1, i_small
        else:
            i_keep, i_drop = i_small + 1, i_small
        # Re-aggregate using the tiles from both
        mask = df["bin"].isin([agg.loc[i_keep, "bin"],
                                agg.loc[i_drop, "bin"]])
        sub = df.loc[mask, "fp_raw"]
        if args._resolved_spread == "sd":
            new_med = float(sub.mean())
            new_mad = float(sub.std())
        elif args._resolved_spread == "mad":
            new_med = float(np.median(sub))
            new_mad = float(1.4826 * np.median(np.abs(sub - new_med)))
        else:  # iqr
            new_med = float(np.median(sub))
            q75, q25 = np.percentile(sub, [75, 25])
            new_mad = float((q75 - q25) / 1.349)
        agg.loc[i_keep, "null_med"] = new_med
        agg.loc[i_keep, "null_mad"] = new_mad if new_mad > 0 else agg["null_mad"].median()
        agg.loc[i_keep, "n_min"] = min(agg.loc[i_keep, "n_min"],
                                        agg.loc[i_drop, "n_min"])
        agg.loc[i_keep, "n_max"] = max(agg.loc[i_keep, "n_max"],
                                        agg.loc[i_drop, "n_max"])
        agg.loc[i_keep, "n_tiles"] = int(agg.loc[i_keep, "n_tiles"]
                                          + agg.loc[i_drop, "n_tiles"])
        # Update the source df bin labels so subsequent merges work
        df.loc[df["bin"] == agg.loc[i_drop, "bin"], "bin"] = (
            agg.loc[i_keep, "bin"]
        )
        agg = agg.drop(index=i_drop).reset_index(drop=True)

    # Sort by n_med for clean lookups
    agg = agg.sort_values("n_med").reset_index(drop=True)
    agg["bin_idx"] = np.arange(len(agg), dtype=int)
    return agg


# ------------------------------------------------------------------
# Pass 2 -- score hits
# ------------------------------------------------------------------

def lookup_bin(n: int, null_model: pd.DataFrame) -> int:
    """Return the bin index whose [n_min, n_max] range contains n.
    For n outside the model range, use first/last bin."""
    n_meds = null_model["n_med"].values
    if n <= n_meds[0]:
        return 0
    if n >= n_meds[-1]:
        return len(n_meds) - 1
    # Find the bin whose [n_min, n_max] contains n; if none, use nearest n_med
    contains = (null_model["n_min"] <= n) & (null_model["n_max"] >= n)
    if contains.any():
        return int(null_model.loc[contains, "bin_idx"].iloc[0])
    # Fall back to nearest n_med (handles gaps from merged bins)
    return int(np.argmin(np.abs(n_meds - n)))


def score_hits_in_region(
    region_str: str,
    hits_df: pd.DataFrame,
    fp_h5: h5py.File,
    scale_mask: np.ndarray,
    n_frags: int,
    null_model: pd.DataFrame,
    args,
) -> list[dict]:
    chrom_rest = region_str.split(":")
    start_end = chrom_rest[1].split("-")
    region_start = int(start_end[0])

    fp_tensor = fp_h5["obsm"][region_str][:].squeeze()
    fp_slice = fp_tensor[scale_mask, :]
    n_pos = fp_slice.shape[1]

    # Coverage -> bin -> null stats (constant across all hits in this region)
    if n_frags == 0:
        bin_idx = -1
        null_med = np.nan
        null_mad = np.nan
        bin_n_med = np.nan
    else:
        bin_idx = lookup_bin(n_frags, null_model)
        null_med = float(null_model.loc[bin_idx, "null_med"])
        null_mad = float(null_model.loc[bin_idx, "null_mad"])
        bin_n_med = float(null_model.loc[bin_idx, "n_med"])

    rows = []
    for row in hits_df.itertuples(index=False):
        center_idx = int(row.hit_center) - region_start
        if center_idx < 0 or center_idx >= n_pos:
            rows.append({
                "hit_id": int(row.hit_id),
                "region_str": region_str,
                "motif_id": row.motif_id,
                "motif_name": row.motif_name,
                "hit_center": int(row.hit_center),
                "fp_raw": np.nan,
                "n_fragments": n_frags,
                "bin_idx": bin_idx,
                "bin_n_med": bin_n_med,
                "bin_null_med": null_med,
                "bin_null_mad": null_mad,
                "z_global": np.nan,
                "fp_present": False,
            })
            continue

        fp_raw = extract_fp_at_window(fp_slice, center_idx,
                                       args.hit_window, n_pos)
        if np.isfinite(fp_raw) and np.isfinite(null_mad) and null_mad > 0:
            z_global = (fp_raw - null_med) / null_mad
        else:
            z_global = np.nan
        fp_present = bool(np.isfinite(z_global)
                          and z_global > args.presence_z)

        rows.append({
            "hit_id": int(row.hit_id),
            "region_str": region_str,
            "motif_id": row.motif_id,
            "motif_name": row.motif_name,
            "hit_center": int(row.hit_center),
            "fp_raw": fp_raw,
            "n_fragments": n_frags,
            "bin_idx": bin_idx,
            "bin_n_med": bin_n_med,
            "bin_null_med": null_med,
            "bin_null_mad": null_mad,
            "z_global": z_global,
            "fp_present": fp_present,
        })
    return rows


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

def main():
    args = parse_args()

    job = load_job_info(args.jobs_tsv, args.task_id)
    sample_id = job["sample_id"]
    cross = job["cross"]
    coord_system = job["coord_system"]
    fp_outdir = Path(job["fp_outdir"])

    print(f"[INFO] Task {args.task_id}: {sample_id} "
          f"({cross}/{coord_system})", flush=True)
    print(f"[INFO] Output root: {args.outdir}")
    print(f"[INFO] Null model: {args.n_bins} bins, "
          f"{'mean+sd' if args.use_mean else 'median+MAD'}")
    print(f"[INFO] Presence threshold: z_global > {args.presence_z}")

    hits_by_region = load_motif_hits(args.hits_dir, cross, coord_system)
    coverage = load_acr_coverage(args.coverage, sample_id)
    print(f"[INFO] Coverage: {len(coverage):,} regions with counts",
          flush=True)

    fp_outdir_rel = str(fp_outdir).split("15_Heterosis/")[-1]
    fp_outdir_local = PROJECT_ROOT / fp_outdir_rel
    # Prefer the new −log10(p) output (return_pval=true); fall back to legacy z-score.
    fp_path = fp_outdir_local / "printer_supp" / "FP_neglog10p__ALL.h5ad"
    if not fp_path.exists():
        fp_path = fp_outdir_local / "printer_supp" / "FP_zscore__ALL.h5ad"
    tfbs_path = fp_outdir_local / "printer_supp" / "TFBS__ALL.h5ad"
    for p in (fp_path, tfbs_path):
        if not p.exists():
            print(f"[ERROR] file not found: {p}", file=sys.stderr)
            sys.exit(1)
    print(f"[INFO] FP input: {fp_path.name}", flush=True)

    fp_h5 = h5py.File(str(fp_path), "r")
    tfbs_h5 = h5py.File(str(tfbs_path), "r")

    scales = fp_h5["uns"]["scales"][:]
    scale_mask = (scales >= args.scales_min) & (scales <= args.scales_max)
    sel_scales = scales[scale_mask]
    print(f"[INFO] Scales: {len(sel_scales)} selected "
          f"({sel_scales[0]}-{sel_scales[-1]}bp)", flush=True)

    fp_region_keys = list(fp_h5["obsm"].keys())

    usable_bed_path: Path | None = None
    if args.usable_peaks and args.usable_peaks_dir:
        print("[ERROR] pass only one of --usable-peaks / --usable-peaks-dir",
              file=sys.stderr)
        sys.exit(1)
    if args.usable_peaks:
        usable_bed_path = Path(args.usable_peaks)
    elif args.usable_peaks_dir:
        import re as _re
        ct_match = _re.search(r"-C(\d+)", sample_id)
        if not ct_match:
            print(f"[ERROR] could not extract cell type from {sample_id}",
                  file=sys.stderr)
            sys.exit(1)
        cell_type = f"C{ct_match.group(1)}"
        usable_bed_path = (Path(args.usable_peaks_dir) / cross / cell_type
                            / f"{coord_system}_coord.bed")
        print(f"[INFO] auto-derived usable-peaks: {usable_bed_path}", flush=True)

    if usable_bed_path is not None:
        if not usable_bed_path.exists():
            print(f"[ERROR] usable-peaks file not found: {usable_bed_path}",
                  file=sys.stderr)
            sys.exit(1)
        usable = set()
        with open(usable_bed_path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                usable.add(f"{parts[0]}:{parts[1]}-{parts[2]}")
        n_before = len(fp_region_keys)
        fp_region_keys = [k for k in fp_region_keys if k in usable]
        print(f"[INFO] usable-peaks filter: {n_before} -> "
              f"{len(fp_region_keys)} regions "
              f"(BED has {len(usable):,} entries)", flush=True)
        if not fp_region_keys:
            print("[ERROR] usable-peaks filter left 0 regions", file=sys.stderr)
            sys.exit(1)

    # ── PASS 1 ──
    print("\n[INFO] === PASS 1: Collecting null tiles ===", flush=True)
    null_tiles = collect_null_tiles(
        fp_h5=fp_h5, tfbs_h5=tfbs_h5, coverage=coverage,
        fp_region_keys=fp_region_keys, scale_mask=scale_mask, args=args,
    )
    tfbs_h5.close()

    print("\n[INFO] === Building null model ===", flush=True)
    null_model = build_null_model(null_tiles, args)
    print(f"[INFO] Null model ({len(null_model)} bins):")
    print(null_model.to_string(index=False))

    # ── PASS 2 ──
    print("\n[INFO] === PASS 2: Scoring motif hits ===", flush=True)
    all_rows = []
    n_total = len(fp_region_keys)
    n_no_cov = 0
    t0 = time.time()
    for ri, region_str in enumerate(fp_region_keys):
        if (ri + 1) % 5000 == 0 or ri == n_total - 1:
            elapsed = time.time() - t0
            rate = (ri + 1) / elapsed if elapsed > 0 else 0
            print(f"[INFO] region {ri+1}/{n_total} ({rate:.0f} reg/s, "
                  f"hits_so_far={len(all_rows):,})", flush=True)
        if region_str not in hits_by_region:
            continue
        n_frags = coverage.get(region_str, 0)
        if n_frags == 0:
            n_no_cov += 1
        rows = score_hits_in_region(
            region_str=region_str, hits_df=hits_by_region[region_str],
            fp_h5=fp_h5, scale_mask=scale_mask, n_frags=n_frags,
            null_model=null_model, args=args,
        )
        all_rows.extend(rows)
    fp_h5.close()

    # ── Write outputs ──
    out_dir = (Path(args.outdir) / cross / f"{coord_system}_coords"
               / sample_id)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Null model dump
    null_path = out_dir / "null_model.tsv"
    null_model.to_csv(null_path, sep="\t", index=False, float_format="%.6g")
    print(f"\n[DONE] null model -> {null_path}")

    # Hit table
    out_path = out_dir / "hit_fp_scores.tsv.gz"
    df_out = pd.DataFrame(all_rows)
    with gzip.open(str(out_path), "wt") as f:
        df_out.to_csv(f, sep="\t", index=False, float_format="%.6g")
    print(f"[DONE] hit_fp_scores -> {out_path}")

    # ── Summary ──
    n_present = int(df_out["fp_present"].sum())
    pct = 100.0 * n_present / max(1, len(df_out))
    print(f"\n  Hits:           {len(df_out):,}")
    print(f"  Regions no cov: {n_no_cov:,}")
    print(f"  Presence:       {n_present:,} / {len(df_out):,} ({pct:.2f}%)")

    if len(df_out) > 0:
        z_vals = df_out["z_global"].dropna()
        if len(z_vals) > 0:
            pcts = np.percentile(z_vals, [0, 25, 50, 75, 100])
            print(f"  z_global: min={pcts[0]:.3f}  Q1={pcts[1]:.3f}  "
                  f"med={pcts[2]:.3f}  Q3={pcts[3]:.3f}  max={pcts[4]:.3f}")

        try:
            df_out["depth_q"] = pd.qcut(
                df_out["n_fragments"].clip(lower=0), q=4,
                labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
            for q, grp in df_out.groupby("depth_q", observed=True):
                qn = int(grp["fp_present"].sum())
                qpct = 100.0 * qn / max(1, len(grp))
                flo = int(grp["n_fragments"].min())
                fhi = int(grp["n_fragments"].max())
                print(f"    {q} (frags {flo}-{fhi}): "
                      f"{qn:,}/{len(grp):,} ({qpct:.2f}%)")
        except ValueError:
            print("    (depth quartile split failed -- too many ties)")

    print(f"  Wall time (pass 2): {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
