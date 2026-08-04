"""
04a_make_regions.py — Resize peaks to 2000 bp and write per-coordinate BEDs.

For each cross, reads the paired-coordinate peaks.tsv, extracts per-parent
coordinates, resizes to 2000 bp centered on the peak midpoint, clips to
chromosome boundaries, and writes:

    regions/{cross}/{parent}_coords/regions_2000bp.bed   (0-based BED)
    regions/{cross}_peak_pairing.tsv                      (peak_id ↔ B73 ↔ P2 regions)

The pairing map is used in Phase 4 to pair FP scores across coordinate systems
for the Wittkopp cis/trans decomposition.

Usage (from project root):
    python 04a_make_regions.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"


def load_chrom_sizes(fai_path: Path) -> dict[str, int]:
    """Read chromosome lengths from a .fai index."""
    sizes = {}
    with open(fai_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            sizes[parts[0]] = int(parts[1])
    return sizes


def resize_regions(
    df: pd.DataFrame, width: int, chrom_sizes: dict[str, int]
) -> pd.DataFrame:
    """Resize regions to `width` bp centered on midpoint, clipped to chrom bounds.

    After clipping to chromosome boundaries, the window is shifted to maintain
    the target width (scPrinter requires uniform-width regions).
    """
    half = width // 2
    mid = (df["start"] + df["end"]) // 2
    new_start = (mid - half).clip(lower=0)
    new_end = mid + half

    # Clip to chromosome length
    max_len = df["chrom"].map(chrom_sizes)
    new_end = new_end.clip(upper=max_len)

    # Enforce uniform width: if clipped at start, extend end; if clipped at end,
    # pull start back.  Chromosomes are >>2000 bp so both adjustments always fit.
    new_end = new_end.clip(lower=new_start + width)
    new_end = new_end.clip(upper=max_len)
    new_start = (new_end - width).clip(lower=0)

    out = df.copy()
    out["start"] = new_start.astype(int)
    out["end"] = new_end.astype(int)
    out["width"] = out["end"] - out["start"]
    return out


def main() -> int:
    with open(CONFIG_PATH) as f:
        config = yaml.safe_load(f)

    region_width = int(config["footprinting"]["region_width"])
    genomes_dir = PROJECT_ROOT / config["paths"]["split_genomes"]
    regions_dir = PROJECT_ROOT / config["paths"]["regions_dir"]
    drop_chroms = set(config["split"]["drop_chromosomes"])

    for cross_name, cross_cfg in config["crosses"].items():
        p1 = cross_cfg["parents"]["P1"]  # always B73
        p2 = cross_cfg["parents"]["P2"]  # Ki3 or Oh43
        peaks_path = PROJECT_ROOT / cross_cfg["peaks_tsv"]

        print(f"\n[INFO] === {cross_name} ===")
        print(f"[INFO] peaks: {peaks_path}")

        # Load chromosome sizes for both parents
        p1_sizes = load_chrom_sizes(genomes_dir / f"{p1}.fa.fai")
        p2_sizes = load_chrom_sizes(genomes_dir / f"{p2}.fa.fai")

        # Read peaks — only need coordinate + ID columns
        peaks = pd.read_csv(peaks_path, sep="\t")
        n_total = len(peaks)

        # Identify coordinate columns
        b73_cols = {"chrom": f"{p1}_chr", "start": f"{p1}_start", "end": f"{p1}_end"}
        p2_cols = {"chrom": f"{p2}_chr", "start": f"{p2}_start", "end": f"{p2}_end"}

        # Drop Mt/Pt peaks (check both coordinate systems)
        mask_drop = (
            peaks[b73_cols["chrom"]].isin(drop_chroms)
            | peaks[p2_cols["chrom"]].isin(drop_chroms)
        )
        n_dropped = mask_drop.sum()
        peaks = peaks[~mask_drop].copy()
        print(f"[INFO] {n_total} total peaks, dropped {n_dropped} on Mt/Pt, kept {len(peaks)}")

        # --- B73 (P1) regions ---
        b73_df = pd.DataFrame({
            "chrom": peaks[b73_cols["chrom"]],
            "start": peaks[b73_cols["start"]].astype(int),
            "end": peaks[b73_cols["end"]].astype(int),
            "peak_id": peaks["peak_id"],
        })
        b73_resized = resize_regions(b73_df, region_width, p1_sizes)

        out_dir = regions_dir / cross_name / f"{p1}_coords"
        out_dir.mkdir(parents=True, exist_ok=True)
        bed_path = out_dir / "regions_2000bp.bed"
        b73_resized[["chrom", "start", "end"]].to_csv(
            bed_path, sep="\t", header=False, index=False
        )
        n_clipped = (b73_resized["width"] < region_width).sum()
        print(f"[OK] {p1} coords: {bed_path} ({len(b73_resized)} regions, {n_clipped} edge-clipped)")

        # --- P2 regions ---
        p2_df = pd.DataFrame({
            "chrom": peaks[p2_cols["chrom"]],
            "start": peaks[p2_cols["start"]].astype(int),
            "end": peaks[p2_cols["end"]].astype(int),
            "peak_id": peaks["peak_id"],
        })
        p2_resized = resize_regions(p2_df, region_width, p2_sizes)

        out_dir = regions_dir / cross_name / f"{p2}_coords"
        out_dir.mkdir(parents=True, exist_ok=True)
        bed_path = out_dir / "regions_2000bp.bed"
        p2_resized[["chrom", "start", "end"]].to_csv(
            bed_path, sep="\t", header=False, index=False
        )
        n_clipped = (p2_resized["width"] < region_width).sum()
        print(f"[OK] {p2} coords: {bed_path} ({len(p2_resized)} regions, {n_clipped} edge-clipped)")

        # --- Peak pairing map ---
        pairing = pd.DataFrame({
            "peak_id": peaks["peak_id"].values,
            f"{p1}_region": (
                b73_resized["chrom"] + ":" +
                b73_resized["start"].astype(str) + "-" +
                b73_resized["end"].astype(str)
            ).values,
            f"{p2}_region": (
                p2_resized["chrom"] + ":" +
                p2_resized["start"].astype(str) + "-" +
                p2_resized["end"].astype(str)
            ).values,
        })
        pairing_path = regions_dir / f"{cross_name}_peak_pairing.tsv"
        pairing.to_csv(pairing_path, sep="\t", index=False)
        print(f"[OK] pairing: {pairing_path} ({len(pairing)} rows)")

    print("\n[DONE] all regions written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
