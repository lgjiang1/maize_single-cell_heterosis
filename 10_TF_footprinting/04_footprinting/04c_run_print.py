"""
04c_run_print.py — scPrinter import + TFBS + NucBS + multi-scale FP scoring.

One invocation per condition (called from the SLURM array wrapper).
Adapted from 13_Arabidopsis.../03_footprinting/03a_run_print.py and
2_PopulationStress.../3_04_score_grid.py.

Usage:
    python 04c_run_print.py \\
        --frag 3_fragments/B-K/B73_coords/B73-C1.frags.tsv.gz \\
        --genome-obj bias_model/B73/B73_genome_OBJ \\
        --regions regions/B-K/B73_coords/regions_2000bp.bed \\
        --outdir 4_footprinting/B-K/B73_coords/B73-C1 \\
        --sample-id B-K_B73_B73-C1
"""
from __future__ import annotations

import argparse
import gc
import os
import pickle
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import scprinter as scp

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"


def safe_key(s: str) -> str:
    return s.replace("-", "_").replace(".", "_").replace("/", "_")


def supp_dir_of(printer_path: str) -> str:
    base_dir = os.path.dirname(printer_path)
    base_stem = os.path.splitext(os.path.basename(printer_path))[0]
    return os.path.join(base_dir, f"{base_stem}_supp")


def wait_for(paths: list[str], tries: int = 60, sleep_s: int = 2):
    """Wait for one of the paths to appear (backed output may take a moment)."""
    for _ in range(tries):
        for p in paths:
            if os.path.exists(p) and os.path.getsize(p) > 0:
                return p
        time.sleep(sleep_s)
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frag", required=True, help="Fragment file (.tsv.gz)")
    ap.add_argument("--genome-obj", required=True, help="Pickled scprinter.genome.Genome")
    ap.add_argument("--regions", required=True, help="Region BED (0-based, 3-col)")
    ap.add_argument("--regions-subset", default=None,
                    help="Optional 3-col BED of regions to score, in the same "
                         "2 kb coordinate system as --regions. When provided, "
                         "the script inner-joins on (chrom, start, end) and "
                         "computes FP only on the intersection. Use to "
                         "pre-filter low-coverage peaks (from 07_0d output) "
                         "and shrink runtime + h5ad size.")
    ap.add_argument("--outdir", required=True, help="Output directory for this condition")
    ap.add_argument("--sample-id", required=True, help="Sample identifier")
    args = ap.parse_args()

    with open(CONFIG_PATH) as f:
        config = yaml.safe_load(f)

    fp_cfg = config["footprinting"]
    region_width = int(fp_cfg["region_width"])
    context_radius = int(fp_cfg["tfbs_context_radius"])
    downsample = int(fp_cfg["tfbs_downsample"])
    fp_modes = np.arange(int(fp_cfg["fp_modes_min"]), int(fp_cfg["fp_modes_max"]) + 1)

    n_jobs = int(os.environ.get("SLURM_CPUS_PER_TASK", "12"))

    # --- Load genome ---
    print(f"[INFO] loading genome: {args.genome_obj}")
    with open(args.genome_obj, "rb") as f:
        genome = pickle.load(f)
    print(f"[INFO] genome chroms: {list(genome.chrom_sizes.keys())[:5]}...")

    # The genome objects were originally built without a GFF. If gff_file is
    # still None, look for the per-parent GFF3 next to the fasta (produced by
    # 01e_split_gffs.py). scPrinter's import_fragments unconditionally calls
    # snap.metrics.tsse() which needs genome.gff_file.
    if genome.gff_file is None:
        fa_path = genome.fa_file if hasattr(genome, "fa_file") else None
        if fa_path:
            candidate = Path(fa_path).with_suffix(".gff3")
            if candidate.exists():
                genome.gff_file = str(candidate)
                print(f"[INFO] genome.gff_file was None — patched with {candidate}")
        if genome.gff_file is None:
            print("[ERROR] genome.gff_file is None and no per-parent GFF3 found", file=sys.stderr)
            print("[ERROR] run 01e_split_gffs.py first", file=sys.stderr)
            return 1

    # --- Load regions ---
    regions = pd.read_csv(
        args.regions, sep="\t", header=None, usecols=[0, 1, 2],
        names=["Chromosome", "Start", "End"],
    )
    regions["Start"] = regions["Start"].astype(np.int64)
    regions["End"] = regions["End"].astype(np.int64)
    regions.drop_duplicates(subset=["Chromosome", "Start", "End"], inplace=True)
    regions.sort_values(["Chromosome", "Start", "End"], inplace=True, ignore_index=True)
    n_full = len(regions)
    print(f"[INFO] regions: {n_full} from {args.regions}")

    # Optional pre-filter to a subset (e.g. coverage-usable peaks from 07_0d).
    if args.regions_subset:
        subset = pd.read_csv(
            args.regions_subset, sep="\t", header=None, usecols=[0, 1, 2],
            names=["Chromosome", "Start", "End"],
        )
        subset["Start"] = subset["Start"].astype(np.int64)
        subset["End"] = subset["End"].astype(np.int64)
        subset.drop_duplicates(subset=["Chromosome", "Start", "End"],
                                inplace=True)
        regions = regions.merge(subset,
                                 on=["Chromosome", "Start", "End"],
                                 how="inner")
        regions.sort_values(["Chromosome", "Start", "End"],
                             inplace=True, ignore_index=True)
        print(f"[INFO] subset filter: {n_full} -> {len(regions)} regions "
              f"({100.0*len(regions)/max(1,n_full):.1f}%) from "
              f"{args.regions_subset}")
        if regions.empty:
            print("[ERROR] empty region set after subset intersection — "
                  "check that --regions and --regions-subset share the same "
                  "coordinate system", file=sys.stderr)
            return 1

    # --- Prepare output ---
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    wl_dir = outdir / "_wl"
    wl_dir.mkdir(exist_ok=True)

    wl_file = str(wl_dir / "bulk.barcodes.txt")
    with open(wl_file, "w") as f:
        f.write("bulk\n")

    out_h5 = str(outdir / "printer.h5ad")
    if os.path.exists(out_h5):
        os.remove(out_h5)

    sample_id = args.sample_id

    # ── Step 1: Import fragments ──────────────────────────────────────────
    print(f"[STEP] import_fragments for {sample_id}")
    t0 = time.time()
    pr = scp.pp.import_fragments(
        path_to_frags=[args.frag],
        barcodes=[wl_file],
        sample_names=[sample_id],
        savename=out_h5,
        genome=genome,
        auto_detect_shift=bool(fp_cfg["auto_detect_shift"]),
        plus_shift=int(fp_cfg["plus_shift"]),
        minus_shift=int(fp_cfg["minus_shift"]),
        min_num_fragments=0,
        min_tsse=0,
        sorted_by_barcode=False,
        n_jobs=1,
    )
    pr.close()
    gc.collect()
    print(f"[OK] import done ({time.time() - t0:.0f}s)")

    # ── Step 2: Load printer + set up grouping ────────────────────────────
    printer = scp.load_printer(out_h5, genome)
    printer.load_disp_model()

    cell_ids = list(printer.obs_names)
    if "bulk" in cell_ids:
        cell_grouping = [["bulk"]]
    else:
        cell_grouping = [[cell_ids[0]]]
    group_names = [sample_id]

    supp_dir = supp_dir_of(out_h5)
    os.makedirs(supp_dir, exist_ok=True)

    # ── Step 3: TFBS ──────────────────────────────────────────────────────
    printer.load_bindingscore_model("TF", scp.datasets.pretrained_TFBS_model)

    tf_key = safe_key(f"TFBS_{sample_id}")
    print(f"[STEP] TFBS key={tf_key} n_jobs={n_jobs}")
    t0 = time.time()
    scp.tl.get_binding_score(
        printer,
        cell_grouping=cell_grouping,
        group_names=group_names,
        regions=regions,
        model_key="TF",
        n_jobs=n_jobs,
        contextRadius=context_radius,
        region_width=region_width,
        downsample=downsample,
        save_key=tf_key,
        backed=True,
        overwrite=True,
    )
    print(f"[OK] TFBS done ({time.time() - t0:.0f}s)")

    # ── Step 4: NucBS ─────────────────────────────────────────────────────
    printer.load_bindingscore_model("Nuc", scp.datasets.pretrained_NucBS_model)

    nuc_key = safe_key(f"NucBS_{sample_id}")
    print(f"[STEP] NucBS key={nuc_key} n_jobs={n_jobs}")
    t0 = time.time()
    scp.tl.get_binding_score(
        printer,
        cell_grouping=cell_grouping,
        group_names=group_names,
        regions=regions,
        model_key="Nuc",
        n_jobs=n_jobs,
        contextRadius=context_radius,
        region_width=region_width,
        downsample=downsample,
        save_key=nuc_key,
        backed=True,
        overwrite=True,
    )
    print(f"[OK] NucBS done ({time.time() - t0:.0f}s)")

    # ── Step 5: Multi-scale FP ───────────────────────────────────────────
    return_pval = bool(fp_cfg.get("fp_return_pval", False))
    fp_label = "FP_neglog10p" if return_pval else "FP_zscore"
    fp_key = safe_key(f"{fp_label}_{sample_id}")
    print(f"[STEP] FP label={fp_label} key={fp_key} modes={fp_modes[0]}-{fp_modes[-1]} return_pval={return_pval} n_jobs={n_jobs}")
    t0 = time.time()
    scp.tl.get_footprint_score(
        printer,
        cell_grouping=cell_grouping,
        group_names=group_names,
        regions=regions,
        modes=fp_modes,
        n_jobs=n_jobs,
        region_width=region_width,
        save_key=fp_key,
        backed=True,
        overwrite=True,
        return_pval=return_pval,
    )
    print(f"[OK] FP {'pval' if return_pval else 'zscore'} done ({time.time() - t0:.0f}s)")

    printer.close()

    # ── Step 6: Copy backed outputs to canonical names ────────────────────
    for label, key in [("TFBS", tf_key), ("NucBS", nuc_key), (fp_label, fp_key)]:
        src = wait_for([
            os.path.join(supp_dir, f"{key}.h5ad"),
            os.path.join(supp_dir, f"{safe_key(key)}.h5ad"),
        ])
        if src is None:
            print(f"[ERROR] {label} backed output not found in {supp_dir}", file=sys.stderr)
            return 1
        dst = os.path.join(supp_dir, f"{label}__ALL.h5ad")
        if src != dst:
            tmp = dst + ".tmp"
            shutil.copy2(src, tmp)
            os.replace(tmp, dst)
        print(f"[OK] {label} -> {dst}")

    print(f"\n[DONE] {sample_id}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
