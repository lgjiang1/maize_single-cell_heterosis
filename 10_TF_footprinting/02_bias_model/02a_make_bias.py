"""
02a_make_bias.py — predict the per-parent Tn5 bias model and build the
scprinter Genome object.

This is a parameterized version of the original `bias_model/make_maize_bias.py`.
It runs once per parent (B73, Ki3, Oh43) and consumes a per-parent fasta from
`1_split/genomes/{parent}.fa` (produced upstream by `01a_split_genomes.py`).

Important: scPrinter's bias model is **predicted, not trained**. A pretrained
Tn5 sequence-bias CNN ships with scPrinter (trained once on naked-DNA Tn5
insertions from human BACs). We just apply it to the maize sequence — no ATAC
data is needed for this step. Mirrors the Arabidopsis approach at:
    13_Arbidopsis_protoplast/5_TF_FP/v1/3_00_Create_bias_and_genome_obj.py

Outputs (under `bias_model/{parent}/`):
    {parent}_bias.h5     : per-position bias prediction (HDF5)
    {parent}_genome_OBJ  : pickled `scprinter.genome.Genome`

Usage (run from project root):
    python 02a_make_bias.py --parent B73
    python 02a_make_bias.py --parent Ki3
    python 02a_make_bias.py --parent Oh43
"""
from __future__ import annotations

import argparse
import os
import pickle
import random
import sys
from pathlib import Path
from typing import Dict, List

import scprinter as scp
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"


def load_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def parse_fai(fai_path: Path) -> Dict[str, int]:
    """Read a samtools .fai file and return {chromosome_name: length}."""
    sizes: Dict[str, int] = {}
    with open(fai_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


def make_chrom_splits(
    chroms: List[str],
    k: int = 5,
    test_size: int = 2,
    valid_size: int = 1,
    seed: int = 2025,
) -> List[Dict[str, List[str]]]:
    """Rotation-based chromosome cross-validation splits.

    Mirrors the helper in `13_Arbidopsis_protoplast/5_TF_FP/v1/3_00_Create_bias_and_genome_obj.py`.
    """
    n = len(chroms)
    if test_size + valid_size >= n:
        raise ValueError("test_size + valid_size must be < total number of chromosomes")
    rng = random.Random(seed)
    order = chroms[:]
    rng.shuffle(order)
    splits = []
    stride = max(1, n // k)
    for fold in range(k):
        rot = order[fold * stride % n:] + order[: fold * stride % n]
        test = rot[:test_size]
        valid = rot[test_size: test_size + valid_size]
        train = rot[test_size + valid_size:]
        splits.append({"test": test, "valid": valid, "train": train})
    return splits


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--parent",
        required=True,
        choices=["B73", "Ki3", "Oh43"],
        help="Which parent to build the bias model for",
    )
    args = parser.parse_args()

    config = load_config()
    parent = args.parent
    bias_params = config["bias"]

    fasta_path = PROJECT_ROOT / config["paths"]["split_genomes"] / f"{parent}.fa"
    fai_path = fasta_path.with_suffix(fasta_path.suffix + ".fai")

    if not fasta_path.exists():
        print(f"[ERROR] per-parent fasta not found: {fasta_path}", file=sys.stderr)
        print(f"[ERROR] run 01a_split_genomes.py first", file=sys.stderr)
        return 1
    if not fai_path.exists():
        print(f"[ERROR] fasta index not found: {fai_path}", file=sys.stderr)
        return 1

    out_dir = PROJECT_ROOT / config["paths"]["bias_model_dir"] / parent
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] parent:           {parent}")
    print(f"[INFO] fasta:            {fasta_path.relative_to(PROJECT_ROOT)}")
    print(f"[INFO] output directory: {out_dir.relative_to(PROJECT_ROOT)}")

    # Read chromosome sizes from .fai (per-parent fastas use chr1..chr10 only,
    # no Mt/Pt — these were dropped at 01a_split_genomes time)
    chrom_sizes = parse_fai(fai_path)
    print(f"[INFO] chromosomes:      {len(chrom_sizes)}")
    for k_, v in chrom_sizes.items():
        print(f"[INFO]   {k_:8s}  {v:>14,} bp")

    # CV splits over the chromosomes (organelles already excluded)
    splits = make_chrom_splits(
        list(chrom_sizes.keys()),
        k=int(bias_params["cv_folds"]),
        test_size=int(bias_params["cv_test_size"]),
        valid_size=int(bias_params["cv_valid_size"]),
        seed=int(bias_params["cv_seed"]),
    )

    # Step 1: predict Tn5 bias from sequence (uses pretrained CNN)
    bias_path = out_dir / f"{parent}_bias.h5"
    print(f"\n[STEP] predicting genome-wide Tn5 bias -> {bias_path.relative_to(PROJECT_ROOT)}")
    scp.genome.predict_genome_tn5_bias(
        fa_file=str(fasta_path),
        save_name=str(bias_path),
        tn5_model=scp.datasets.pretrained_Tn5_bias_model,
        context_radius=int(bias_params["context_radius"]),
        device=str(bias_params["device"]),
        batch_size=int(bias_params["batch_size"]),
    )
    print(f"[OK]   bias HDF5 written: {bias_path.relative_to(PROJECT_ROOT)}")

    # Step 2: wrap into a scprinter.genome.Genome object
    print("\n[STEP] building scprinter.genome.Genome object")

    # GFF is optional — at scaffold time we don't need it for FP scoring, but
    # if a per-parent GFF exists alongside the fasta we'll wire it in.
    gff_candidate = fasta_path.with_suffix(".gff3")
    gff_file = str(gff_candidate) if gff_candidate.exists() else None
    if gff_file:
        print(f"[INFO] using GFF: {gff_candidate.relative_to(PROJECT_ROOT)}")
    else:
        print("[INFO] no per-parent GFF found — Genome.gff_file=None (OK for FP scoring)")

    genome = scp.genome.Genome(
        name=f"Zm_{parent}",
        chrom_sizes=chrom_sizes,
        gff_file=gff_file,
        fa_file=str(fasta_path),
        bias_file=str(bias_path),
        bg=str(bias_path),
        blacklist_file=None,
        splits=splits,
    )

    genome_save_path = out_dir / f"{parent}_genome_OBJ"
    with open(genome_save_path, "wb") as f:
        pickle.dump(genome, f)

    print(f"[OK]   saved genome object: {genome_save_path.relative_to(PROJECT_ROOT)}")
    print(f"[OK]   genome chrom keys:   {list(genome.chrom_sizes.keys())}")
    print(f"\n[DONE] {parent} bias model + genome object ready")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
