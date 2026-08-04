"""
03a_make_eq_jobs.py — Build the equalization job table for Phase 2.

For each of 168 split BAMs, decide whether to symlink (already at or below
parental floor) or subsample (F1 allele exceeds parental floor).

Reads:
    qc/1_split/depth_table.tsv      — per-BAM read counts
    qc/1_split/parental_floor.tsv   — per-(cross, cell_type, coord_system) floor

Writes:
    _logs/03a_jobs.tsv  — 168 rows, one per SLURM array task
        columns: input_bam  output_bam  action  seed_fraction  input_reads  floor_reads

Usage (from project root):
    python 03a_make_eq_jobs.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"


def main() -> int:
    with open(CONFIG_PATH) as f:
        config = yaml.safe_load(f)

    seed = int(config["equalize"]["seed"])
    split_dir = PROJECT_ROOT / config["paths"]["split_dir"]
    eq_dir = PROJECT_ROOT / config["paths"]["equalized_dir"]
    logs_dir = PROJECT_ROOT / config["paths"]["logs_dir"]

    dt = pd.read_csv(PROJECT_ROOT / "qc/1_split/depth_table.tsv", sep="\t")
    pf = pd.read_csv(PROJECT_ROOT / "qc/1_split/parental_floor.tsv", sep="\t")

    # Build lookup: (cross, cell_type, coord_system) -> parental floor reads
    floor_lookup = {}
    for _, row in pf.iterrows():
        key = (row["cross"], row["cell_type"], row["coord_system"])
        floor_lookup[key] = int(row["parent_reads"])

    jobs = []
    for _, row in dt.iterrows():
        cross = row["cross"]
        coord = row["coord_system"]
        cell_type = row["cell_type"]
        role = row["role"]
        output_reads = int(row["output_reads"])
        input_bam = str(row["output_bam"])  # the split BAM is our input

        # Derive output path under 2_equalized/
        # Input: 1_split/B-K/B73_coords/B73-C1.bam
        # Output: 2_equalized/B-K/B73_coords/B73-C1.bam
        rel = Path(input_bam).relative_to(
            Path(input_bam).parents[2]  # strip up to 1_split/
        )
        # Handle both absolute and relative paths
        parts = Path(input_bam).parts
        try:
            idx = parts.index("1_split")
            rel = Path(*parts[idx + 1:])
        except ValueError:
            print(f"[ERROR] cannot parse path: {input_bam}", file=sys.stderr)
            return 1

        output_bam = str(eq_dir / rel)

        floor_key = (cross, cell_type, coord)
        floor_reads = floor_lookup.get(floor_key)
        if floor_reads is None:
            print(f"[ERROR] no floor for {floor_key}", file=sys.stderr)
            return 1

        if role in ("P1", "P2"):
            action = "symlink"
            seed_fraction = ""
        elif output_reads <= floor_reads:
            action = "symlink"
            seed_fraction = ""
        else:
            action = "subsample"
            frac = floor_reads / output_reads
            # samtools format: SEED.DECIMAL (e.g., 42.735614)
            frac_str = f"{frac:.6f}"
            seed_fraction = f"{seed}.{frac_str[2:]}"  # strip "0."

        jobs.append({
            "input_bam": input_bam,
            "output_bam": output_bam,
            "action": action,
            "seed_fraction": seed_fraction,
            "input_reads": output_reads,
            "floor_reads": floor_reads,
        })

    df = pd.DataFrame(jobs)
    logs_dir.mkdir(parents=True, exist_ok=True)
    out_path = logs_dir / "03a_jobs.tsv"
    df.to_csv(out_path, sep="\t", index=False)

    n_sym = (df["action"] == "symlink").sum()
    n_sub = (df["action"] == "subsample").sum()
    print(f"[OK] wrote {out_path} ({len(df)} rows: {n_sym} symlinks, {n_sub} subsamples)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
