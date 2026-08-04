"""
04b_make_fp_jobs.py — Build the fragment extraction + scPrinter job table.

Walks 2_equalized/ to discover all 168 BAMs, pairs each with the correct
genome object and region BED, and writes a single job table used by both
04b_bam_to_frags.sh and 04c_run_print.sh.

Writes:
    _logs/04b_jobs.tsv — 168 rows, columns:
        cross  coord_system  parent  sample_id  eq_bam  frag_out  genome_obj  region_bed  fp_outdir

Usage (from project root):
    python 04b_make_fp_jobs.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"


def main() -> int:
    with open(CONFIG_PATH) as f:
        config = yaml.safe_load(f)

    eq_dir = PROJECT_ROOT / config["paths"]["split_dir"]  # read directly from 1_split/ (no equalization)
    frag_dir = PROJECT_ROOT / config["paths"]["fragments_dir"]
    fp_dir = PROJECT_ROOT / config["paths"]["footprinting_dir"]
    bias_dir = PROJECT_ROOT / config["paths"]["bias_model_dir"]
    regions_dir = PROJECT_ROOT / config["paths"]["regions_dir"]
    logs_dir = PROJECT_ROOT / config["paths"]["logs_dir"]

    rows = []
    for cross_name, cross_cfg in config["crosses"].items():
        parents = [cross_cfg["parents"]["P1"], cross_cfg["parents"]["P2"]]

        for parent in parents:
            coord_dir = eq_dir / cross_name / f"{parent}_coords"
            if not coord_dir.exists():
                print(f"[WARN] missing: {coord_dir}", file=sys.stderr)
                continue

            genome_obj = str(bias_dir / parent / f"{parent}_genome_OBJ")
            region_bed = str(regions_dir / cross_name / f"{parent}_coords" / "regions_2000bp.bed")

            bams = sorted(coord_dir.glob("*.bam"))
            # Filter out .bam.bai files and resolve symlinks
            bams = [b for b in bams if b.suffix == ".bam"]

            for bam_path in bams:
                stem = bam_path.stem  # e.g., B73-C1, BxK-C1.B73allele
                sample_id = f"{cross_name}_{parent}_{stem}"

                frag_out = str(
                    frag_dir / cross_name / f"{parent}_coords" / f"{stem}.frags.tsv.gz"
                )
                fp_outdir = str(
                    fp_dir / cross_name / f"{parent}_coords" / stem
                )

                rows.append({
                    "cross": cross_name,
                    "coord_system": parent,
                    "parent": parent,
                    "sample_id": sample_id,
                    "eq_bam": str(bam_path),
                    "frag_out": frag_out,
                    "genome_obj": genome_obj,
                    "region_bed": region_bed,
                    "fp_outdir": fp_outdir,
                })

    if not rows:
        print("[ERROR] no BAMs found in 2_equalized/", file=sys.stderr)
        return 1

    logs_dir.mkdir(parents=True, exist_ok=True)
    out_path = logs_dir / "04b_jobs.tsv"
    with open(out_path, "w") as f:
        header = "\t".join(rows[0].keys())
        f.write(header + "\n")
        for row in rows:
            f.write("\t".join(row.values()) + "\n")

    print(f"[OK] wrote {out_path} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
