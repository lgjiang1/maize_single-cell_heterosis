"""
01b_make_jobs.py — emit the SLURM array job table for 01b_split_and_rename.

Reads `config/samples.tsv` and `config/config.yaml` to derive the 168 (input,
target_parent, output) tuples that 01b_split_and_rename.py will be invoked on.

Each parental BAM produces 1 invocation. Each F1 BAM produces 2 invocations
(one per parental allele). Total: 56 + 112 = 168 rows.

Output (default): `_logs/01b_jobs.tsv`. Columns: input_bam, target_parent, output_bam.

Usage (run from project root):
    python 00_scripts/01_split/01b_make_jobs.py
    python 00_scripts/01_split/01b_make_jobs.py --output _logs/01b_jobs.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"


def load_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def load_samples(samples_tsv: Path) -> list[dict]:
    rows: list[dict] = []
    with open(samples_tsv) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            values = line.rstrip("\n").split("\t")
            rows.append(dict(zip(header, values)))
    return rows


def derive_jobs(config: dict, samples: list[dict]) -> list[tuple[str, str, str]]:
    """Build (input_bam, target_parent, output_bam) tuples for all 168 jobs."""
    jobs: list[tuple[str, str, str]] = []

    for row in samples:
        cross = row["cross"]
        role = row["role"]
        name = row["name"]
        cell_type = row["cell_type"]
        input_bam = row["input_bam"]

        cross_cfg = config["crosses"][cross]
        p1 = cross_cfg["parents"]["P1"]
        p2 = cross_cfg["parents"]["P2"]

        if role in ("P1", "P2"):
            # Parental BAM: 1 invocation, target = the parent itself
            target_parent = name
            output_bam = f"1_split/{cross}/{target_parent}_coords/{target_parent}-{cell_type}.bam"
            jobs.append((input_bam, target_parent, output_bam))
        elif role in ("F1_fwd", "F1_rev"):
            # F1 BAM: 2 invocations, one per parental allele
            for target_parent in (p1, p2):
                output_bam = (
                    f"1_split/{cross}/{target_parent}_coords/"
                    f"{name}-{cell_type}.{target_parent}allele.bam"
                )
                jobs.append((input_bam, target_parent, output_bam))
        else:
            print(f"[WARN] unknown role: {role} (row skipped)", file=sys.stderr)

    return jobs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="_logs/01b_jobs.tsv",
        help="Output TSV path (default _logs/01b_jobs.tsv)",
    )
    args = parser.parse_args()

    config = load_config()
    samples = load_samples(PROJECT_ROOT / config["paths"]["sample_sheet"])
    jobs = derive_jobs(config, samples)

    out_path = PROJECT_ROOT / args.output
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, "w") as f:
        for input_bam, target_parent, output_bam in jobs:
            f.write(f"{input_bam}\t{target_parent}\t{output_bam}\n")

    print(f"[INFO] derived {len(jobs)} jobs from {len(samples)} samples")
    print(f"[OK]   wrote {out_path.relative_to(PROJECT_ROOT)}")
    print(f"[INFO] sbatch with --array=0-{len(jobs) - 1}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
