"""
00a_build_sample_sheet.py — build the project sample sheet for `15_Heterosis`.

Walks the read-only input directories under `0_bams/`, parses BAM filenames
against the known patterns, and emits a 112-row TSV at `config/samples.tsv`.

Usage (run from project root):
    python 00a_build_sample_sheet.py

Output schema (`config/samples.tsv`):
    cross       B-K | B-O
    role        P1 | P2 | F1_fwd | F1_rev
    name        B73 | Ki3 | Oh43 | BxK | KxB | BxO | OxB
    cell_type   C1 .. C14
    input_bam   path relative to project root

Filename conventions handled:
    {P1|P2}-C{n}.renamed.bam   → parental BAMs (header has both parents' SQs,
                                  reads only on this parent's chrs)
    {F1}-C{n}.bam              → F1 BAMs (mapped to merged genome, both alleles)
"""
from __future__ import annotations

import sys
from pathlib import Path

import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"
OUTPUT_TSV = PROJECT_ROOT / "config" / "samples.tsv"


def load_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def parental_bam_path(cross: str, parent_name: str, cell_type: str) -> Path:
    """Path to a parental (.renamed.bam) BAM under 0_bams/."""
    return PROJECT_ROOT / "0_bams" / f"0_{cross}_input" / f"{parent_name}-{cell_type}.renamed.bam"


def f1_bam_path(cross: str, f1_name: str, cell_type: str) -> Path:
    """Path to an F1 BAM under 0_bams/."""
    return PROJECT_ROOT / "0_bams" / f"0_{cross}_input" / f"{f1_name}-{cell_type}.bam"


def build_rows(config: dict) -> list[dict]:
    rows: list[dict] = []
    cell_types: list[str] = config["cell_types"]

    for cross, cross_cfg in config["crosses"].items():
        p1_name = cross_cfg["parents"]["P1"]
        p2_name = cross_cfg["parents"]["P2"]
        f1_fwd_name = cross_cfg["f1_directions"]["forward"]
        f1_rev_name = cross_cfg["f1_directions"]["reverse"]

        for cell_type in cell_types:
            rows.append({
                "cross": cross,
                "role": "P1",
                "name": p1_name,
                "cell_type": cell_type,
                "input_bam": str(parental_bam_path(cross, p1_name, cell_type).relative_to(PROJECT_ROOT)),
            })
            rows.append({
                "cross": cross,
                "role": "P2",
                "name": p2_name,
                "cell_type": cell_type,
                "input_bam": str(parental_bam_path(cross, p2_name, cell_type).relative_to(PROJECT_ROOT)),
            })
            rows.append({
                "cross": cross,
                "role": "F1_fwd",
                "name": f1_fwd_name,
                "cell_type": cell_type,
                "input_bam": str(f1_bam_path(cross, f1_fwd_name, cell_type).relative_to(PROJECT_ROOT)),
            })
            rows.append({
                "cross": cross,
                "role": "F1_rev",
                "name": f1_rev_name,
                "cell_type": cell_type,
                "input_bam": str(f1_bam_path(cross, f1_rev_name, cell_type).relative_to(PROJECT_ROOT)),
            })

    return rows


def verify_inputs(rows: list[dict]) -> tuple[int, int]:
    """Check that every input_bam resolves on disk."""
    found = 0
    missing = 0
    for row in rows:
        full = PROJECT_ROOT / row["input_bam"]
        if full.exists():
            found += 1
        else:
            missing += 1
            print(f"  [MISSING] {row['cross']:5s} {row['role']:6s} {row['name']:5s} {row['cell_type']:4s} → {row['input_bam']}", file=sys.stderr)
    return found, missing


def write_tsv(rows: list[dict], path: Path) -> None:
    columns = ["cross", "role", "name", "cell_type", "input_bam"]
    with open(path, "w") as f:
        f.write("\t".join(columns) + "\n")
        for row in rows:
            f.write("\t".join(row[c] for c in columns) + "\n")


def main() -> int:
    config = load_config()
    rows = build_rows(config)

    print(f"[INFO] generated {len(rows)} sample sheet rows")
    print(f"[INFO] verifying input files exist…")

    found, missing = verify_inputs(rows)
    print(f"[INFO]   found:   {found}")
    print(f"[INFO]   missing: {missing}")

    if missing > 0:
        print(f"[ERROR] {missing} input BAMs are missing — aborting before writing sample sheet", file=sys.stderr)
        return 1

    OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(rows, OUTPUT_TSV)
    print(f"[OK]   wrote {OUTPUT_TSV.relative_to(PROJECT_ROOT)} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
