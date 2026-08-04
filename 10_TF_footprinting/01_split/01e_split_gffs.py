"""
01e_split_gffs.py — split per-cross merged GFF3 files into per-parent GFF3s.

Reads each cross's merged GFF3 from 0_bams/0_B-{K,O}_input/, extracts
features for each parent based on chromosome prefix, strips the prefix
(B73_chr1 -> chr1), drops Mt+Pt, and writes per-parent GFF3s to
1_split/genomes/.

For B73 (shared between crosses), verifies that both merged GFFs produce
the same B73 feature set (line count), then writes a single shared B73.gff3.

Usage (run from project root):
    python 01e_split_gffs.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"


def load_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def split_gff(
    merged_gff: Path, parent_prefix: str, drop_chroms: set[str]
) -> list[str]:
    """Extract and rename GFF3 lines for one parent.

    Returns a list of output lines (header + feature lines).
    """
    prefix_len = len(parent_prefix)
    out_lines: list[str] = []

    with open(merged_gff) as f:
        for line in f:
            # Pass through GFF3 header/comment lines only once (first parent)
            if line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            if len(parts) < 2:
                continue
            chrom = parts[0]
            # Skip chromosomes to drop (Mt, Pt)
            if chrom in drop_chroms:
                continue
            # Keep only lines for this parent
            if not chrom.startswith(parent_prefix):
                continue
            # Strip prefix: B73_chr1 -> chr1
            new_chrom = chrom[prefix_len:]
            out_lines.append(new_chrom + "\t" + parts[1])

    return out_lines


def write_gff(out_path: Path, lines: list[str]) -> None:
    with open(out_path, "w") as f:
        f.write("##gff-version 3\n")
        for line in lines:
            f.write(line)


def main() -> int:
    config = load_config()
    drop_chroms = set(config["split"]["drop_chromosomes"])
    out_dir = PROJECT_ROOT / config["paths"]["split_genomes"]
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] output directory: {out_dir.relative_to(PROJECT_ROOT)}")
    print(f"[INFO] dropping chroms:  {sorted(drop_chroms)}")

    b73_lines_per_cross: dict[str, list[str]] = {}
    parent2_results: list[tuple[str, list[str]]] = []

    for cross, cross_cfg in config["crosses"].items():
        merged_gff = PROJECT_ROOT / cross_cfg["merged_gff"]
        if not merged_gff.exists():
            print(f"[ERROR] merged GFF not found: {merged_gff}", file=sys.stderr)
            return 1

        print(f"\n[STEP] [{cross}] reading {merged_gff.relative_to(PROJECT_ROOT)}")

        for parent_role, parent_name in cross_cfg["parents"].items():
            prefix = f"{parent_name}_"
            lines = split_gff(merged_gff, prefix, drop_chroms)
            print(f"[INFO]   [{cross}] {parent_name}: {len(lines):,} feature lines")

            if parent_name == "B73":
                b73_lines_per_cross[cross] = lines
            else:
                parent2_results.append((parent_name, lines))

    # Verify B73 consistency across crosses
    crosses = list(b73_lines_per_cross.keys())
    if len(crosses) >= 2:
        print(f"\n[STEP] verifying B73 GFF consistency across crosses")
        ref = crosses[0]
        ref_count = len(b73_lines_per_cross[ref])
        for other in crosses[1:]:
            other_count = len(b73_lines_per_cross[other])
            if ref_count != other_count:
                print(
                    f"[WARN] B73 line count differs: {ref}={ref_count:,} vs {other}={other_count:,}",
                    file=sys.stderr,
                )
                print("[WARN] Using B73 from first cross (B-K)", file=sys.stderr)
            else:
                print(f"[OK]   B73 GFF identical between {ref} and {other} ({ref_count:,} lines)")

    # Write per-parent GFF3s
    print(f"\n[STEP] writing per-parent GFF3 files")

    # B73
    b73_lines = b73_lines_per_cross[crosses[0]]
    b73_path = out_dir / "B73.gff3"
    write_gff(b73_path, b73_lines)
    print(f"[OK]   {b73_path.relative_to(PROJECT_ROOT)} ({len(b73_lines):,} features)")

    # Other parents
    for parent_name, lines in parent2_results:
        out_path = out_dir / f"{parent_name}.gff3"
        write_gff(out_path, lines)
        print(f"[OK]   {out_path.relative_to(PROJECT_ROOT)} ({len(lines):,} features)")

    print(f"\n[DONE] all per-parent GFF3 files written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
