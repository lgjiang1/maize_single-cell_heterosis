"""
01a_split_genomes.py — split per-cross merged fastas into per-parent fastas.

Reads each cross's merged pseudo-diploid fasta from `0_bams/0_B-{K,O}_input/`,
extracts the chromosomes for each parent, drops Mt+Pt, strips the parent prefix
(`B73_chr1` → `chr1`), and writes per-parent fastas to `1_split/genomes/`.

For B73 (which appears in both crosses), the script verifies byte-identity of
the sequences across the two merged fastas via per-chromosome MD5, then writes
a single shared `B73.fa`. Mismatch is a hard error.

Usage (run from project root):
    python 01a_split_genomes.py
    python 01a_split_genomes.py --verify-only
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import pysam
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"


def load_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def md5_hex(seq: str) -> str:
    return hashlib.md5(seq.encode("ascii")).hexdigest()


def write_fasta(out_path: Path, records: list[tuple[str, str]], line_width: int = 80) -> None:
    """Write a multi-record fasta with consistent line wrapping."""
    with open(out_path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i + line_width] + "\n")


def extract_per_parent(
    merged_fa: Path, parent_prefix: str, drop_chroms: set[str]
) -> tuple[list[tuple[str, str]], dict[str, str], dict[str, int]]:
    """Extract chromosomes for one parent from a merged fasta.

    Returns:
        records:  list of (renamed_name, sequence) tuples, in fasta order
        md5s:     dict of {renamed_name: md5_hex_of_sequence}
        lengths:  dict of {renamed_name: bp_length}
    """
    records: list[tuple[str, str]] = []
    md5s: dict[str, str] = {}
    lengths: dict[str, int] = {}

    fa = pysam.FastaFile(str(merged_fa))
    try:
        for ref in fa.references:
            if ref in drop_chroms:
                continue
            if not ref.startswith(parent_prefix):
                continue
            seq = fa.fetch(ref).upper()
            new_name = ref[len(parent_prefix):]  # strip "B73_" or "Ki3_" or "Oh43_"
            records.append((new_name, seq))
            md5s[new_name] = md5_hex(seq)
            lengths[new_name] = len(seq)
    finally:
        fa.close()

    return records, md5s, lengths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Only verify B73 byte-identity across crosses; do not write output fastas",
    )
    args = parser.parse_args()

    config = load_config()
    drop_chroms = set(config["split"]["drop_chromosomes"])
    out_dir = PROJECT_ROOT / config["paths"]["split_genomes"]
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] project root:      {PROJECT_ROOT}")
    print(f"[INFO] output directory:  {out_dir.relative_to(PROJECT_ROOT)}")
    print(f"[INFO] dropping chroms:   {sorted(drop_chroms)}")

    # Track B73 results across crosses for byte-identity check
    b73_records_per_cross: dict[str, list[tuple[str, str]]] = {}
    b73_md5_per_cross: dict[str, dict[str, str]] = {}
    b73_lengths_per_cross: dict[str, dict[str, int]] = {}

    # Non-B73 parents are unique per cross — collect them in order
    parent2_results: list[tuple[str, list[tuple[str, str]], dict[str, int]]] = []

    for cross, cross_cfg in config["crosses"].items():
        merged_fa = PROJECT_ROOT / cross_cfg["merged_fasta"]
        if not merged_fa.exists():
            print(f"[ERROR] merged fasta not found: {merged_fa}", file=sys.stderr)
            return 1

        print(f"\n[STEP] [{cross}] reading {merged_fa.relative_to(PROJECT_ROOT)}")

        for parent_role, parent_name in cross_cfg["parents"].items():
            prefix = f"{parent_name}_"
            records, md5s, lengths = extract_per_parent(merged_fa, prefix, drop_chroms)
            total_bp = sum(lengths.values())
            print(f"[INFO]   [{cross}] {parent_name}: {len(records)} sequences, {total_bp:>14,} bp")
            for name in sorted(lengths):
                print(f"[INFO]     {name:8s}  {lengths[name]:>12,} bp  md5={md5s[name][:12]}…")

            if parent_name == "B73":
                b73_records_per_cross[cross] = records
                b73_md5_per_cross[cross] = md5s
                b73_lengths_per_cross[cross] = lengths
            else:
                parent2_results.append((parent_name, records, lengths))

    # Verify B73 byte-identity across crosses
    print("\n[STEP] verifying B73 byte-identity across crosses")
    crosses = list(b73_md5_per_cross.keys())
    if len(crosses) < 2:
        print("[WARN] only one cross provides B73 — skipping cross-cross verification")
    else:
        ref_cross = crosses[0]
        ref_md5 = b73_md5_per_cross[ref_cross]
        ref_len = b73_lengths_per_cross[ref_cross]
        all_ok = True
        for other_cross in crosses[1:]:
            other_md5 = b73_md5_per_cross[other_cross]
            other_len = b73_lengths_per_cross[other_cross]
            if set(ref_md5.keys()) != set(other_md5.keys()):
                print(f"[ERROR] B73 chromosome SET differs between {ref_cross} and {other_cross}", file=sys.stderr)
                print(f"  {ref_cross}: {sorted(ref_md5.keys())}", file=sys.stderr)
                print(f"  {other_cross}: {sorted(other_md5.keys())}", file=sys.stderr)
                all_ok = False
                continue
            for chrom in sorted(ref_md5):
                if ref_len[chrom] != other_len[chrom]:
                    print(
                        f"[ERROR] {chrom} length differs: "
                        f"{ref_cross}={ref_len[chrom]:,} vs {other_cross}={other_len[chrom]:,}",
                        file=sys.stderr,
                    )
                    all_ok = False
                elif ref_md5[chrom] != other_md5[chrom]:
                    print(
                        f"[ERROR] {chrom} md5 differs: "
                        f"{ref_cross}={ref_md5[chrom]} vs {other_cross}={other_md5[chrom]}",
                        file=sys.stderr,
                    )
                    all_ok = False
            if all_ok:
                print(f"[OK]   B73 byte-identical between {ref_cross} and {other_cross}")
        if not all_ok:
            print("[ERROR] B73 sequences differ between crosses — aborting", file=sys.stderr)
            return 1

    if args.verify_only:
        print("\n[INFO] --verify-only mode: not writing output fastas")
        return 0

    # Write per-parent fastas
    print("\n[STEP] writing per-parent fastas")

    # B73 — use the first cross's records (already verified identical)
    b73_records = b73_records_per_cross[crosses[0]]
    b73_path = out_dir / "B73.fa"
    print(f"[STEP] writing {b73_path.relative_to(PROJECT_ROOT)} ({len(b73_records)} sequences)")
    write_fasta(b73_path, b73_records)
    print(f"[STEP] indexing {b73_path.name}")
    pysam.faidx(str(b73_path))

    # Other parents
    for parent_name, records, lengths in parent2_results:
        out_path = out_dir / f"{parent_name}.fa"
        total_bp = sum(lengths.values())
        print(f"[STEP] writing {out_path.relative_to(PROJECT_ROOT)} ({len(records)} sequences, {total_bp:,} bp)")
        write_fasta(out_path, records)
        print(f"[STEP] indexing {out_path.name}")
        pysam.faidx(str(out_path))

    print("\n[DONE] all per-parent fastas written and indexed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
