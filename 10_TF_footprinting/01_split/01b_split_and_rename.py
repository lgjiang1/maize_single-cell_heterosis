"""
01b_split_and_rename.py — split a single BAM by parent + rename chromosomes.

Generic pysam-based BAM rewriter. Used both for parental BAMs (1 invocation per
input BAM, target_parent matches the parent's own name) and for F1 BAMs (2
invocations per input BAM, one per allele).

Behavior (one input → one output):
    1. Read input BAM
    2. Build a new header keeping only @SQ entries that:
        (a) start with the target_parent prefix, e.g. "B73_"
        (b) are not in the drop list (Mt, Pt by default)
       Kept SQs are renamed by stripping the prefix (B73_chr1 → chr1).
    3. Stream reads, applying these filters in order:
        - drop unmapped
        - drop secondary / supplementary
        - drop reads on chromosomes outside the target set (other parent, Mt, Pt)
        - drop reads whose mate is on a chromosome outside the target set
        - drop reads with MAPQ < --mapq
    4. Rewrite RNAME and RNEXT to the new chromosome names
    5. Write the output BAM

Per-BAM stats are written to {output}.stats.json next to the output BAM.

Usage (single BAM, runnable directly or from a SLURM array task):
    python 01b_split_and_rename.py \\
        --input  0_bams/0_B-K_input/B73-C1.renamed.bam \\
        --output 1_split/B-K/B73_coords/B73-C1.bam \\
        --target-parent B73
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pysam
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"


def load_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def build_new_header(
    old_header_dict: dict, target_prefix: str, drop_chroms: set[str]
) -> tuple[dict, dict[str, str]]:
    """Build the new header keeping only target-parent SQs (renamed).

    Returns:
        new_header_dict: pass to pysam.AlignmentFile(..., header=new_header_dict)
        old_to_new_name: dict mapping old_full_name → new_short_name
    """
    new_sq: list[dict] = []
    old_to_new_name: dict[str, str] = {}

    for sq in old_header_dict.get("SQ", []):
        name = sq["SN"]
        if name in drop_chroms:
            continue
        if not name.startswith(target_prefix):
            continue
        new_name = name[len(target_prefix):]
        new_sq.append({"SN": new_name, "LN": sq["LN"]})
        old_to_new_name[name] = new_name

    new_header: dict = {
        "HD": old_header_dict.get("HD", {"VN": "1.6", "SO": "coordinate"}),
        "SQ": new_sq,
    }
    # Pass through @PG / @RG / @CO so we keep upstream provenance
    for key in ("PG", "RG", "CO"):
        if key in old_header_dict:
            new_header[key] = old_header_dict[key]

    return new_header, old_to_new_name


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="Input BAM (path relative to project root or absolute)")
    parser.add_argument("--output", required=True, help="Output BAM")
    parser.add_argument("--target-parent", required=True, choices=["B73", "Ki3", "Oh43"])
    parser.add_argument("--mapq", type=int, default=None, help="MAPQ threshold (overrides config)")
    parser.add_argument("--qc-mapq-threshold", type=int, default=None, help="Higher MAPQ tracked for QC (overrides config)")
    args = parser.parse_args()

    config = load_config()
    drop_chroms = set(config["split"]["drop_chromosomes"])
    mapq_filter = args.mapq if args.mapq is not None else int(config["split"]["mapq_filter"])
    qc_mapq = args.qc_mapq_threshold if args.qc_mapq_threshold is not None else int(config["split"]["qc_mapq_threshold"])
    qc_max_q30_loss = float(config["split"]["qc_max_q30_loss"])
    target_prefix = f"{args.target_parent}_"

    # Resolve paths (allow both relative-to-project-root and absolute)
    def resolve(p: str) -> Path:
        path = Path(p)
        if not path.is_absolute():
            path = PROJECT_ROOT / path
        return path

    input_bam = resolve(args.input)
    output_bam = resolve(args.output)

    if not input_bam.exists():
        print(f"[ERROR] input BAM not found: {input_bam}", file=sys.stderr)
        return 1

    output_bam.parent.mkdir(parents=True, exist_ok=True)
    stats_path = output_bam.with_suffix(output_bam.suffix + ".stats.json")

    print(f"[INFO] input:           {input_bam}")
    print(f"[INFO] output:          {output_bam}")
    print(f"[INFO] target_parent:   {args.target_parent} (prefix={target_prefix})")
    print(f"[INFO] mapq filter:     >= {mapq_filter}")
    print(f"[INFO] qc mapq:         >= {qc_mapq}")
    print(f"[INFO] drop chroms:     {sorted(drop_chroms)}")

    in_bam = pysam.AlignmentFile(str(input_bam), "rb")
    old_header_dict = in_bam.header.to_dict()

    new_header, old_to_new_name = build_new_header(old_header_dict, target_prefix, drop_chroms)
    if not new_header["SQ"]:
        print(f"[ERROR] new header has no SQ entries — check target_parent prefix '{target_prefix}'", file=sys.stderr)
        in_bam.close()
        return 1

    print(f"[INFO] kept {len(new_header['SQ'])} SQ entries (renamed)")

    out_bam = pysam.AlignmentFile(str(output_bam), "wb", header=new_header)

    # Counters
    total_reads = 0
    mapped_reads = 0
    primary_reads = 0
    target_parent_reads = 0
    mate_ok_reads = 0
    mapq10_pass = 0
    mapq30_pass = 0
    output_reads = 0

    for read in in_bam:
        total_reads += 1

        if read.is_unmapped:
            continue
        mapped_reads += 1

        if read.is_secondary or read.is_supplementary:
            continue
        primary_reads += 1

        ref_name = read.reference_name
        if ref_name is None:
            continue
        if ref_name in drop_chroms or not ref_name.startswith(target_prefix):
            continue
        target_parent_reads += 1

        # Mate check (proper-pair filter was applied upstream, so this should
        # almost always pass; here we just defend against the rare case where
        # the mate landed on a chromosome we're dropping)
        if read.is_paired and not read.mate_is_unmapped:
            mate_ref = read.next_reference_name
            if mate_ref is None:
                continue
            if mate_ref != "=" and (mate_ref in drop_chroms or not mate_ref.startswith(target_prefix)):
                continue
        mate_ok_reads += 1

        # MAPQ counters (computed on the post-mate-check pool, so the q30
        # loss fraction is well-defined)
        passes_q10 = read.mapping_quality >= mapq_filter
        passes_q30 = read.mapping_quality >= qc_mapq
        if passes_q10:
            mapq10_pass += 1
        if passes_q30:
            mapq30_pass += 1

        if not passes_q10:
            continue

        # Rewrite reference fields via the SAM-dict round-trip
        d = read.to_dict()
        old_ref = d["ref_name"]
        d["ref_name"] = old_to_new_name[old_ref]
        old_next = d.get("next_ref_name", "*")
        if old_next not in ("=", "*"):
            d["next_ref_name"] = old_to_new_name[old_next]

        new_read = pysam.AlignedSegment.from_dict(d, out_bam.header)
        out_bam.write(new_read)
        output_reads += 1

    out_bam.close()
    in_bam.close()

    # Stats
    if mapq10_pass > 0:
        q30_loss_fraction = 1.0 - (mapq30_pass / mapq10_pass)
    else:
        q30_loss_fraction = None

    stats = {
        "input_bam": str(input_bam),
        "output_bam": str(output_bam),
        "target_parent": args.target_parent,
        "mapq_threshold": mapq_filter,
        "qc_mapq_threshold": qc_mapq,
        "total_reads": total_reads,
        "mapped_reads": mapped_reads,
        "primary_reads": primary_reads,
        "target_parent_reads": target_parent_reads,
        "mate_ok_reads": mate_ok_reads,
        "mapq10_pass": mapq10_pass,
        "mapq30_pass": mapq30_pass,
        "q30_loss_fraction": q30_loss_fraction,
        "output_reads": output_reads,
    }
    with open(stats_path, "w") as f:
        json.dump(stats, f, indent=2)

    print(f"[OK]   wrote {output_reads:,} reads -> {output_bam.relative_to(PROJECT_ROOT)}")
    print(f"[OK]   stats -> {stats_path.relative_to(PROJECT_ROOT)}")
    print(
        f"[INFO] total={total_reads:,} mapped={mapped_reads:,} primary={primary_reads:,} "
        f"target_parent={target_parent_reads:,} mate_ok={mate_ok_reads:,} "
        f"mapq{mapq_filter}={mapq10_pass:,} output={output_reads:,}"
    )
    if q30_loss_fraction is not None:
        flag = "  [FLAG: > %.0f%%]" % (qc_max_q30_loss * 100) if q30_loss_fraction > qc_max_q30_loss else ""
        print(f"[INFO] q30_loss_fraction = {q30_loss_fraction:.4f}{flag}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
