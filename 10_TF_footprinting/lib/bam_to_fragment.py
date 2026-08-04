#!/usr/bin/env python3
"""
bam_to_fragment.py — Extract paired-end fragments from a NAME-SORTED BAM.

Output: 0-based fragment TSV (chr, start, end, barcode), gzipped.
Adapted from 13_Arabidopsis_protoplast/5_TF_FP/final/lib/bam_to_fragment.py.

For our pseudobulk pipeline, all fragments get barcode="bulk" (one pseudo-cell
per condition). MAPQ filtering is skipped because it was already applied at
the BAM split step (Phase 1, q>=10).

Usage:
    samtools sort -n input.bam -o input.namesort.bam
    python bam_to_fragment.py --bam input.namesort.bam --out frags.tsv.gz
"""
import sys
import gzip
import argparse
from pathlib import Path
from collections import Counter

import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True, help="NAME-SORTED BAM (samtools sort -n)")
    ap.add_argument("--out", required=True, help="Output fragments.tsv.gz")
    ap.add_argument("--barcode", default="bulk", help="Barcode label for all fragments (default: bulk)")
    ap.add_argument("--min-mapq", type=int, default=0, help="MAPQ filter (default: 0, already filtered)")
    ap.add_argument("--threads", type=int, default=4, help="pysam reading threads")
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb", check_sq=False, threads=args.threads)
    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)

    drops = Counter()
    n_pairs = n_kept = 0
    cur_name = None
    first = None
    min_mapq = args.min_mapq
    bc = args.barcode

    def ok(r):
        return (
            not r.is_unmapped
            and not r.is_secondary
            and not r.is_supplementary
            and r.mapping_quality >= min_mapq
            and r.is_paired
        )

    with gzip.open(outp, "wt") as out:
        for rec in bam:
            if rec.query_name != cur_name:
                cur_name = rec.query_name
                first = rec
                continue

            second = rec
            n_pairs += 1

            if not ok(first) or not ok(second):
                drops["filtered"] += 1
                cur_name = None
                first = None
                continue

            if first.reference_id != second.reference_id:
                drops["different_chrom"] += 1
                cur_name = None
                first = None
                continue

            rname = bam.get_reference_name(first.reference_id)
            s1, e1 = first.reference_start, first.reference_end
            s2, e2 = second.reference_start, second.reference_end
            if None in (s1, e1, s2, e2):
                drops["no_coords"] += 1
                cur_name = None
                first = None
                continue

            start = min(s1, s2)
            end = max(e1, e2)
            if end <= start:
                drops["invalid_len"] += 1
                cur_name = None
                first = None
                continue

            out.write(f"{rname}\t{start}\t{end}\t{bc}\n")
            n_kept += 1
            cur_name = None
            first = None

    bam.close()

    sys.stderr.write(f"[STATS] pairs={n_pairs:,} kept={n_kept:,} dropped={sum(drops.values()):,}\n")
    for k, v in drops.most_common():
        sys.stderr.write(f"[DROP] {k}: {v:,}\n")
    sys.stderr.write(f"[OUT] {outp}\n")


if __name__ == "__main__":
    main()
