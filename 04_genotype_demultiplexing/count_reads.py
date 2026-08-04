#!/usr/bin/env python

import pysam
import sys
from collections import defaultdict

# === Input arguments ===
bam_file = sys.argv[1]        # Input BAM file (aligned reads with barcode tag "BC")
vcf_file = sys.argv[2]        # Input VCF file containing SNPs for two parents
output_file = sys.argv[3]     # Output file: barcode read counts for each parent

# === Initialize data structure ===
# For each barcode: [P1_read_count, P2_read_count]
barcode_counts = defaultdict(lambda: [0, 0])

# === Open BAM file ===
bam = pysam.AlignmentFile(bam_file, "rb")

# === Parse VCF and process each informative SNP ===
with open(vcf_file) as vcf:
    for line in vcf:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        chrom = fields[0]
        pos = int(fields[1]) - 1
        ref = fields[3]
        alt = fields[4]
        format_fields = fields[8].split(":")
        samples = fields[9:]

        # Extract GTs for P1 and P2
        if len(samples) < 2:
            continue  # Not enough samples
        p1_gt = samples[0].split(":")[0]
        p2_gt = samples[1].split(":")[0]

        # Skip if genotypes are the same (no informative difference)
        if p1_gt == p2_gt:
            continue

        # Function to convert GT to allele
        def get_allele(gt):
            if gt == "0/0":
                return ref
            elif gt == "1/1":
                return alt
            elif gt == "0/1":
                return None  # Ambiguous heterozygous
            else:
                return None

        p1_allele = get_allele(p1_gt)
        p2_allele = get_allele(p2_gt)

        # Require both parents to be homozygous and different
        if p1_allele is None or p2_allele is None or p1_allele == p2_allele:
            continue

        # Pileup at this SNP
        for pileupcolumn in bam.pileup(chrom, pos, pos + 1, truncate=True, stepper="samtools"):
            if pileupcolumn.pos != pos:
                continue
            for pileupread in pileupcolumn.pileups:
                if pileupread.is_del or pileupread.is_refskip:
                    continue
                read = pileupread.alignment
                base = read.query_sequence[pileupread.query_position]
                if not read.has_tag("BC"):
                    continue
                barcode = read.get_tag("BC")

                if base == p1_allele:
                    barcode_counts[barcode][0] += 1
                if base == p2_allele:
                    barcode_counts[barcode][1] += 1

# === Write result ===
with open(output_file, "w") as out:
    out.write("barcode\tP1_reads\tP2_reads\n")
    for bc in sorted(barcode_counts):
        p1, p2 = barcode_counts[bc]
        out.write(f"{bc}\t{p1}\t{p2}\n")

