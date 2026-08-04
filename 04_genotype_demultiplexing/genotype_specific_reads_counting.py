#!/usr/bin/env python

import pysam
import sys
from collections import defaultdict

# === Input arguments ===
bam_file = sys.argv[1]        # Input BAM file (aligned reads with barcode tag "BC")
vcf_file = sys.argv[2]        # Input VCF file containing SNPs for one parent
output_file = sys.argv[3]     # Output file: barcode read stats

# === Initialize data structure ===
# barcode_counts[barcode] = [matching_reads, total_reads]
barcode_counts = defaultdict(lambda: [0, 0])

# === Open BAM file ===
bam = pysam.AlignmentFile(bam_file, "rb")

# === Parse VCF and process each SNP ===
with open(vcf_file) as vcf:
    for line in vcf:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        chrom = fields[0]
        pos = int(fields[1]) - 1  # Convert to 0-based for pileup
        ref = fields[3]
        alt = fields[4]
        samples = fields[9:]

        if len(samples) < 1:
            continue  # No genotype info

        gt = samples[0].split(":")[0]

        # Convert GT to allele (must be homozygous)
        def get_allele(gt):
            if gt == "0/0":
                return ref
            elif gt == "1/1":
                return alt
            else:
                return None  # Skip heterozygous or missing

        parent_allele = get_allele(gt)
        if parent_allele is None:
            continue

        # Pileup at SNP
        for pileupcolumn in bam.pileup(chrom, pos, pos + 1, truncate=True, stepper="samtools"):
            if pileupcolumn.pos != pos:
                continue
            for pileupread in pileupcolumn.pileups:
                if pileupread.is_del or pileupread.is_refskip:
                    continue
                read = pileupread.alignment
                if not read.has_tag("BC"):
                    continue
                barcode = read.get_tag("BC")
                base = read.query_sequence[pileupread.query_position]
                barcode_counts[barcode][1] += 1  # total_reads
                if base == parent_allele:
                    barcode_counts[barcode][0] += 1  # matching_reads

# === Write result ===
with open(output_file, "w") as out:
    out.write("barcode\tmatching_reads\ttotal_reads\tfraction\n")
    for bc in sorted(barcode_counts):
        matching, total = barcode_counts[bc]
        frac = matching / total if total > 0 else 0
        out.write(f"{bc}\t{matching}\t{total}\t{frac:.3f}\n")

### Barcodes exhibiting >5% matching to the non-target parent were flagged as putative contaminants and excluded.
