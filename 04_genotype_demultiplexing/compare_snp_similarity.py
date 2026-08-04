#!/usr/bin/env python

import pandas as pd
import numpy as np
import gzip
import argparse

# -------- Get all sample names from VCF header --------
def get_sample_names_from_vcf(vcf_file):
    opener = gzip.open if vcf_file.endswith(".gz") else open
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith("#CHROM"):
                header = line.strip().split('\t')
                return header[9:]

# -------- Read GT values from VCF file for all samples --------
def read_vcf_genotypes(vcf_file, sample_names):
    opener = gzip.open if vcf_file.endswith(".gz") else open
    data = {name: [] for name in sample_names}
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            for name, idx in zip(sample_names, range(9, 9 + len(sample_names))):
                gt = fields[idx].split(":")[0]
                data[name].append("NA" if gt == "./." else gt)
    return pd.DataFrame(data)

# -------- Compute similarity matrix per cell --------
def compute_similarity_matrix(tsv_df, vcf_df, min_valid_snp=100):
    results = {}
    for cell in tsv_df.columns:
        cell_vals = tsv_df[cell].values

        ratios = []
        for sample in vcf_df.columns:
            sample_vals = vcf_df[sample].values
            valid_mask = (cell_vals != "NA") & (sample_vals != "NA")
            total = np.sum(valid_mask)
            if total == 0 or total < min_valid_snp:
                ratios.append(np.nan)
            else:
                match = np.sum(cell_vals[valid_mask] == sample_vals[valid_mask])
                ratios.append(match / total)
        results[cell] = ratios

    return pd.DataFrame(results, index=vcf_df.columns).T

# -------- Process all cells in column-wise chunks --------
def process_in_chunks(tsv_file, vcf_df, chunk_size=1000, output_file="output.tsv", min_valid_snp=100):
    print(">> Reading full SNP × Cell matrix...")
    full_df = pd.read_csv(tsv_file, sep="\t", index_col=0, dtype=str)
    full_df = full_df.fillna("NA")

    all_cells = full_df.columns
    num_cells = len(all_cells)
    first = True

    for start in range(0, num_cells, chunk_size):
        end = min(start + chunk_size, num_cells)
        cell_chunk = all_cells[start:end]
        chunk = full_df[cell_chunk]

        print(f">> Processing cell columns {start + 1} to {end} of {num_cells} ...")
        result = compute_similarity_matrix(chunk, vcf_df, min_valid_snp=min_valid_snp).round(2)
        result.to_csv(output_file, sep="\t", mode='w' if first else 'a', header=first)
        first = False

# -------- Main function --------
def main():
    parser = argparse.ArgumentParser(description="Compare SNPs of single cells to VCF genotypes.")
    parser.add_argument("--vcf_file", required=True, help="Input VCF file (.vcf or .vcf.gz)")
    parser.add_argument("--tsv_file", required=True, help="Input SNP × Cell matrix TSV file")
    parser.add_argument("--output_file", required=True, help="Output file for similarity results")
    parser.add_argument("--chunk_size", type=int, default=1000, help="Number of cells to process per chunk (default: 1000)")
    parser.add_argument("--min_snps", type=int, default=100, help="Minimum number of valid SNPs per cell to include (default: 100)")
    args = parser.parse_args()

    print(">> Extracting sample names from VCF...")
    sample_names = get_sample_names_from_vcf(args.vcf_file)
    print(f">> Found {len(sample_names)} samples: {sample_names}")

    print(">> Reading VCF genotypes...")
    vcf_df = read_vcf_genotypes(args.vcf_file, sample_names)
    vcf_snp_count = vcf_df.shape[0]
    print(f">> SNP count in VCF: {vcf_snp_count}")

    # Count SNPs in TSV (lines - 1 header)
    with open(args.tsv_file) as f:
        tsv_lines = sum(1 for _ in f)
    tsv_snp_count = tsv_lines - 1
    print(f">> SNP count in TSV: {tsv_snp_count}")

    if vcf_snp_count != tsv_snp_count:
        print(f"[ERROR] SNP row count mismatch: VCF has {vcf_snp_count}, TSV has {tsv_snp_count}")
        exit(1)

    print(">> Processing cells in column-wise chunks...")
    process_in_chunks(
        tsv_file=args.tsv_file,
        vcf_df=vcf_df,
        chunk_size=args.chunk_size,
        output_file=args.output_file,
        min_valid_snp=args.min_snps
    )

    print(f">> Done! Similarity matrix saved to: {args.output_file}")

if __name__ == "__main__":
    main()

