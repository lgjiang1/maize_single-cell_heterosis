#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
import numpy as np
import scipy.io

# =================== Parse command-line argument =======================
parser = argparse.ArgumentParser(description="Compute allele ratio matrix from vartrix output.")
parser.add_argument('--input_dir', required=True, help='Directory with alt.mtx, ref.mtx, barcodes.tsv, variants.tsv')
parser.add_argument('--sample_name', required=True, help='Sample name used as output prefix')
args = parser.parse_args()

# =================== Data Path =======================
alt_path = os.path.join(args.input_dir, "alt.mtx")
ref_path = os.path.join(args.input_dir, "ref.mtx")
barcode_path = os.path.join(args.input_dir, "barcode.tsv")
variant_path = os.path.join(args.input_dir, "variants.tsv")

# ======================= Load matrices ===========================
alt = scipy.io.mmread(alt_path).tocoo().tocsr()
ref = scipy.io.mmread(ref_path).tocoo().tocsr()
total = alt + ref
dense_total = total.toarray()

barcodes = pd.read_csv(barcode_path, header=None)[0].tolist()
snp_ids = pd.read_csv(variant_path, header=None)[0].tolist()

# ===================== Compute allele ratio ======================
allele_ratio = alt.multiply(1 / (dense_total + 1e-6)).toarray()
allele_ratio[dense_total == 0] = np.nan

# ===================== Write output ======================
df = pd.DataFrame(allele_ratio, index=snp_ids, columns=barcodes)
output_path = os.path.join(args.input_dir, f"{args.sample_name}_allele_ratio.csv")
df.to_csv(output_path)

print(f"[DONE] Wrote allele ratio matrix to {output_path}")
print(f"[INFO] Matrix shape: {df.shape[0]} SNPs × {df.shape[1]} cells")
