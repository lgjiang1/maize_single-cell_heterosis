#!/usr/bin/env python
import sys
import glob
import pandas as pd
from intervaltree import IntervalTree
import os

# === 1. Read chromosome name from command line (e.g., "chr1") ===
chrom = sys.argv[1]

# === 2. Recursively find all peak files under the target folder (for example B73, Ki3 or BxK, etc.) ===
files = glob.glob("5_ACR_union/B73/0_ACR_per_celltype_replicate/*.reproducible.bed")

# === 3. Merge all input files into a single DataFrame ===
all_peaks = pd.concat(
    [pd.read_csv(f, sep="\t", header=None, usecols=[0,1,2,3,4,5,6]) for f in files],
    ignore_index=True
)

# === 4. Assign column names ===
all_peaks.columns = ["chr", "start", "end", "name", "qValue", "summit", "source"]

# Ensure types are valid
all_peaks["start"]  = pd.to_numeric(all_peaks["start"],  errors="coerce").astype("Int64")
all_peaks["end"]    = pd.to_numeric(all_peaks["end"],    errors="coerce").astype("Int64")
all_peaks["qValue"] = pd.to_numeric(all_peaks["qValue"], errors="coerce")
all_peaks = all_peaks.dropna(subset=["start","end","qValue"])
all_peaks = all_peaks[all_peaks["start"] < all_peaks["end"]]

# === 5. Filter for the specified chromosome only ===
chr_peaks = all_peaks[all_peaks["chr"].astype(str) == str(chrom)].copy()

# === 6. Sort peaks by qValue descending (most significant first) ===
chr_peaks.sort_values(by="qValue", ascending=False, inplace=True)

# === 7. Remove overlapping peaks, keeping the most significant ===
selected = []
tree = IntervalTree()

for _, row in chr_peaks.iterrows():
    s, e = int(row["start"]), int(row["end"])
    if not tree.overlaps(s, e):
        selected.append([
            row["chr"], s, e, float(row["qValue"]), int(row["summit"]), row["source"]
        ])
        tree[s:e] = True  # mark region as used

# === 8. Save the output for this chromosome ===
os.makedirs("5_ACR_union/B73/1_union_peak_set/filtered_by_chr", exist_ok=True)
out_path = f"5_ACR_union/B73/1_union_peak_set/filtered_by_chr/non_overlapping_{chrom}.bed"
pd.DataFrame(selected, columns=["chr","start","end","qValue","summit","source"]) \
  .to_csv(out_path, sep="\t", index=False, header=False)


