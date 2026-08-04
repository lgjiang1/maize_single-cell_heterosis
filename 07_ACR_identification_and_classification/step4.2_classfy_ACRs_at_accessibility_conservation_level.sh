#!/bin/bash

DIR="6_ACR_category/B73-Ki3-BxK-KxB/2_cross_genome_analysis/1_output"
cd $DIR

# Haplotype names
haps=("B73" "Ki3")

#############################################################################
# STEP1) Merge the three non-ambiguous ACR structure classes per haplotype
#############################################################################
for h in "${haps[@]}"; do
    files=(1_${h}_Conserved-ACR.tsv 2_${h}_PAV-ACR.tsv 3_${h}_Insertion-shift-ACR.tsv)
    outfile="${h}_Conserved_PAV_IS_merged.tsv"

    {
        head -n 1 "${files[0]}"
        for f in "${files[@]}"; do
            tail -n +2 "$f"
        done | sort -k1,1 -k2,2n
    } > "$outfile"

done


#############################################################################
# STEP2) Classify ACRs as shared / variable / PAV
#############################################################################
h1="${haps[0]}"
h2="${haps[1]}"
f1="${h1}_Conserved_PAV_IS_merged.tsv"
f2="${h2}_Conserved_PAV_IS_merged.tsv"

# -------- h1 view: use h1 coords (cols 1-3), compare to h1 coords in h2 file (cols 6-8) --------
awk 'NR>1{print $1,$2,$3,NR-1}' OFS="\t" "$f1" > ref.bed
awk 'NR>1 && $6!="."{print $6,$7,$8}' OFS="\t" "$f2" > target.bed

bedtools intersect -a ref.bed -b target.bed -wa | cut -f4 | sort -u > shared.ids

awk -v OFS="\t" '
    NR==FNR { shared[$1]=1; next }        # first file: shared.ids
    FNR==1  { print $0, "peak_level"; next }  # header line of merged file
    {
        id = FNR-1
        if($10=="PAV")            lab="PAV";
        else if(id in shared)     lab="shared";
        else                      lab="variable";
        print $0, lab
    }
' shared.ids "$f1" > Final_${h1}_Conserved_PAV_IS_with_peak_conservation_status.tsv


# -------- h2 view: use h2 coords (cols 6-8), compare to h2 coords in h1 file (cols 1-3) --------
awk 'NR>1 && $6!="."{print $6,$7,$8,NR-1}' OFS="\t" "$f2" > ref.bed
awk 'NR>1{print $1,$2,$3}' OFS="\t" "$f1" > target.bed

bedtools intersect -a ref.bed -b target.bed -wa | cut -f4 | sort -u > shared.ids

awk -v OFS="\t" '
    NR==FNR { shared[$1]=1; next }        # first file: shared.ids
    FNR==1  { print $0, "peak_level"; next }  # header line of merged file
    {
        id = FNR-1
        if($10=="PAV")            lab="PAV";
        else if(id in shared)     lab="shared";
        else                      lab="variable";
        print $0, lab
    }
' shared.ids "$f2" > Final_${h2}_Conserved_PAV_IS_with_peak_conservation_status.tsv

rm -f ref.bed target.bed shared.ids
