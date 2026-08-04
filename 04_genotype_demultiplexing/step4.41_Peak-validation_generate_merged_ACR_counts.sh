#!/bin/bash
#SBATCH --job-name=run_multicov
#SBATCH --output=logs/run_multicov_%A_%a.out
#SBATCH --error=logs/run_multicov_%A_%a.err
#SBATCH --time=05:00:00
#SBATCH --mem=6gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard

cd 4_genotyping_QC/4_peak_validation

ACR_bed="all_ACR.bed"
ACR_count="all_ACR_counts.txt"
bam_files=(*.bam) 

# ====== 1.Generate the union of all narrowPeak files (merged overlapping or adjacent regions)
cat *.narrowPeak | sort -k1,1 -k2,2n | bedtools merge > ${ACR_bed}

# ====== 2.Construct header line: chr start end sample1 sample2 ...
header="chr\tstart\tend"

for bam in "${bam_files[@]}"; do
    sample_name=$(basename "$bam" .bam)
    header="${header}\t${sample_name}"
done

# ====== 3.Count reads from each BAM file over the union regions
(echo -e "$header"; bedtools multicov -bams "${bam_files[@]}" -bed "${ACR_bed}") > "${ACR_count}"

