#!/bin/bash
#SBATCH --job-name=subset-bam
#SBATCH --output=logs/subset-bam_%A_%a.out
#SBATCH --error=logs/subset-bam_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=15gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --array=1-18

bam_dir="4_genotyping_QC/0_bam_file"
barcode_dir="4_genotyping_QC/0_high_quality_and_genotyped_barcodes"
subset_bam_dir="4_genotyping_QC/1_subset_bam_file"

mkdir -p ${subset_bam_dir}

#=========== step1: parse input line to get BAM path, barcode path, and sample ID components ===================
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${subset_bam_dir}/subset_input_list.txt)
bam_file=$(echo "$line" | cut -f1)
bc_file=$(echo "$line" | cut -f2)

bam="${bam_dir}/${bam_file}"
bc="${barcode_dir}/${bc_file}"

basename_noext="${bc_file%_barcode.txt}"
prefix=${basename_noext%%_*}
geno=${basename_noext#*_}


#================== step2: run subset-bam, add read group for each genotype and index ========================
raw_bam="${subset_bam_dir}/${basename_noext}_subset.bam"
rg_bam="${subset_bam_dir}/${basename_noext}_subset.rg.bam"

# Extract reads by barcode
echo "[INFO] Running subset-bam on ${bam} with ${bc}"
subset-bam \
    --bam $bam \
    --bam-tag BC \
    --cell-barcodes $bc \
    --cores 4 \
    --out-bam $raw_bam 

# Add read group
echo "[INFO] Adding read group to $raw_bam"
samtools addreplacerg \
    -r "@RG\tID:${prefix}_${geno}\tSM:${prefix}_${geno}\tPL:ILLUMINA" \
    -o $rg_bam \
    -@ 4 \
    $raw_bam

# Build BAM index
echo "[INFO] Indexing $rg_bam"
samtools index $rg_bam

rm $raw_bam

