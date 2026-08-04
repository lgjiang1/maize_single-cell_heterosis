#!/bin/bash
#SBATCH --job-name=append_bcs
#SBATCH --output=logs/append_bcs_%A_%a.out
#SBATCH --error=logs/append_bcs_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=1gb
#SBATCH --time=4-00:00:00
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --array=1-2        ###Adjust the range based on the number of files

##This script is used to add 10x cell barcodes and Tn5 barcodes to each read name.

dir="01_raw_fastqs"
outdir="02_clean_reads"

mkdir -p $outdir

# retrieve file prefixes
prefixes=($(ls ${dir}/*_R1_001.fastq.gz | sed 's|.*/||' | sed 's|_R1_001.fastq.gz||'))
prefix=${prefixes[$SLURM_ARRAY_TASK_ID-1]}

# input_file
tenxBC=${dir}/${prefix}_R2_001.fastq.gz
R1=${dir}/${prefix}_R1_001.fastq.gz
R2=${dir}/${prefix}_R3_001.fastq.gz

# output_file
R1tenx=${outdir}/${prefix}_R1_001.bc1.fastq.gz
R2tenx=${outdir}/${prefix}_R3_001.bc1.fastq.gz
R1tenxtn5=${outdir}/${prefix}_R1_001.bc1.bc2.fastq.gz
R2tenxtn5=${outdir}/${prefix}_R3_001.bc1.bc2.fastq.gz


#==================== append 10x barcodes to read_name ==========================
# append 10x barcode to R1 read-header
umi_tools extract \
    --bc-pattern=NNNNNNNNNNNNNNNN \
    --stdin=$tenxBC \
    --read2-in=$R1 \
    --stdout=$R1tenx \
    --read2-stdout

# append 10x barcode to R3 read-header
umi_tools extract \
    --bc-pattern=NNNNNNNNNNNNNNNN \
    --stdin=$tenxBC \
    --read2-in=$R2 \
    --stdout=$R2tenx \
    --read2-stdout


#=========== append Tn5 barcode from R1 and R3 to read-header ====================
cutadapt -e 0.2 \
         --pair-filter=any \
         -j 10 \
         --rename='{id}_{r1.cut_prefix}_{r2.cut_prefix} {comment}' \
         -u 5 \
         -U 5 \
         -g AGATGTGTATAAGAGACAG \
         -G AGATGTGTATAAGAGACAG \
         -o $R1tenxtn5 \
         -p $R2tenxtn5 \
         $R1tenx $R2tenx

