#!/bin/bash
#SBATCH --job-name=peakcalling
#SBATCH --output=logs/peakcalling_%A_%a.out
#SBATCH --error=logs/peakcalling_%A_%a.err
#SBATCH --time=1-15:00:00
#SBATCH --mem=25gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --array=1-30

ml Bioinformatics py-macs2

#======================= step1: make Tn5 bed file ========================
bam_dir="4_genotyping_QC/2_peak_calling/bam_files/"
bam_list=($(ls $bam_dir/*.bam))
bam=${bam_list[$((SLURM_ARRAY_TASK_ID - 1))]}
basename="${bam##*/}"; basename="${basename%%.*}"

tn5_bed=${bam_dir}/${basename}.tn5.bed

perl ./makeTn5bed.pl --bam $bam | sort -k1,1 -k2,2n - | uniq - > $tn5_bed

gzip $tn5_bed


#========================= step2: peak calling ===========================
peak_dir="4_genotyping_QC/2_peak_calling/macs2_result"
mkdir -p $peak_dir

macs2 callpeak -t $tn5_bed.gz -f BED -g 1.6e9 -n $basename \
               --nomodel --extsize 150 --shift -75 \
               --keep-dup all --qvalue 0.05 --outdir $peak_dir/$basename

