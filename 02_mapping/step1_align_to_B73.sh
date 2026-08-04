#!/bin/bash
#SBATCH --job-name=bwa_align_sort
#SBATCH --output=logs/bwa_align_sort_%A_%a.out
#SBATCH --error=logs/bwa_align_sort_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=15
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --array=1-14           ######Adjust the range based on the number of read pairs

ml Bioinformatics bwa

input_dir="02_clean_reads/Final_reads_for_mapping"
output_dir="03_mapping/1_raw_bam"
reference="../../../reference_genomes/zea/B73/B73v5_bwa"
stats_dir="03_mapping/1_raw_bam/alignment_stats"

mkdir -p $output_dir
mkdir -p $stats_dir

# get read pairs (R1 and R3 read files)
file_pairs=($(ls ${input_dir}/*_R1.fastq.gz | sed 's/_R1.fastq.gz//' | sort))
pair=${file_pairs[$SLURM_ARRAY_TASK_ID-1]}

# input_file
r1_file=${pair}_R1.fastq.gz
r3_file=${pair}_R3.fastq.gz

# output file prefixes
sample_name=$(basename ${pair})
sam_file=${output_dir}/${sample_name}.sam
bam_file=${output_dir}/${sample_name}.sorted.bam

#===================================================================================
#===================          Main Command           ==============================
#===================================================================================

#1) bwa mapping
bwa mem -M -t 15 ${reference} ${r1_file} ${r3_file} > ${sam_file}

#2) convert SAM to BAM and sort
samtools view -@ 15 -bS ${sam_file} | samtools sort -@ 15 -o ${bam_file}

#3) indexing
samtools index -@ 15 ${bam_file}

#4) generate bw file
bamCoverage --bam ${bam_file} \
            --outFileName ${bam_file}.bw \
            --outFileFormat bigwig \
            --binSize 10 \
            --numberOfProcessors 15 \
            --normalizeUsing CPM

#5) summarize alignment information in bam files
samtools flagstat -@ 15 ${bam_file} > ${stats_dir}/alignment_stats_${sample_name}.txt

#6) remove intermediate SAM file
rm ${sam_file}

echo "Alignment and BAM file generation completed for sample ${sample_name}"

#====================================================================================
