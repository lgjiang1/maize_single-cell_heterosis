#!/bin/bash
#SBATCH --job-name=process_bam
#SBATCH --output=logs/process_bam_%A_%a.out
#SBATCH --error=logs/process_bam_%A_%a.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=30gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --array=1-14          ######Adjust the range based on the number of read pairs

input_dir="03_mapping/1_raw_bam"
intermediate_dir="03_mapping/2_intermediate_bam"
output_dir="03_mapping/3_final_bam_and_tn5bed"

mkdir -p $intermediate_dir
mkdir -p $output_dir

input_files=($(ls $input_dir/*.sorted.bam))
input_file=${input_files[$SLURM_ARRAY_TASK_ID-1]}

basename=$(basename $input_file .sorted.bam)
mq10_bam=$intermediate_dir/${basename}.mq10.bam
barcodes_txt=$intermediate_dir/${basename}.mq10.barcodes.txt
barcodes_corrected_txt=$intermediate_dir/${basename}.mq10.barcodes.corrected.txt
bc_bam=$intermediate_dir/${basename}.mq10.BC.bam
bc_rmdup_bam=$intermediate_dir/${basename}.mq10.BC.rmdup.bam
bc_rmdup_mm_bam=$output_dir/${basename}.mq10.BC.rmdup.mm.bam
tn5_bed=$output_dir/${basename}.tn5.bed


#============================== Main Function ======================================
# Step 1: Modify BC flag and filter low MQ
perl ./modify_BC_flag.pl $input_file | samtools view -@ 15 -hbq 10 -f 3 - > $mq10_bam
echo ".. done bc tag and low mq for $basename .."

# Step 2: Get barcode counts
perl ./countBCs.BAM.pl $mq10_bam | awk '{ if ($2 > 49) { print }}' - > $barcodes_txt
echo ".. done barcode counts for $basename .."

# Step 3: Correct barcodes
cat $barcodes_txt | parallel --pipe -k -j 30 -N 1000 perl ./correctBCs.10x.v2.pl > $barcodes_corrected_txt
echo ".. done barcode correction for $basename .."

# Step 4: Update barcodes
perl ./correctBAM.pl $barcodes_corrected_txt $mq10_bam | samtools view -@ 15 -bhS -f 3 - > $bc_bam
echo ".. done barcode update for $basename .."

# Step 5: Remove duplicates
picard MarkDuplicates \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    REMOVE_DUPLICATES=true \
    METRICS_FILE=$intermediate_dir/${basename}.metrics \
    I=$bc_bam \
    O=$bc_rmdup_bam \
    BARCODE_TAG=BC \
    ASSUME_SORT_ORDER=coordinate
echo ".. done deduplication for $basename .."

# Step 6: Fix barcodes and remove multi-mapped reads
perl ./fixBC.pl $bc_rmdup_bam $basename | samtools view -@ 15 -bhS - > $bc_rmdup_mm_bam
echo ".. done fixing barcodes for $basename .."

# Step 7: Make Tn5 bed files
perl ./makeTn5bed.pl --bam $bc_rmdup_mm_bam | sort -k1,1 -k2,2n - | uniq - > $tn5_bed
gzip $tn5_bed
echo ".. done making Tn5 bed file for $basename .."

