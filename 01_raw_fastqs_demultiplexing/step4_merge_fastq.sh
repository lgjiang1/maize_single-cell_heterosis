#!/bin/bash
#SBATCH --job-name=fastq-merge
#SBATCH --output=logs/fastq-merge_%j.out
#SBATCH --error=logs/fastq-merge_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5gb
#SBATCH --time=01:00:00
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard

input_dir="02_clean_reads/Assign"
output_dir="02_clean_reads/Final_reads_for_mapping"
log_file=${output_dir}/merge_log.txt

mkdir -p $output_dir
: > $log_file

#========================================== Merge and Rename ====================================================
# Loop through all fastq files that contain 'part_001' in the filename
for first_file in $input_dir/*part_001_*.fastq.gz; do
    # Example: LSC2_S10_R1_001.bc1.bc2.part_001_B73-Oh43-Ki3.fastq.gz
    prefix=$(basename $first_file | sed -E 's/\.part_001_(.*)\.fastq\.gz$/\.part_/')
    sample=$(basename $first_file | sed -E 's/.*\.part_001_(.*)\.fastq\.gz/\1/')

    # Define intermediate merged file name
    merged_file="${prefix%.*}_$sample.fastq.gz"

    # Create list of the 10 part files to merge
    merged_list=""
    for i in $(seq -w 1 10); do
        merged_list+="$input_dir/${prefix}${i}_$sample.fastq.gz"$'\n'
    done

    # Merge the 10 part files into one file
    cat \
        $input_dir/${prefix}001_$sample.fastq.gz \
        $input_dir/${prefix}002_$sample.fastq.gz \
        $input_dir/${prefix}003_$sample.fastq.gz \
        $input_dir/${prefix}004_$sample.fastq.gz \
        $input_dir/${prefix}005_$sample.fastq.gz \
        $input_dir/${prefix}006_$sample.fastq.gz \
        $input_dir/${prefix}007_$sample.fastq.gz \
        $input_dir/${prefix}008_$sample.fastq.gz \
        $input_dir/${prefix}009_$sample.fastq.gz \
        $input_dir/${prefix}010_$sample.fastq.gz \
        > $output_dir/$merged_file

    # Rename the merged file to final name format: LSC2_B73-Oh43-Ki3_R1.fastq.gz
    final_name=$(echo "$merged_file" | sed -E 's/^(LSC[0-9]+)_S[0-9]+_R([0-9])_.*_(.+)\.fastq\.gz$/\1_\3_R\2.fastq.gz/')

    # Rename the merged file
    mv $output_dir/$merged_file $output_dir/$final_name

    echo -e "[Merged] $final_name\n[Input files]" >> $log_file
    echo "$merged_list" >> $log_file
    echo "[Merged into] $final_name [Done]" >> $log_file
    echo "" >> $log_file
done
