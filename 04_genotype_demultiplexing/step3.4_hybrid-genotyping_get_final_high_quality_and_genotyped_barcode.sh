#!/bin/bash
cd 3_hybrid_genotyping_based_on_read_counts

count_dir="3_genotype-specific-read_count_file"
bayes_dir="2_bayes_results"
out_dir="4_final_high_quality_and_genotyped_barcode"

mkdir -p "$out_dir"

for count_file in ${count_dir}/*_count.txt
do
    base=$(basename "$count_file" _count.txt)

    bayes_file="${bayes_dir}/${base}_perfect_F1_cells.txt"
    filtered_count="${out_dir}/${base}_filtered_count.txt"
    low_fraction_barcode="${out_dir}/${base}_barcode.txt"

    # Step 1:
    # Keep only the rows in the count file whose barcodes are present
    # in the corresponding perfect_F1_cells file.
    awk '
        BEGIN { FS=OFS="\t" }
        NR==FNR {
            if (FNR > 1) keep[$1] = 1
            next
        }

        FNR == 1 {
            print
            next
        }

        ($1 in keep) {
            print
        }
    ' "$bayes_file" "$count_file" > "$filtered_count"

    # Step 2:
    # From the filtered count file, keep only rows with fraction < 0.05.
    awk '
        BEGIN { FS=OFS="\t" }

        FNR > 1 && $4 <= 0.05 {
            print $1
        }
    ' "$filtered_count" > "$low_fraction_barcode"

    echo "Done: $base"
done

rm ${out_dir}/*_filtered_count.txt

