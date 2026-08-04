#!/bin/bash
#SBATCH --job-name=processing
#SBATCH --output=logs/processing.out
#SBATCH --error=logs/processing.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard

base_dir="4_genotyping_QC/5_Cell_validation"
csv_list="${base_dir}/csv_files_to_merge.txt"
output="${base_dir}/merged_allele_ratio_matrix.csv"
final_output="${base_dir}/final_merged_allele_ratio_matrix.tsv"

find "$base_dir" -type f -name "*_allele_ratio.csv" | sort > ${csv_list}

# =========================== Step 1: Build merged header ===============================
echo "Building header..."
header=""
mapfile -t csvs < ${csv_list}

for i in "${!csvs[@]}"; do
  f="${csvs[i]}"
  h=$(head -n1 "$f")
  if [ "$i" -eq 0 ]; then
    header="$h"
  else
    h_trimmed=$(echo "$h" | cut -d',' -f2-)
    header="$header,$h_trimmed"
  fi
done

echo $header > $output

# ========================== Step 2: Merge data rows ====================================
echo "Merging data..."
paste_cmd="paste -d','"

for i in "${!csvs[@]}"; do
  f="${csvs[i]}"
  if [ "$i" -eq 0 ]; then
    paste_cmd="$paste_cmd <(tail -n +2 \"$f\")"
  else
    paste_cmd="$paste_cmd <(tail -n +2 \"$f\" | cut -d',' -f2-)"
  fi
done

# =========================== Step 3: Execute paste and write to output =================
eval ${paste_cmd} >> $output

sed -i '1s/^ *,/SNP,/' $output   #Replace an empty or space-only header in the first column with "SNP"

# ========================== Step 4: Recode CSV data with missing values and genotype labels ================
echo "Starting CSV processing: recoding values and converting to tab-delimited format..."
awk -F',' -v OFS='\t' -v out="$final_output" '
NR==1 {
    printf "%s", $1 > out;
    for (i=2; i<=NF; i++) printf "\t%s", $i >> out;
    printf "\n" >> out;
    next;
}
{
    # Progress display every 1000 data rows
    if ((NR - 1) % 1000 == 0) {
        printf("[INFO] Processed %d rows...\n", NR - 1) > "/dev/stderr";
    }

    for (i=1; i<=NF; i++) {
        val = $i;
        gsub(/^ +| +$/, "", val);

        if (val == "") {
            val = ".";
        } else if (val ~ /^-?[0-9.]+([eE][-+]?[0-9]+)?$/) {
            v = val + 0;

            ## Recode based on value range
            if (v >= 0.0 && v <= 0.2) val = "0/0";
            else if (v > 0.2 && v < 0.8) val = "0/1";
            else if (v >= 0.8 && v <= 1.0) val = "1/1";
        }
        $i = val;
    }

    printf "%s", $1 >> out;
    for (i=2; i<=NF; i++) printf "\t%s", $i >> out;
    printf "\n" >> out;
}
' "$output"


input=$final_output
prefix="output_part_"
chunk_size=2000
outdir="${base_dir}/chunks" 
mkdir $outdir

total_cols=$(head -n1 "$input" | awk -F'\t' '{print NF}')
data_cols=$((total_cols - 1))

start=2
part=1

# Loop through the data columns in chunks
while [ $start -le $total_cols ]; do
    end=$((start + chunk_size - 1))
    if [ $end -gt $total_cols ]; then
        end=$total_cols
    fi

    echo "Processing columns $start to $end ..."

    cut_fields="1,${start}-${end}"
    out="${outdir}/${prefix}${part}.txt"

    # Extract the selected columns and write to output file
    cut -f$cut_fields "$input" > "$out"

    echo " --> Written to $out"

    # Move to next chunk
    start=$((end + 1))
    part=$((part + 1))
done
