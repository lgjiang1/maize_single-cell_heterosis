#!/bin/bash
#SBATCH --job-name=filter
#SBATCH --output=logs/filter_by_replicates_%A_%a.out
#SBATCH --error=logs/filter_by_replicates_%A_%a.err
#SBATCH --time=00:10:00
#SBATCH --mem=200mb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --array=1-210     #210 pairs (rep1 & rep2)

input_dir="4_filter_by_replicates/all_filtered_peaks"
output_dir="4_filter_by_replicates/final_fltered_peaks"

mkdir -p ${output_dir}

# ------- Build pair keys (normalize _rep1/_rep2 -> _repX) ------------
# example：OxK_rep1-C3_Oh43.filtered.bed ↔ OxK_rep2-C3_Oh43.filtered.bed
mapfile -t keys < <(
  ls "${input_dir}"/*_rep[12]-*.filtered.bed \
    | sed -E 's/_rep[12]-/_repX-/' \
    | xargs -n1 basename \
    | sort -u
)

key="${keys[$((SLURM_ARRAY_TASK_ID - 1))]}"
rep1_name="${key/_repX-/_rep1-}"
rep2_name="${key/_repX-/_rep2-}"

rep1_peak="${input_dir}/${rep1_name}"
rep2_peak="${input_dir}/${rep2_name}"

rep1_base="${rep1_name%.filtered.bed}"
rep2_base="${rep2_name%.filtered.bed}"

# ---- Output file names (500bp-expanded & reproducible) ----
out_rep1="${output_dir}/${rep1_base}.reproducible.bed"
out_rep2="${output_dir}/${rep2_base}.reproducible.bed"


# ---- Pairwise reproducibility filter (both directions, once per pair) ----
# Keep peaks that overlap between replicates (>1bp overlap)
bedtools intersect -u -a ${rep1_peak} -b ${rep2_peak}  > ${out_rep1}
bedtools intersect -u -a ${rep2_peak} -b ${rep1_peak}  > ${out_rep2}

# ------ Replace summit offset with absolute summit position ------
awk 'BEGIN{OFS="\t"} { $6 = $2 + 1 + $6; print }' ${out_rep1} > ${out_rep1}.tmp && mv ${out_rep1}.tmp ${out_rep1}
awk 'BEGIN{OFS="\t"} { $6 = $2 + 1 + $6; print }' ${out_rep2} > ${out_rep2}.tmp && mv ${out_rep2}.tmp ${out_rep2}

# --- Replace cols 2-3 with summit±250bp (exact 501bp window) ---
awk 'BEGIN{OFS="\t"} {start=$6-251; if(start<0)start=0; end=$6+250; $2=start; $3=end; print}' ${out_rep1} > ${out_rep1}.tmp && mv ${out_rep1}.tmp ${out_rep1}
awk 'BEGIN{OFS="\t"} {start=$6-251; if(start<0)start=0; end=$6+250; $2=start; $3=end; print}' ${out_rep2} > ${out_rep2}.tmp && mv ${out_rep2}.tmp ${out_rep2}

echo "Pair $SLURM_ARRAY_TASK_ID done:"
echo "  $out_rep1"
echo "  $out_rep2"
