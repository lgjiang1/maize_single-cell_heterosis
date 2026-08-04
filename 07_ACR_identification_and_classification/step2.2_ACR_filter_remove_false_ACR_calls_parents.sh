#!/bin/bash
#SBATCH --job-name=filter_ACRs
#SBATCH --output=logs/filter_false_ACRs_%A_%a.out
#SBATCH --error=logs/filter_false_ACRs_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=1gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --array=1-28

ml Bioinformatics bedops

BASE="3_remove_false_ACR_calls/filtered_peaks_in_parents/B73"   # Replace B73 (Ki3 or Oh43) with target genotype if needed
PEAK_DIR=${BASE}
TN5_DIR=${BASE}
GENOME=${BASE}/B73_genome_size           # Replace B73 (Ki3 or Oh43) with target genotype if needed
MAP=${BASE}/Zm-B73.mappable_genomic_regions.bed   # Replace B73 (Ki3 or Oh43) with target genotype if needed
WORK=${BASE}/work

PREFIX=$(printf '%s\n' "${PEAK_DIR}"/*_peaks.narrowPeak | sed -n "${SLURM_ARRAY_TASK_ID}{s#.*/##; s/_peaks\.narrowPeak$//; p}")
echo "[INFO] Processing ${PREFIX}"
mkdir -p "${WORK}/${PREFIX}"
RAW_PEAK=${PEAK_DIR}/${PREFIX}_peaks.narrowPeak 
TN5=${TN5_DIR}/${PREFIX}.tn5.bed

# beta filter parameters (can be adjusted as needed)
min_prop=0.01
min_cell=1
fdr_cut=0.01

# =========================== Step 1: extract and sort the peak file ============================
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$9,$10}' "${RAW_PEAK}" > "${WORK}/${PREFIX}/${PREFIX}.peaks6.bed"
sort -k1,1 -k2,2n "${WORK}/${PREFIX}/${PREFIX}.peaks6.bed" > "${WORK}/${PREFIX}/${PREFIX}.peaks.sorted.bed"

# =================== Step 2: generate control (background) regions =================
bedtools shuffle \
  -i ${WORK}/${PREFIX}/${PREFIX}.peaks.sorted.bed \
  -g ${GENOME} \
  -incl ${MAP} \
  -excl ${WORK}/${PREFIX}/${PREFIX}.peaks.sorted.bed \
  -seed 2025 \
| sort -k1,1 -k2,2n > ${WORK}/${PREFIX}/${PREFIX}.control.bed

# ========= Step 3: count number of unique cells =============
N_CELLS=$(cut -f4 ${TN5} | sort -u | wc -l)
echo "[INFO] ${PREFIX}: ${N_CELLS} cells"

# ========= Step 4: compute accessibility fraction (p) for background regions ========
bedmap --echo --echo-map-id-uniq --delim '\t' \
  ${WORK}/${PREFIX}/${PREFIX}.control.bed ${TN5} \
| awk -v N="${N_CELLS}" -F'\t' 'BEGIN{OFS="\t"}{
     n = ($NF=="."||$NF=="") ? 0 : split($NF,a,";");
     p = (N>0)? n/N : 0;
     print $1,$2,$3,p;
   }' > ${WORK}/${PREFIX}/${PREFIX}.bg.p.bed


# ============ Step 5: compute accessibility fraction (p) for candidate peaks ==========
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"peak_"NR,$0}' \
  ${WORK}/${PREFIX}/${PREFIX}.peaks.sorted.bed \
> ${WORK}/${PREFIX}/${PREFIX}.peaks.withID.bed

bedmap --echo --echo-map-id-uniq --delim '\t' \
  ${WORK}/${PREFIX}/${PREFIX}.peaks.withID.bed ${TN5} \
| awk -v N="${N_CELLS}" 'BEGIN{OFS="\t"}{
     m=NF; split($(m),a,";"); n=($(m)=="."||$(m)=="")?0:length(a);
     p=(N>0)?n/N:0; print $0,p;
   }' > ${WORK}/${PREFIX}/${PREFIX}.peaks.p.tsv

# =========== Step 6: run R script to fit beta distribution and filter peaks ===========
Rscript ./ecdf_filter.R ${PREFIX} ${BASE} ${fdr_cut} ${min_prop} ${min_cell}


