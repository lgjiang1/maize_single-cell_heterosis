#!/bin/bash
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --job-name=07b_08_downstream
#SBATCH --output=_logs/07b_08_downstream_%A.log
#SBATCH --error=_logs/07b_08_downstream_%A.log

# Phase 3 downstream wrapper: runs 07b -> 08a -> 08b sequentially in one
# SLURM job for a given (cross, cell_type) on a given FP extraction root.
#
# Defaults: B-K C5 quad on the canonical p-value extraction root.
# Override via env vars on submission.
#
# Required (or rely on default):
#   CROSS      e.g. B-K  (default B-K)
#   CELL_TYPE  e.g. C5   (default C5)
#
# Optional:
#   FP_DIR        root with per-sample hit_fp_scores.tsv.gz
#                 (default 7_fp_extraction)
#   MODEL         {v1, v2} -- passed to 08a and 08b (default v2)
#   PRESENCE_Z    presence threshold override -- currently the threshold is
#                 baked into 07a outputs (`fp_present` column). 07b/08a/08b
#                 read `fp_present` and continuous z_global; they do not
#                 re-threshold here. If you want a different threshold,
#                 re-run 07a_v2 with --presence-z.
#
# Examples:
#   # Default: B-K C5 quad on 7_fp_extraction, model v2
#   sbatch 07b_08_run_downstream.sh
#
#   # Same quad but pointed at an alternate extraction root
#   sbatch --export=ALL,FP_DIR=7_fp_extraction_alt \
#          07b_08_run_downstream.sh
#
#   # Different cell type
#   sbatch --export=ALL,CELL_TYPE=C1 \
#          07b_08_run_downstream.sh


PROJ="10_TF_Footprinting"
cd "$PROJ"

# --- Resolve env vars / defaults --------------------------------------
CROSS="${CROSS:-B-K}"
CELL_TYPE="${CELL_TYPE:-C5}"
FP_DIR="${FP_DIR:-7_fp_extraction}"
MODEL="${MODEL:-v2}"

JOINED_DIR="${FP_DIR}/joined"
CIS_TRANS_DIR="${FP_DIR}/cis_trans"
MOI_DIR="${FP_DIR}/moi"

echo "Job ID:        ${SLURM_JOB_ID}"
echo "Node:          $(hostname)"
echo "Start:         $(date)"
echo ""
echo "CROSS:         ${CROSS}"
echo "CELL_TYPE:     ${CELL_TYPE}"
echo "FP_DIR:        ${FP_DIR}"
echo "MODEL:         ${MODEL}"
echo "JOINED_DIR:    ${JOINED_DIR}"
echo "CIS_TRANS_DIR: ${CIS_TRANS_DIR}"
echo "MOI_DIR:       ${MOI_DIR}"
echo ""

# --- 07b: join quad + pair across coord systems -----------------------
echo "============================================================"
echo "[$(date)] STEP 07b: join_quad_and_pair"
echo "============================================================"
python -u 00_scripts/07_fp_extraction/07b_join_quad_and_pair.py \
    --cross "${CROSS}" \
    --cell-type "${CELL_TYPE}" \
    --fp-dir "${FP_DIR}" \
    --outdir "${JOINED_DIR}"

# --- 08a: cis/trans decomposition -------------------------------------
echo ""
echo "============================================================"
echo "[$(date)] STEP 08a: cis_trans_decomposition (model=${MODEL})"
echo "============================================================"
python -u 00_scripts/07_fp_extraction/08a_cis_trans_decomposition.py \
    --cross "${CROSS}" \
    --cell-type "${CELL_TYPE}" \
    --model "${MODEL}" \
    --joined-dir "${JOINED_DIR}" \
    --outdir "${CIS_TRANS_DIR}"

# --- 08b: mode of inheritance -----------------------------------------
echo ""
echo "============================================================"
echo "[$(date)] STEP 08b: mode_of_inheritance (model=${MODEL})"
echo "============================================================"
python -u 00_scripts/07_fp_extraction/08b_mode_of_inheritance.py \
    --cross "${CROSS}" \
    --cell-type "${CELL_TYPE}" \
    --model "${MODEL}" \
    --cis-trans-dir "${CIS_TRANS_DIR}" \
    --outdir "${MOI_DIR}"

echo ""
echo "============================================================"
echo "[$(date)] DONE: 07b -> 08a -> 08b for ${CROSS} ${CELL_TYPE}"
echo "============================================================"
echo ""
echo "Outputs:"
echo "  07b -> ${JOINED_DIR}/${CROSS}/${CELL_TYPE}/"
echo "  08a -> ${CIS_TRANS_DIR}/${CROSS}/${CELL_TYPE}/"
echo "  08b -> ${MOI_DIR}/${CROSS}/${CELL_TYPE}/"
echo ""
echo "End:           $(date)"
