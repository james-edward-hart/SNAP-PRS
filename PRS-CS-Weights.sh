#!/bin/bash
#SBATCH --job-name=PRS_CS_Weights
#SBATCH --output=Logs/%x_%A_%a.out
#SBATCH --error=Logs/%x_%A_%a.err
#SBATCH --time=5-00:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --partition=normal
#SBATCH --array=1-110

# ---------------------------------------------------------------------------
# PRS-CS weight estimation — parallelized by chromosome.
# Set --array=1-$((N_TRAITS * 22)) where N_TRAITS = rows in traits_prscs.csv
# (excluding header). Default: 5 traits × 22 = 110.
#
# Task ID layout:  TRAIT_IDX = ceil(TASK_ID / 22)
#                  CHROM     = ((TASK_ID - 1) % 22) + 1
# ---------------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
log_versions

cd "$SLURM_SUBMIT_DIR"

export PYTHONUNBUFFERED=1
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-10}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK:-10}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-10}

# Decode task ID → trait + chromosome
TRAIT_IDX=$(( (SLURM_ARRAY_TASK_ID - 1) / N_CHR + 1 ))
CHROM=$(( (SLURM_ARRAY_TASK_ID - 1) % N_CHR + 1 ))

# Read trait from CSV (format: OUT_NAME,SST_FILE,N_GWAS,POP)
LINE=$(sed -n "$((TRAIT_IDX + 1))p" "$TRAITS_PRSCS")
if [[ -z "$LINE" ]]; then
  echo "Task $SLURM_ARRAY_TASK_ID: trait index $TRAIT_IDX exceeds CSV rows. Exiting." >&2
  exit 0
fi

OUT_NAME=$(echo "$LINE" | awk -F',' '{print $1}')
SST_FILE=$(echo "$LINE" | awk -F',' '{print $2}')
N_GWAS=$(echo "$LINE" | awk -F',' '{print $3}')
POP=$(echo "$LINE" | awk -F',' '{print $4}')

SST_PATH="${GWAS_DIR}/${SST_FILE}"

# Select LD reference by population
LD_REF="${REF_DIR_PRSCS}/ldblk_ukbb_$(echo "$POP" | tr 'A-Z' 'a-z')"

# Pre-flight checks
[[ -f "$SST_PATH" ]] || { echo "Missing GWAS file: $SST_PATH" >&2; exit 1; }
[[ -d "$LD_REF" ]] || { echo "Missing LD ref: $LD_REF" >&2; exit 1; }
[[ -f "${BIM_PREFIX}.bim" ]] || { echo "Missing BIM: ${BIM_PREFIX}.bim" >&2; exit 1; }

echo "PRS-CS: ${OUT_NAME} | chr${CHROM} | pop=${POP} | task ${SLURM_ARRAY_TASK_ID}"

# Build PRS-CS command
PRSCS_CMD=("$PYTHON_BIN" "${PRSCS_DIR}/PRScs.py" \
    --ref_dir="$LD_REF" \
    --bim_prefix="$BIM_PREFIX" \
    --sst_file="$SST_PATH" \
    --n_gwas="$N_GWAS" \
    --chrom="$CHROM" \
    --out_dir="$WEIGHTS_DIR" \
    --out_name="$OUT_NAME" \
    --seed="$SEED")

if [[ "$PHI" != "auto" ]]; then
  PRSCS_CMD+=(--phi="$PHI")
fi

"${PRSCS_CMD[@]}"

echo "Completed: ${OUT_NAME} chr${CHROM}"
