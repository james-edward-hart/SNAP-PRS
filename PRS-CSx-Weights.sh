#!/bin/bash
#SBATCH --job-name=PRS_CSx_Weights
#SBATCH --output=Logs/%x_%A_%a.out
#SBATCH --error=Logs/%x_%A_%a.err
#SBATCH --time=5-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --partition=normal
#SBATCH --array=1-220

# ---------------------------------------------------------------------------
# PRS-CSx weight estimation — parallelized by chromosome.
# Set --array=1-$((N_TRAITS * 22)) where N_TRAITS = rows in traits_prscsx.csv
# (excluding header). Default: 10 traits × 22 = 220.
#
# IMPORTANT: With a 220-job HPC cap, submit this AFTER PRS-CS jobs finish.
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
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

# Decode task ID → trait + chromosome
TRAIT_IDX=$(( (SLURM_ARRAY_TASK_ID - 1) / N_CHR + 1 ))
CHROM=$(( (SLURM_ARRAY_TASK_ID - 1) % N_CHR + 1 ))

# Read trait from CSV (format: OUT_NAME,"SST_FILES","N_GWAS","POP")
LINE=$(sed -n "$((TRAIT_IDX + 1))p" "$TRAITS_PRSCSX")
if [[ -z "$LINE" ]]; then
  echo "Task $SLURM_ARRAY_TASK_ID: trait index $TRAIT_IDX exceeds CSV rows. Exiting." >&2
  exit 0
fi

OUT_NAME=$(echo "$LINE" | awk -F',' '{print $1}')
SST_FILES=$(echo "$LINE" | awk -F'"' '{print $2}')
N_GWAS=$(echo "$LINE" | awk -F'"' '{print $4}')
POP=$(echo "$LINE" | awk -F'"' '{print $6}')

if [[ -z "$SST_FILES" || -z "$N_GWAS" || -z "$POP" ]]; then
  echo "ERROR: Failed to parse trait $TRAIT_IDX from $TRAITS_PRSCSX" >&2
  echo "  Ensure multi-value fields are double-quoted in the CSV." >&2
  echo "  LINE: $LINE" >&2
  exit 1
fi

# Add full path to each summary stats file
SST_PATHS=$(echo "$SST_FILES" | sed "s|[^,]*|${GWAS_DIR}/&|g")

# Pre-flight checks
[[ -d "$REF_DIR_PRSCSX" ]] || { echo "Missing LD ref dir: $REF_DIR_PRSCSX" >&2; exit 1; }
[[ -f "${BIM_PREFIX}.bim" ]] || { echo "Missing BIM: ${BIM_PREFIX}.bim" >&2; exit 1; }

echo "PRS-CSx: ${OUT_NAME} | chr${CHROM} | pop=${POP} | task ${SLURM_ARRAY_TASK_ID}"
echo "Files: $SST_PATHS"
echo "N: $N_GWAS"

# Build PRS-CSx command
PRSCSX_CMD=("$PYTHON_BIN" "${PRSCSX_DIR}/PRScsx.py" \
    --ref_dir="$REF_DIR_PRSCSX" \
    --bim_prefix="$BIM_PREFIX" \
    --sst_file="$SST_PATHS" \
    --n_gwas="$N_GWAS" \
    --pop="$POP" \
    --chrom="$CHROM" \
    --out_dir="$WEIGHTS_DIR" \
    --out_name="$OUT_NAME" \
    --meta="$META_FLAG" \
    --seed="$SEED")

if [[ "$PHI" != "auto" ]]; then
  PRSCSX_CMD+=(--phi="$PHI")
fi

"${PRSCSX_CMD[@]}"

echo "Completed: ${OUT_NAME} chr${CHROM}"
