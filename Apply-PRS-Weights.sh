#!/bin/bash
#SBATCH --job-name=ApplyPRSWeights
#SBATCH --output=Logs/%x_%A_%a.out
#SBATCH --error=Logs/%x_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=104G
#SBATCH --partition=normal
#SBATCH --array=1-20

# ---------------------------------------------------------------------------
# Score participants using PLINK2 and genome-wide PRS weight files.
# Set --array=1-N where N = number of *_allchr.txt files in COMBINED_DIR.
# Over-allocating is safe — extra tasks exit gracefully.
# ---------------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
log_versions

PGEN_PREFIX="${BIM_PREFIX}"

# Discover genome-wide weight files
mapfile -t WEIGHT_FILES < <(
  find "$COMBINED_DIR" -maxdepth 1 -name "*_allchr.txt" -printf "%f\n" | sort
)

if [ "${#WEIGHT_FILES[@]}" -eq 0 ]; then
  echo "ERROR: No *_allchr.txt files found in ${COMBINED_DIR}" >&2
  echo "Run Concat-PRS-Weights.sh first." >&2
  exit 1
fi

INDEX=$((SLURM_ARRAY_TASK_ID - 1))
if [ "$INDEX" -ge "${#WEIGHT_FILES[@]}" ]; then
  echo "Task ${SLURM_ARRAY_TASK_ID} exceeds file count (${#WEIGHT_FILES[@]})." >&2
  exit 0
fi

WEIGHT_FILE="${WEIGHT_FILES[$INDEX]}"
TRAIT="${WEIGHT_FILE%_allchr.txt}"
COMBINED="${COMBINED_DIR}/${WEIGHT_FILE}"

echo "Scoring trait: $TRAIT"

# Pre-flight checks
[[ -f "${PGEN_PREFIX}.pgen" ]] || { echo "Missing: ${PGEN_PREFIX}.pgen" >&2; exit 1; }
[[ -f "$COMBINED" ]] || { echo "Missing: $COMBINED" >&2; exit 1; }

echo "  Weight file SNPs: $(wc -l < "$COMBINED")"
head -3 "$COMBINED" | sed 's/^/    /'

# Score with PLINK2
echo "Starting PLINK2 scoring at $(date)"

"$PLINK2_EXEC" \
  --pfile "$PGEN_PREFIX" \
  --score "$COMBINED" $SCORE_COLS cols=+scoresums \
  --threads "$SLURM_CPUS_PER_TASK" \
  --memory 100000 \
  --out "$SCORES_DIR/$TRAIT"

if [ -f "${SCORES_DIR}/${TRAIT}.sscore" ]; then
  N_SCORED=$(tail -n +2 "${SCORES_DIR}/${TRAIT}.sscore" | wc -l)
  echo "  Output: ${SCORES_DIR}/${TRAIT}.sscore (${N_SCORED} samples)"
fi

echo "Completed: $TRAIT at $(date)"
