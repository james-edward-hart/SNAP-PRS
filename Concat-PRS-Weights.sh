#!/bin/bash
#SBATCH --job-name=ConcatPRSWeights
#SBATCH --output=Logs/%x_%A_%a.out
#SBATCH --error=Logs/%x_%A_%a.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=normal
#SBATCH --array=1-20

# ---------------------------------------------------------------------------
# Concatenate per-chromosome PRS weight files into genome-wide files.
# Run twice — once for PRS-CS, once for PRS-CSx META weights:
#
#   sbatch Concat-PRS-Weights.sh        # PRS-CS outputs (default)
#   sbatch Concat-PRS-Weights.sh true    # PRS-CSx META outputs
#
# Set --array=1-N where N = number of trait weight sets.
# Over-allocating is safe — extra tasks exit gracefully.
# ---------------------------------------------------------------------------

set -euo pipefail

resolve_config() {
  local script_dir candidate
  local searched=()

  if [[ -n "${SNAP_PRS_CONFIG:-}" ]]; then
    searched+=("$SNAP_PRS_CONFIG")
    if [[ -f "$SNAP_PRS_CONFIG" ]]; then
      CONFIG_PATH="$(cd "$(dirname "$SNAP_PRS_CONFIG")" && pwd)/$(basename "$SNAP_PRS_CONFIG")"
      return 0
    fi
  fi

  script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  for candidate in \
    "${SLURM_SUBMIT_DIR:-}/config.sh" \
    "${PWD}/config.sh" \
    "${script_dir}/config.sh"; do
    [[ "$candidate" == "/config.sh" ]] && continue
    searched+=("$candidate")
    if [[ -f "$candidate" ]]; then
      CONFIG_PATH="$(cd "$(dirname "$candidate")" && pwd)/$(basename "$candidate")"
      return 0
    fi
  done

  echo "ERROR: Could not locate config.sh." >&2
  echo "Set SNAP_PRS_CONFIG=/path/to/config.sh or submit from the repository root." >&2
  echo "Searched:" >&2
  printf '  %s\n' "${searched[@]}" >&2
  exit 1
}

resolve_config
source "$CONFIG_PATH"

# META_ONLY: true = PRS-CSx META files only, false = all PRS-CS files
META_ONLY="${1:-false}"

if [ "$META_ONLY" = true ]; then
  PATTERN="*_META_*_chr1.txt"
  EXCLUDE=""
else
  PATTERN="*_chr1.txt"
  EXCLUDE="*_META_*_chr1.txt"
fi

# Each task discovers the trait list independently (no race condition)
mapfile -t PREFIXES < <(
  find "$WEIGHTS_DIR" -maxdepth 1 -name "$PATTERN" \
    ${EXCLUDE:+! -name "$EXCLUDE"} \
    -printf "%f\n" | sed 's/_chr1\.txt$//' | sort
)

if [ "${#PREFIXES[@]}" -eq 0 ]; then
  echo "ERROR: No $PATTERN files found in ${WEIGHTS_DIR}" >&2
  exit 1
fi

echo "Found ${#PREFIXES[@]} weight set(s): ${PREFIXES[*]}"
echo "Tip: set --array=1-${#PREFIXES[@]} to match."

# Select trait for this task
INDEX=$((SLURM_ARRAY_TASK_ID - 1))
if [ "$INDEX" -ge "${#PREFIXES[@]}" ]; then
  echo "Task ${SLURM_ARRAY_TASK_ID} exceeds trait count (${#PREFIXES[@]}). Exiting." >&2
  exit 0
fi

PREFIX="${PREFIXES[$INDEX]}"
echo "Processing: $PREFIX"

# Concatenate per-chromosome files
mkdir -p "$COMBINED_DIR"
COMBINED="${COMBINED_DIR}/${PREFIX}_allchr.txt"
rm -f "$COMBINED"

FOUND=0
for CHR in $(seq 1 "$N_CHR"); do
  CHR_FILE="${WEIGHTS_DIR}/${PREFIX}_chr${CHR}.txt"
  if [ ! -f "$CHR_FILE" ]; then
    echo "WARNING: Missing ${CHR_FILE}" >&2
    continue
  fi
  cat "$CHR_FILE" >> "$COMBINED"
  FOUND=$((FOUND + 1))
done

if [ "$FOUND" -eq 0 ]; then
  echo "ERROR: No chromosome files found for ${PREFIX}." >&2
  rm -f "$COMBINED"
  exit 1
fi

if [ ! -s "$COMBINED" ]; then
  echo "ERROR: Output file ${COMBINED} is empty." >&2
  exit 1
fi

echo "Concatenated ${FOUND}/${N_CHR} chromosome files → ${COMBINED} ($(wc -l < "$COMBINED") lines)"
