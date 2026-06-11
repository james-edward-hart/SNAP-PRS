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

preflight_prscsx_alleles() {
  local chrom="$1"
  local bim_file="$2"
  shift 2

  local gwas_file
  echo "Running PRS-CSx allele preflight for chr${chrom}"
  for gwas_file in "$@"; do
    [[ -f "$gwas_file" ]] || { echo "Missing GWAS file: $gwas_file" >&2; exit 1; }

    awk -v chrom="$chrom" -v bim_file="$bim_file" -v gwas_file="$gwas_file" '
      function complement(allele) {
        if (allele == "A") return "T"
        if (allele == "T") return "A"
        if (allele == "C") return "G"
        if (allele == "G") return "C"
        return "N"
      }

      function compatible(bim_a1, bim_a2, gwas_a1, gwas_a2) {
        if (bim_a1 == gwas_a1 && bim_a2 == gwas_a2) return 1
        if (bim_a1 == gwas_a2 && bim_a2 == gwas_a1) return 1
        if (complement(bim_a1) == gwas_a1 && complement(bim_a2) == gwas_a2) return 1
        if (complement(bim_a1) == gwas_a2 && complement(bim_a2) == gwas_a1) return 1
        return 0
      }

      NR == FNR {
        if ($1 == chrom || $1 == "chr" chrom) {
          bim_a1 = toupper($5)
          bim_a2 = toupper($6)
          if (bim_a1 ~ /^[ACGT]$/ && bim_a2 ~ /^[ACGT]$/) {
            bim_allele1[$2] = bim_a1
            bim_allele2[$2] = bim_a2
            bim_count++
          }
        }
        next
      }

      FNR == 1 {
        header_seen = 1
        for (i = 1; i <= NF; i++) {
          header = toupper($i)
          if (header == "SNP") snp_col = i
          if (header == "A1") a1_col = i
          if (header == "A2") a2_col = i
        }
        if (!snp_col || !a1_col || !a2_col) {
          printf "ERROR: %s must have SNP, A1, and A2 columns in the header.\n", gwas_file > "/dev/stderr"
          fatal = 1
          exit 2
        }
        next
      }

      {
        snp = $snp_col
        if (snp in bim_allele1) {
          overlap++
          gwas_a1 = toupper($a1_col)
          gwas_a2 = toupper($a2_col)
          bad_alleles = (gwas_a1 !~ /^[ACGT]$/ || gwas_a2 !~ /^[ACGT]$/)
          if (bad_alleles || !compatible(bim_allele1[snp], bim_allele2[snp], gwas_a1, gwas_a2)) {
            mismatches++
            if (!example_seen) {
              example_seen = 1
              example_snp = snp
              example_bim = bim_allele1[snp] "/" bim_allele2[snp]
              example_gwas = gwas_a1 "/" gwas_a2
            }
          }
        }
      }

      END {
        if (fatal) exit 2
        if (!header_seen) {
          printf "ERROR: %s is empty or missing a header.\n", gwas_file > "/dev/stderr"
          exit 2
        }

        printf "Allele preflight: file=%s chr=%s bim_snps=%d overlap=%d mismatches=%d\n", \
          gwas_file, chrom, bim_count, overlap, mismatches

        if (bim_count == 0) {
          printf "ERROR: No chr%s A/C/G/T SNPs found in BIM: %s\n", chrom, bim_file > "/dev/stderr"
          exit 3
        }
        if (overlap == 0) {
          printf "ERROR: Zero overlapping SNPs for file=%s against chr%s BIM SNPs.\n", gwas_file, chrom > "/dev/stderr"
          exit 4
        }
        if (mismatches > 0) {
          printf "ERROR: %d allele mismatch(es) for file=%s against chr%s BIM SNPs.\n", mismatches, gwas_file, chrom > "/dev/stderr"
          printf "Example mismatch: SNP=%s file=%s BIM=%s GWAS=%s\n", example_snp, gwas_file, example_bim, example_gwas > "/dev/stderr"
          exit 5
        }
      }
    ' "$bim_file" "$gwas_file"
  done
}

resolve_config
source "$CONFIG_PATH"
log_versions prscsx
check_prscs_python_deps

cd "${SLURM_SUBMIT_DIR:-$PIPELINE_DIR}"

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
IFS=',' read -r -a SST_FILE_ARRAY <<< "$SST_FILES"
SST_PATH_ARRAY=()
for SST_FILE in "${SST_FILE_ARRAY[@]}"; do
  SST_FILE="${SST_FILE#"${SST_FILE%%[![:space:]]*}"}"
  SST_FILE="${SST_FILE%"${SST_FILE##*[![:space:]]}"}"
  SST_PATH_ARRAY+=("${GWAS_DIR}/${SST_FILE}")
done
SST_PATHS=$(IFS=','; echo "${SST_PATH_ARRAY[*]}")

# Pre-flight checks
[[ -d "$REF_DIR_PRSCSX" ]] || { echo "Missing LD ref dir: $REF_DIR_PRSCSX" >&2; exit 1; }
[[ -f "${BIM_PREFIX}.bim" ]] || { echo "Missing BIM: ${BIM_PREFIX}.bim" >&2; exit 1; }
preflight_prscsx_alleles "$CHROM" "${BIM_PREFIX}.bim" "${SST_PATH_ARRAY[@]}"

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
