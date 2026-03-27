#!/bin/bash
#SBATCH --job-name=GenoQC
#SBATCH --output=Logs/%x_%j.out
#SBATCH --error=Logs/%x_%j.err
#SBATCH --time=5-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --partition=normal

# Genotype QC: merge VCFs (if needed), filter MAF + R2, remove ambiguous SNPs, output pgen.
#
# Controlled by QC_MODE in config.sh:
#   "full"    — Concatenate per-chromosome VCFs, then QC
#   "qc_only" — Use an existing merged VCF, just QC
#
# If your genotypes are already QC'd pgen files, skip this script entirely
# and set BIM_PREFIX in config.sh to point to your existing files.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
log_versions

# Uncomment and adjust if bcftools needs non-standard shared libraries
# export LD_LIBRARY_PATH="${HOME}/local/lib:${LD_LIBRARY_PATH:-}"

# ============================================================
# Step 1: Get genome-wide VCF (mode-dependent)
# ============================================================
if [[ "$QC_MODE" == "full" ]]; then
  echo "QC_MODE=full — concatenating per-chromosome VCFs"

  VCF_LIST="$CLEANED_DIR/vcf_list.txt"
  : > "$VCF_LIST"
  for chr in $(seq 1 "$N_CHR"); do
    vcf=$(eval echo "$VCF_PATTERN")
    [[ -f "$VCF_DIR/$vcf" ]] || { echo "Missing VCF: $VCF_DIR/$vcf" >&2; exit 1; }
    echo "$VCF_DIR/$vcf" >> "$VCF_LIST"
  done

  GENOME_VCF="$CLEANED_DIR/genome_wide.vcf.gz"
  "$BCFTOOLS_EXEC" concat \
    -f "$VCF_LIST" \
    -Oz \
    -o "$GENOME_VCF"
  "$BCFTOOLS_EXEC" index "$GENOME_VCF"

elif [[ "$QC_MODE" == "qc_only" ]]; then
  echo "QC_MODE=qc_only — using existing merged VCF"
  [[ -n "$MERGED_VCF" ]] || { echo "MERGED_VCF is empty in config.sh" >&2; exit 1; }
  [[ -f "$MERGED_VCF" ]] || { echo "MERGED_VCF not found: $MERGED_VCF" >&2; exit 1; }
  GENOME_VCF="$MERGED_VCF"

else
  echo "Unknown QC_MODE: $QC_MODE (expected 'full' or 'qc_only')" >&2
  exit 1
fi

# ============================================================
# Step 2: Temporary BED for ambiguous SNP detection
# ============================================================
"$PLINK2_EXEC" \
  --vcf "$GENOME_VCF" \
  --maf "$MAF_THRESHOLD" \
  --extract-if-info "R2 > $R2_THRESHOLD" \
  --threads "${SLURM_CPUS_PER_TASK:-8}" \
  --make-bed \
  --out "$CLEANED_DIR/genome_wide"

awk '($5=="A" && $6=="T") || ($5=="T" && $6=="A") || \
     ($5=="C" && $6=="G") || ($5=="G" && $6=="C") {print $2}' \
  "$CLEANED_DIR/genome_wide.bim" > "$CLEANED_DIR/ambiguous_snps.txt"

# ============================================================
# Step 3: Final QC pgen
# ============================================================
"$PLINK2_EXEC" \
  --vcf "$GENOME_VCF" dosage=DS \
  --exclude "$CLEANED_DIR/ambiguous_snps.txt" \
  --maf "$MAF_THRESHOLD" \
  --extract-if-info "R2 > $R2_THRESHOLD" \
  --threads "${SLURM_CPUS_PER_TASK:-8}" \
  --make-pgen \
  --out "$CLEANED_DIR/genome_wide.QC"

# Generate .bim for PRS-CS/CSx (needs BIM-format variant list)
"$PLINK2_EXEC" \
  --pfile "$CLEANED_DIR/genome_wide.QC" \
  --threads "${SLURM_CPUS_PER_TASK:-8}" \
  --make-just-bim \
  --out "$CLEANED_DIR/genome_wide.QC"

echo "Genotype QC complete: $CLEANED_DIR/genome_wide.QC"
