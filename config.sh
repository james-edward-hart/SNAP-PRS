#!/bin/bash
# ============================================================================
# PRS Pipeline — Centralized Configuration
# ============================================================================
# Edit the variables below to match your environment.
# All shell scripts in this pipeline source this file.
# ============================================================================

# -- Project root directory
PROJECT_DIR="/path/to/project"

# -- Directory containing formatted GWAS summary statistics
#    Expected format: tab-separated, columns: SNP  A1  A2  BETA  SE (or P)
GWAS_DIR="${PROJECT_DIR}/GWAS-PRS-Ready"

# -- Directory where VCF genotype files are stored
#    Expected naming: chr1_dose.vcf.gz ... chr22_dose.vcf.gz (see VCF_PATTERN)
VCF_DIR="${PROJECT_DIR}/Genotype-Data"

# -- Output directory for cleaned/QC'd genotype data
CLEANED_DIR="${PROJECT_DIR}/Cleaned-Genetic-Data"

# -- Output directory for per-chromosome PRS weights (PRS-CS/CSx output)
WEIGHTS_DIR="${PROJECT_DIR}/PRS-Weights"

# -- Output directory for genome-wide concatenated weight files
COMBINED_DIR="${PROJECT_DIR}/Genome-Wide-Weights"

# -- Output directory for PLINK2 PRS scores
SCORES_DIR="${PROJECT_DIR}/PRS-Scores"

# -- Log directory (created automatically)
LOGS_DIR="${PROJECT_DIR}/Logs"

# -- Prefix for the final genotype files (without .pgen/.bim extension)
#    Must point to the files Apply-PRS-Weights.sh will score (pgen) and
#    PRS-CS/CSx will reference (bim). Genotype-QC.sh produces these
#    automatically. Skip-QC users: set this to your existing file prefix
#    and ensure a .bim file exists (plink2 --pfile X --make-just-bim --out X).
BIM_PREFIX="${CLEANED_DIR}/genome_wide.QC"

# -- Trait configuration CSVs
TRAITS_PRSCS="${PROJECT_DIR}/PRS-Pipeline/traits_prscs.csv"
TRAITS_PRSCSX="${PROJECT_DIR}/PRS-Pipeline/traits_prscsx.csv"

# ============================================================================
# Software paths
# ============================================================================

# -- PLINK2 executable
PLINK2_EXEC="/path/to/plink2"

# -- bcftools executable
BCFTOOLS_EXEC="/path/to/bcftools"

# -- Python binary (with PRS-CS/CSx dependencies installed)
PYTHON_BIN="/path/to/python"

# -- PRS-CS source directory (contains PRScs.py)
PRSCS_DIR="/path/to/PRScs"

# -- PRS-CSx source directory (contains PRScsx.py)
PRSCSX_DIR="/path/to/PRScsx"

# -- LD reference panel directory for PRS-CS
#    Should contain subdirectories like ldblk_ukbb_eur, ldblk_ukbb_afr, etc.
REF_DIR_PRSCS="/path/to/ld-reference/PRS-CS"

# -- LD reference panel directory for PRS-CSx
REF_DIR_PRSCSX="/path/to/ld-reference/PRS-CSx"

# ============================================================================
# Genotype QC mode
# ============================================================================
# QC_MODE controls how Genotype-QC.sh runs:
#   "full"    — Per-chromosome VCFs → concat → QC → pgen  (uses VCF_DIR + VCF_PATTERN)
#   "qc_only" — Already-merged VCF  → QC → pgen           (uses MERGED_VCF)
#   If your genotypes are already QC'd pgen files, skip Genotype-QC.sh entirely
#   and just set BIM_PREFIX above to point to your existing files. Note you still
#   need to provide a BIM file, create this with plink: --make-just-bim if needed
QC_MODE="full"

# -- Path to an already-merged genome-wide VCF (only used when QC_MODE="qc_only")
MERGED_VCF=""

# ============================================================================
# QC parameters
# ============================================================================

# -- Minor allele frequency threshold
MAF_THRESHOLD=0.01

# -- Imputation quality (INFO R2) threshold
R2_THRESHOLD=0.3

# ============================================================================
# Pipeline parameters
# ============================================================================

# -- Random seed for PRS-CS/CSx reproducibility
SEED=42

# -- Global shrinkage parameter (phi) for PRS-CS/CSx
#    "auto" = fully Bayesian (recommended default)
#    See https://github.com/getian107/PRScsx for guidance on setting this
PHI="auto"

# -- PRS-CSx meta flag: combine populations via meta-analysis
#    "True" = produce META output (recommended default)
#    See https://github.com/getian107/PRScsx for details
META_FLAG="True"

# -- Number of autosomes
N_CHR=22

# -- VCF file naming pattern (use ${chr} as placeholder for chromosome number)
#    Example: chr1_dose.vcf.gz
VCF_PATTERN='chr${chr}_dose.vcf.gz'

# -- PRS-CS(x) output score columns for PLINK2 --score
#    Default: col2=SNP, col4=A1, col6=BETA
SCORE_COLS="2 4 6"

# ============================================================================
# Setup (runs on source)
# ============================================================================
mkdir -p "$CLEANED_DIR" "$WEIGHTS_DIR" "$COMBINED_DIR" "$SCORES_DIR" "$LOGS_DIR"

# -- Log software versions to SLURM output
log_versions() {
  echo "=== Software Versions ==="
  "$PLINK2_EXEC" --version 2>/dev/null | head -1 || echo "plink2: not found"
  "$BCFTOOLS_EXEC" --version 2>/dev/null | head -1 || echo "bcftools: not found"
  "$PYTHON_BIN" --version 2>&1 || echo "python: not found"
  R --version 2>/dev/null | head -1 || echo "R: not found"
  echo "========================="
}
