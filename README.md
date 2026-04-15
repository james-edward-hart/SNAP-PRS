# SNAP-PRS 

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

**Author:** James E. Hart (Jimmy) — [James.Hart@VCUhealth.org](mailto:James.Hart@VCUhealth.org)

SNAP-PRS (SLURM Native Automated Pipeline - Polygenic Risk Scores) is a robust, light-weight, SLURM based, config-driven, pipeline for computing Polygenic Scores (PRS) using PRS-CS and PRS-CSx, scoring with PLINK2, and building a QC'd final dataset. The pipeline allows for simple setup and manual job submission, avoiding bulky workflows. If you are comfortable writing scripts directly and submitting them in SLURM, you'll feel at home and find this tool helpful.

Please read the documentation below carefully before proceeding with analysis. 

## Citations

If you use this pipeline, please cite the following:

**This pipeline:**

Cite the Zenodo DOI for this repository: [![DOI](https://zenodo.org/badge/1193776964.svg)](https://doi.org/10.5281/zenodo.19258447) Alternatively, contact me at James.Hart@vcuhealth.org for consultation. 


**Please also cite the tools and data this pipeline depends on:**
- [PRS-CS](https://github.com/getian107/PRScs)
- [PRS-CSx](https://github.com/getian107/PRScsx)
- [UK Biobank LD reference panels](https://github.com/getian107/PRScsx#ld-reference-panels)
- [PLINK2](https://www.cog-genomics.org/plink/2.0/)
- [bcftools](http://www.htslib.org/download/)
- [Python 3](https://www.python.org/) with [scipy](https://scipy.org/), [h5py](https://www.h5py.org/), [numpy](https://numpy.org/)
- [R](https://www.r-project.org/) with [data.table](https://rdatatable.gitlab.io/data.table/), [ggplot2](https://ggplot2.tidyverse.org/), [scales](https://scales.r-lib.org/)

## Overview

1. **Genotype QC** — (Optional) Merge and/or filter imputed VCFs by MAF, imputation quality, and strand ambiguity
2. **PRS-CS** — Estimate posterior SNP effect weights for single-ancestry traits
3. **PRS-CSx** — Estimate posterior SNP effect weights for multi-ancestry traits
4. **Concatenate** — Combine per-chromosome weights into genome-wide files
5. **Score** — Apply weights to genotype data using PLINK2
6. **Final QC** — Inspect distributions, filter to unrelated individuals, standardize, output final dataset

## Prerequisites

### Software
- [PLINK2](https://www.cog-genomics.org/plink/2.0/)
- [bcftools](http://www.htslib.org/download/)
- [PRS-CS](https://github.com/getian107/PRScs) and/or [PRS-CSx](https://github.com/getian107/PRScsx)
- Python 3 with scipy, h5py, numpy
- R with `data.table`, `ggplot2`, `scales`

### Data
- **LD reference panels** — UKBB LD matrices for PRS-CS/CSx (download from their GitHub repos)
- **Genotype data** — Imputed VCFs with INFO/R2 field and dosage (DS) format
- **GWAS summary statistics** — Formatted to the 5-column spec below

### Genome Build

This pipeline assumes **GRCh37 (hg19)**. If your data is on GRCh38, you must either liftover your GWAS summary stats or use GRCh38-compatible LD reference panels.

## Quick Start

Clone this full repository to your HPC (only 224kb)\
`git clone https://github.com/james-edward-hart/SNAP-PRS`

1. **Edit `config.sh`** — Set all paths and software locations
2. **Prepare GWAS summary stats** in the required format (see below)
3. **Populate trait CSVs** — `traits_prscs.csv` and `traits_prscsx.csv` (shipped with examples — replace with your own)
4. **Submit jobs** in the order documented below. **IMPORTANT** -Adjust --array at the top of each script prior to submission to match your trait count. [See Calculating '--arrary' Sizes below]
5. *(Optional)* **Edit paths at the top of `prs-qc.r`** and run it for score inspection, standardization, and plots

**`config.sh` is the only file you need to edit for the core pipeline (plus --array sizes in scripts).** The `prs-qc.r` script is an example post-processing step with its own path variables at the top.

## GWAS Input Format

For more detailed info see: https://github.com/getian107/PRScsx
Tab-separated, one row per SNP, with header:

```
SNP     A1      A2      BETA    SE
```

| Column | Description |
|--------|-------------|
| SNP | rsID (e.g. rs12345) |
| A1 | Effect allele (uppercase) |
| A2 | Other allele (uppercase) |
| BETA | Effect size (log-OR, beta, or Z-score) |
| SE | Standard error (alternatively, replace SE with P for p-value) |

Only `rs`-prefixed SNPs are supported. The included `gwas-formatting-prs.r` in `Scripts/` is an ABCD-specific example — adapt for your own GWAS sources.

## Trait CSV Format

### traits_prscs.csv (PRS-CS — single-ancestry traits)

```
OUT_NAME,SST_FILE,N_GWAS,POP
ADHD2022,ADHD2022_EUR.txt,128214,EUR
```

| Column | Description |
|--------|-------------|
| OUT_NAME | Output label for this trait |
| SST_FILE | Filename in GWAS_DIR |
| N_GWAS | GWAS sample size |
| POP | Population (EUR, AFR, EAS, etc.) — used to select LD reference |

### traits_prscsx.csv (PRS-CSx — multi-ancestry traits)

```
OUT_NAME,SST_FILES,N_GWAS,POP
AlcDep2018,"AlcDep2018_EUR.txt,AlcDep2018_AFR.txt","8922,5732","EUR,AFR"
```

Comma-separated values within quotes for multi-ancestry fields.

## Job Submission Order

Submit steps sequentially. If your cluster has a job submission cap, stagger accordingly.

```bash
# Step 1: Genotype QC (1 job — skip if your pgen files are already QC'd)
sbatch PRS-Pipeline/Genotype-QC.sh

# Step 2: PRS-CS weights (N Analyses × 22 jobs, e.g. 5 × 22 = 110)
# Wait for genotype QC to finish first.
sbatch PRS-Pipeline/PRS-CS-Weights.sh

# Step 3: PRS-CSx weights (N Analyses × 22 jobs, e.g. 10 × 22 = 220)
# *** Submit ONLY after all PRS-CS jobs complete to stay under 220-job cap ***
sbatch PRS-Pipeline/PRS-CSx-Weights.sh

# Step 4: Concatenate per-chromosome weights
sbatch PRS-Pipeline/Concat-PRS-Weights.sh         # PRS-CS outputs (default)
sbatch PRS-Pipeline/Concat-PRS-Weights.sh true     # PRS-CSx META outputs

# Step 5: Score participants
sbatch PRS-Pipeline/Apply-PRS-Weights.sh

# Step 6: QC and final dataset
Rscript PRS-Pipeline/prs-qc.r
```

### Calculating `--array` sizes

For PRS-CS/CSx scripts, set `--array=1-$((N_TRAITS * 22))`:
- Count data rows in your CSV: `tail -n +2 traits_prscs.csv | wc -l` → N_TRAITS
- Multiply by 22: e.g., 5 traits → `--array=1-110`

For Concat and Apply scripts, set `--array=1-N` where N is the number of trait weight sets. Over-allocating is safe — extra tasks exit gracefully.

## Genotype QC Entry Points

Set `QC_MODE` in `config.sh` based on your starting data:

| Starting data | `QC_MODE` | What to set | What happens |
|---|---|---|---|
| Per-chromosome VCFs | `full` | `VCF_DIR`, `VCF_PATTERN` | Concat all chroms → MAF/R2 filter → remove ambiguous SNPs → pgen |
| Already-merged genome-wide VCF | `qc_only` | `MERGED_VCF` | MAF/R2 filter → remove ambiguous SNPs → pgen |
| Already QC'd pgen files | *(skip script)* | `BIM_PREFIX` | Don't submit `Genotype-QC.sh` — just point `BIM_PREFIX` to your files |

> **Note:** PRS-CS/CSx requires a `.bim` file at the `BIM_PREFIX` path. Genotype-QC.sh generates this automatically. Skip-QC users must ensure a `.bim` exists: `plink2 --pfile X --make-just-bim --out X`

### Assumptions

- VCFs have an INFO/R2 field for imputation quality filtering and DS (dosage) format field
- Chromosomes 1-22 (autosomes only)
- Default thresholds: MAF > 0.01, INFO R2 > 0.3, strand-ambiguous SNPs excluded (all configurable)
- If your VCFs lack an R2 field, remove `--extract-if-info` from `Genotype-QC.sh`
- If you have BED/BIM/FAM instead of pgen, convert first: `plink2 --bfile X --make-pgen --out X`

## Unrelated Individuals File

The R QC script optionally filters to unrelated individuals. The file should be:
- Single column, no header
- One individual ID (IID) per line

Set `unrelated_file <- NULL` at the top of `prs-qc.r` to skip this filter.

## Configuration Reference

### config.sh

| Variable | Description |
|----------|-------------|
| `PROJECT_DIR` | Root project directory |
| `GWAS_DIR` | Formatted GWAS summary stats |
| `VCF_DIR` | Raw imputed VCF genotype files |
| `CLEANED_DIR` | QC'd genotype output |
| `WEIGHTS_DIR` | Per-chromosome PRS weight output |
| `COMBINED_DIR` | Genome-wide concatenated weights |
| `SCORES_DIR` | PLINK2 score output |
| `LOGS_DIR` | SLURM log files |
| `BIM_PREFIX` | Final genotype prefix path (without extension); Genotype-QC.sh sets this automatically |
| `TRAITS_PRSCS` | Path to traits_prscs.csv |
| `TRAITS_PRSCSX` | Path to traits_prscsx.csv |
| `PLINK2_EXEC` | PLINK2 binary |
| `BCFTOOLS_EXEC` | bcftools binary |
| `PYTHON_BIN` | Python binary with PRS-CS/CSx deps |
| `PRSCS_DIR` | PRS-CS source directory |
| `PRSCSX_DIR` | PRS-CSx source directory |
| `REF_DIR_PRSCS` | LD reference panels for PRS-CS |
| `REF_DIR_PRSCSX` | LD reference panels for PRS-CSx |
| `MAF_THRESHOLD` | MAF cutoff (default: 0.01) |
| `R2_THRESHOLD` | Imputation R² cutoff (default: 0.3) |
| `SEED` | Random seed (default: 42) |
| `PHI` | Global shrinkage parameter for PRS-CS/CSx (default: `auto`). See [PRScsx docs](https://github.com/getian107/PRScsx) |
| `META_FLAG` | PRS-CSx meta-analysis flag (default: `True`). See [PRScsx docs](https://github.com/getian107/PRScsx) |
| `N_CHR` | Number of chromosomes (default: 22) |
| `VCF_PATTERN` | VCF filename template (default: `chr${chr}_dose.vcf.gz`) |
| `QC_MODE` | Genotype QC mode: `full` or `qc_only` (default: `full`) |
| `MERGED_VCF` | Path to existing merged VCF (only for `qc_only` mode) |
| `SCORE_COLS` | PLINK2 --score column indices (default: `2 4 6`) |

### prs-qc.r (example post-processing)

Edit the four path variables at the top of the script:

| Variable | Description |
|----------|-------------|
| `score_dir` | Directory with .sscore files |
| `weights_dir` | Directory with _allchr.txt weight files |
| `unrelated_file` | Path to unrelated IID list (or `NULL` to skip) |
| `output_path` | Output directory for final dataset + plots |

## Output Files

| File | Step |
|------|------|
| `genome_wide.vcf.gz` | QC (full mode) |
| `genome_wide.bed/bim/fam` | QC |
| `genome_wide.QC.pgen/pvar/psam` | QC |
| `genome_wide.QC.bim` | QC |
| `*_chr*.txt` (per-chrom weights) | PRS-CS/CSx |
| `*_allchr.txt` (genome-wide weights) | Concat |
| `*.sscore` | Apply |
| `prs-final-dataset.txt` | prs-qc.r |
| `prs_distribution.png` | prs-qc.r |
| `prs_correlations.png` | prs-qc.r |
| `snps_per_trait.png` | prs-qc.r |
| `prs_correlation_matrix.txt` | prs-qc.r |

## License

This pipeline is released under the [MIT License](LICENSE). You are free to use, modify, and distribute it with attribution.

All third-party tools called by this pipeline (PRS-CS, PRS-CSx, PLINK2, bcftools, R packages) are subject to their own respective licenses. See the Citations section for links.
