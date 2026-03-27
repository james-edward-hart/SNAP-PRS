#!/usr/bin/env Rscript
# ============================================================
# PRS Score Inspection & Final Dataset
# ============================================================
# 1. PRS distributions and variant missingness
# 2. Filter to unrelated individuals
# 3. Standardize PRS (mean=0, SD=1), build wide-format dataset
# 4. Correlation heatmap
# ============================================================

library(data.table)
library(ggplot2)

# ============================================================
# Configuration — edit these paths for your environment
# ============================================================
score_dir      <- "/path/to/PRS-Scores"          # Directory containing PLINK2 .sscore files
weights_dir    <- "/path/to/Genome-Wide-Weights"  # Directory containing *_allchr.txt weight files
unrelated_file <- "/path/to/unrelated.txt"        # Set to NULL to skip relatedness filter
output_path    <- "/path/to/PRS-Final-Output"     # Output directory for final dataset + plots

dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Input validation
# ============================================================
score_files <- list.files(score_dir, pattern = "\\.sscore$", full.names = TRUE)
if (length(score_files) == 0) stop("No .sscore files found in: ", score_dir)
cat("Found", length(score_files), "score files\n")

weight_files <- list.files(weights_dir, pattern = "_allchr\\.txt$", full.names = TRUE)
if (length(weight_files) == 0) stop("No _allchr.txt files found in: ", weights_dir)

# ============================================================
# Load scores
# ============================================================
scores_list <- lapply(score_files, function(f) {
  dt <- fread(f)
  setnames(dt, gsub("^#", "", names(dt)))
  setnames(dt,
    old = c("NAMED_ALLELE_DOSAGE_SUM", "SCORE1_AVG"),
    new = c("DOSAGE_SUM", "PRS"),
    skip_absent = TRUE
  )
  dt[, TRAIT := gsub("\\.sscore$", "", basename(f))]
})

scores <- rbindlist(scores_list, fill = TRUE)

# ============================================================
# Count variants in weight files
# ============================================================
weight_counts <- rbindlist(lapply(weight_files, function(f) {
  trait <- gsub("_allchr\\.txt$", "", basename(f))
  n <- as.integer(system(paste("wc -l <", shQuote(f)), intern = TRUE))
  data.table(TRAIT = trait, N_WEIGHTS = n)
}))

# ============================================================
# Distribution summary
# ============================================================
cat("\n=== PRS DISTRIBUTION SUMMARY ===\n")

dist_summary <- scores[, .(
  N = .N, MEAN = mean(PRS), SD = sd(PRS),
  MIN = min(PRS), Q25 = quantile(PRS, 0.25),
  MEDIAN = median(PRS), Q75 = quantile(PRS, 0.75), MAX = max(PRS)
), by = TRAIT]

print(dist_summary)

# ============================================================
# Variant missingness
# ============================================================
cat("\n=== VARIANT MISSINGNESS ===\n")

missingness <- scores[, .(
  MEAN_AC = mean(ALLELE_CT), MEDIAN_AC = median(ALLELE_CT),
  MIN_AC = min(ALLELE_CT), MAX_AC = max(ALLELE_CT)
), by = TRAIT]

missingness <- merge(missingness, weight_counts, by = "TRAIT", all.x = TRUE)
missingness[, EXPECTED_AC := N_WEIGHTS * 2]
missingness[, PCT_SCORED := round(100 * MEAN_AC / EXPECTED_AC, 2)]
missingness[, PCT_MISSING := round(100 - PCT_SCORED, 2)]

print(missingness[, .(TRAIT, N_WEIGHTS, EXPECTED_AC, MEAN_AC, PCT_SCORED, PCT_MISSING)])

# ============================================================
# Per-sample missingness
# ============================================================
cat("\n=== PER-SAMPLE MISSINGNESS ===\n")

sample_miss <- merge(
  scores[, .(TRAIT, IID, ALLELE_CT)],
  weight_counts, by = "TRAIT", all.x = TRUE
)
sample_miss[, PCT_SCORED := 100 * ALLELE_CT / (N_WEIGHTS * 2)]

per_sample <- sample_miss[, .(
  MEAN_SCORED = round(mean(PCT_SCORED), 2),
  MIN_SCORED = round(min(PCT_SCORED), 2)
), by = .(IID)]

cat("Sample-level scoring rate:\n")
cat("  Mean across samples:", round(mean(per_sample$MEAN_SCORED), 2), "%\n")
cat("  Worst sample:       ", round(min(per_sample$MIN_SCORED), 2), "%\n")
cat("  Best sample:        ", round(max(per_sample$MEAN_SCORED), 2), "%\n")
cat("  Samples <90% scored:", sum(per_sample$MIN_SCORED < 90), "\n")
cat("  Samples <80% scored:", sum(per_sample$MIN_SCORED < 80), "\n")

# ============================================================
# Plots
# ============================================================
p1 <- ggplot(scores, aes(x = PRS)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~TRAIT, scales = "free") +
  labs(title = "PRS Score Distributions", x = "PRS", y = "Count") +
  theme_bw()
ggsave(file.path(output_path, "prs_distribution.png"), p1, width = 8, height = 5)

p2 <- ggplot(weight_counts, aes(x = reorder(TRAIT, N_WEIGHTS), y = N_WEIGHTS)) +
  geom_col(fill = "darkorange", alpha = 0.7) +
  geom_text(aes(label = scales::comma(N_WEIGHTS)), hjust = -0.1, size = 3) +
  coord_flip() +
  labs(title = "Number of SNPs Used per Trait", x = NULL, y = "N SNPs") +
  theme_bw()
ggsave(file.path(output_path, "snps_per_trait.png"), p2, width = 8, height = 5)

# ============================================================
# Filter to unrelated individuals
# ============================================================
if (!is.null(unrelated_file) && file.exists(unrelated_file)) {
  cat("\n=== FILTERING TO UNRELATED ===\n")
  unrelated <- fread(unrelated_file, header = FALSE)
  setnames(unrelated, old = "V1", new = "IID")

  n_before <- uniqueN(scores$IID)
  scores <- merge(scores, unrelated, by = "IID")
  n_after <- uniqueN(scores$IID)

  cat("Samples before:", n_before, "\n",
      "Unrelated kept:", n_after, "\n",
      "Removed:       ", n_before - n_after, "\n")
} else {
  cat("\nNo unrelated_file specified — using all samples.\n")
}

# ============================================================
# Build final dataset
# ============================================================
cat("\n=== BUILDING FINAL DATASET ===\n")

scores[, PRS := scale(PRS), by = TRAIT]

final <- dcast(scores, IID ~ TRAIT, value.var = "PRS")

cat("Final dataset:", nrow(final), "participants x", ncol(final) - 1, "traits\n")
print(head(final))

# ============================================================
# PRS correlations
# ============================================================
cat("\n=== PRS CORRELATIONS ===\n")

prs_cols <- setdiff(names(final), "IID")
cor_mat <- cor(final[, ..prs_cols], use = "pairwise.complete.obs")

cat("Correlation matrix:\n")
print(round(cor_mat, 3))

cor_melted <- as.data.table(as.table(cor_mat))
setnames(cor_melted, c("Trait1", "Trait2", "r"))

p3 <- ggplot(cor_melted, aes(x = Trait1, y = Trait2, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(r, 2)), size = 3) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, limits = c(-1, 1), name = "r"
  ) +
  labs(title = "PRS Inter-Trait Correlations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank())
ggsave(file.path(output_path, "prs_correlations.png"), p3, width = 10, height = 8)

fwrite(
  as.data.table(cor_mat, keep.rownames = "TRAIT"),
  file.path(output_path, "prs_correlation_matrix.txt"), sep = "\t"
)

# ============================================================
# Save final dataset
# ============================================================
fwrite(final, file.path(output_path, "prs-final-dataset.txt"))
cat("Saved to:", file.path(output_path, "prs-final-dataset.txt"), "\n")
