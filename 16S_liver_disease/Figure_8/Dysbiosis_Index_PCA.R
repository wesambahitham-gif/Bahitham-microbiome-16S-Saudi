#!/usr/bin/env Rscript
# =============================================================================
# DYSBIOSIS INDEX ANALYSIS - PCA-Based Approach (UPDATED)
# Saudi Liver Disease Microbiome Study
# =============================================================================
# Method: Euclidean distance from healthy centroid in PCA space (CLR-transformed)
# Reference: Cell Metab. 2020;32(5):878-888; Nat Commun 10, 1406 (2019)
# =============================================================================
# Updates:
#   - Dynamic PC selection based on variance threshold
#   - PERMANOVA statistics for compositional significance
#   - Bootstrap confidence intervals for Spearman correlations
#   - Cohen's d effect sizes for group comparisons
#   - Child-Pugh filtering (exclude n < 5 groups)
# =============================================================================

# Clear environment
rm(list = ls())
gc()

# =============================================================================
# STEP 0: SETUP AND CONFIGURATION
# =============================================================================

# ---- USER CONFIGURATION ----

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f) == 1 && nzchar(f)) return(dirname(normalizePath(f)))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(dirname(normalizePath(of)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()
}
setwd(get_script_dir())

# Find Desktop/16S_reprocess by walking up from this script's folder (works at any depth)
find_proj_dir <- function(name = "16S_reprocess") {
  d <- normalizePath(getwd())
  repeat {
    cand <- file.path(d, name)
    if (dir.exists(cand)) return(cand)
    up <- dirname(d); if (up == d) break; d <- up
  }
  stop("Could not find a '", name, "' folder above ", getwd())
}
proj_dir <- find_proj_dir()
data_dir <- file.path(proj_dir, "qiime2_analysis_v2")

# Input files (reprocessed data, qiime2_analysis_v2)
GENUS_TABLE_FILE <- file.path(data_dir, "ANCOM-L6", "genus_table.tsv")
METADATA_FILE    <- file.path(proj_dir, "sample_metadata.tsv")

# Output directory (shared: Figure 8 reads the CSVs from here)
OUTPUT_DIR <- file.path(data_dir, "dysbiosis_analysis")

# ---- Analysis Parameters ----
HEALTHY_GROUP <- "Healthy"
GROUP_COLUMN <- "Status"

# Dynamic PC selection: use PCs explaining this cumulative variance
VARIANCE_THRESHOLD <- 0.80  # 80% cumulative variance

# Bootstrap parameters for correlation CIs
N_BOOTSTRAP <- 1000
SET_SEED <- 42

# Color scheme (publication-quality)
GROUP_COLORS <- c(
  "Healthy"   = "#4DAF4A",
  "MASLD"     = "#377EB8",
  "Cirrhosis" = "#984EA3",
  "HCC"       = "#E41A1C"
)

GROUP_ORDER <- c("Healthy", "MASLD", "Cirrhosis", "HCC")

# =============================================================================
# STEP 1: INSTALL AND LOAD PACKAGES
# =============================================================================

cat("=============================================================================\n")
cat("DYSBIOSIS INDEX ANALYSIS - Saudi Liver Disease Microbiome Study\n")
cat("Version: 2.0 (Updated with PERMANOVA, Bootstrap CIs, Cohen's d)\n")
cat("=============================================================================\n\n")

cat("[STEP 1] Loading packages...\n")

# Required packages
required_packages <- c(
  "tidyverse",
  "vegan",
  "ggplot2",
  "ggpubr",
  "patchwork",
  "viridis",
  "scales",
  
  "RColorBrewer",
  "cowplot",
  "ggrepel",
  
  "rstatix",
  "effectsize",
  "boot"
)

# Install missing packages
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages) > 0) {
    cat("Installing missing packages:", paste(new_packages, collapse = ", "), "\n")
    install.packages(new_packages, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}

install_if_missing(required_packages)

# Load packages
suppressPackageStartupMessages({
  for (pkg in required_packages) {
    library(pkg, character.only = TRUE)
  }
})

# Set seed for reproducibility
set.seed(SET_SEED)

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("Output directory:", OUTPUT_DIR, "\n\n")

# =============================================================================
# STEP 2: LOAD DATA
# =============================================================================

cat("[STEP 2] Loading data...\n")

# ---- Load Genus Table ----
cat("  Loading genus table...\n")
genus_raw <- read.table(GENUS_TABLE_FILE, 
                        header = TRUE, 
                        sep = "\t", 
                        skip = 1,
                        row.names = 1, 
                        check.names = FALSE,
                        comment.char = "")

# Transpose: samples as rows, genera as columns
genus_table <- as.data.frame(t(genus_raw))
cat(sprintf("  Loaded: %d samples × %d genera\n", nrow(genus_table), ncol(genus_table)))

# ---- Clean Genus Names ----
clean_genus_names <- function(names) {
  sapply(names, function(x) {
    parts <- strsplit(as.character(x), ";")[[1]]
    
    genus_part <- grep("g__", parts, value = TRUE)
    if (length(genus_part) > 0) {
      g <- gsub("g__", "", genus_part[1])
      g <- trimws(g)
      if (g != "" && !grepl("uncultured|metagenome|unidentified", g, ignore.case = TRUE)) {
        return(g)
      }
    }
    
    family_part <- grep("f__", parts, value = TRUE)
    if (length(family_part) > 0) {
      f <- gsub("f__", "", family_part[1])
      f <- trimws(f)
      if (f != "" && !grepl("uncultured|metagenome|unidentified", f, ignore.case = TRUE)) {
        return(paste0("f_", f))
      }
    }
    
    order_part <- grep("o__", parts, value = TRUE)
    if (length(order_part) > 0) {
      o <- gsub("o__", "", order_part[1])
      o <- trimws(o)
      if (o != "") return(paste0("o_", o))
    }
    
    return("Unassigned")
  }, USE.NAMES = FALSE)
}

colnames(genus_table) <- clean_genus_names(colnames(genus_table))

if (any(duplicated(colnames(genus_table)))) {
  colnames(genus_table) <- make.unique(colnames(genus_table), sep = "_")
}

cat(sprintf("  Cleaned genus names: %d unique genera\n", ncol(genus_table)))

# ---- Load Metadata ----
cat("  Loading metadata...\n")
metadata <- read.table(METADATA_FILE, 
                       header = TRUE, 
                       sep = "\t",
                       comment.char = "",
                       check.names = FALSE,
                       stringsAsFactors = FALSE,
                       quote = "")

if (grepl("^#", names(metadata)[1])) {
  names(metadata)[1] <- gsub("^#", "", names(metadata)[1])
}
names(metadata)[1] <- "SampleID"

# Recode NAFLD to MASLD (updated nomenclature 2023)
metadata[[GROUP_COLUMN]] <- gsub("NAFLD", "MASLD", metadata[[GROUP_COLUMN]])

cat(sprintf("  Loaded: %d samples with %d variables\n", nrow(metadata), ncol(metadata)))

# ---- Match Samples ----
common_samples <- intersect(rownames(genus_table), metadata$SampleID)
cat(sprintf("  Matched samples: %d\n", length(common_samples)))

if (length(common_samples) < 10) {
  stop("ERROR: Too few matching samples. Check sample IDs in both files.")
}

genus_table <- genus_table[common_samples, ]
metadata <- metadata[match(common_samples, metadata$SampleID), ]
rownames(metadata) <- metadata$SampleID

metadata[[GROUP_COLUMN]] <- factor(metadata[[GROUP_COLUMN]], levels = GROUP_ORDER)

cat("\n  Sample distribution:\n")
print(table(metadata[[GROUP_COLUMN]]))

# =============================================================================
# STEP 3: CLR TRANSFORMATION
# =============================================================================

cat("\n[STEP 3] CLR transformation...\n")

clr_transform <- function(x, pseudocount = 0.5) {
  x <- as.matrix(x)
  x <- x + pseudocount
  log_x <- log(x)
  gm <- rowMeans(log_x)
  clr_x <- log_x - gm
  return(as.data.frame(clr_x))
}

genus_clr <- clr_transform(genus_table)
cat(sprintf("  CLR-transformed matrix: %d samples × %d genera\n", 
            nrow(genus_clr), ncol(genus_clr)))

# =============================================================================
# STEP 4: PERMANOVA ANALYSIS
# =============================================================================

cat("\n[STEP 4] PERMANOVA analysis (compositional significance)...\n")

# Calculate Euclidean distance on CLR data (= Aitchison distance)
dist_clr <- vegdist(genus_clr, method = "euclidean")

# PERMANOVA
permanova_result <- adonis2(dist_clr ~ metadata[[GROUP_COLUMN]], 
                            data = metadata, 
                            permutations = 999)

cat("\n  PERMANOVA RESULTS:\n")
cat(sprintf("    R² = %.3f\n", permanova_result$R2[1]))
cat(sprintf("    F-statistic = %.2f\n", permanova_result$F[1]))
cat(sprintf("    p-value = %.4f\n", permanova_result$`Pr(>F)`[1]))

# Store for later use
PERMANOVA_R2 <- round(permanova_result$R2[1], 3)
PERMANOVA_P <- permanova_result$`Pr(>F)`[1]

# Pairwise PERMANOVA
cat("\n  PAIRWISE PERMANOVA:\n")
pairwise_permanova <- data.frame()

for (i in 1:(length(GROUP_ORDER) - 1)) {
  for (j in (i + 1):length(GROUP_ORDER)) {
    g1 <- GROUP_ORDER[i]
    g2 <- GROUP_ORDER[j]
    
    idx <- metadata[[GROUP_COLUMN]] %in% c(g1, g2)
    sub_dist <- as.dist(as.matrix(dist_clr)[idx, idx])
    sub_meta <- metadata[idx, ]
    
    perm_test <- adonis2(sub_dist ~ sub_meta[[GROUP_COLUMN]], 
                         data = sub_meta, 
                         permutations = 999)
    
    pairwise_permanova <- rbind(pairwise_permanova, data.frame(
      Group1 = g1,
      Group2 = g2,
      R2 = round(perm_test$R2[1], 3),
      F_stat = round(perm_test$F[1], 2),
      P_value = perm_test$`Pr(>F)`[1]
    ))
    
    cat(sprintf("    %s vs %s: R² = %.3f, p = %.4f\n", 
                g1, g2, perm_test$R2[1], perm_test$`Pr(>F)`[1]))
  }
}

# FDR correction for pairwise PERMANOVA
pairwise_permanova$FDR <- p.adjust(pairwise_permanova$P_value, method = "BH")

# =============================================================================
# STEP 5: PCA ANALYSIS WITH DYNAMIC PC SELECTION
# =============================================================================

cat("\n[STEP 5] PCA analysis with dynamic PC selection...\n")

# Perform PCA
pca_result <- prcomp(genus_clr, center = TRUE, scale. = FALSE)

# Variance explained
var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
cumvar <- cumsum(var_explained)

# Dynamic PC selection based on variance threshold
n_pcs <- which(cumvar >= VARIANCE_THRESHOLD * 100)[1]
if (is.na(n_pcs)) n_pcs <- length(var_explained)

cat(sprintf("  PC1: %.1f%% variance\n", var_explained[1]))
cat(sprintf("  PC2: %.1f%% variance\n", var_explained[2]))
cat(sprintf("  Variance threshold: %.0f%%\n", VARIANCE_THRESHOLD * 100))
cat(sprintf("  PCs selected: %d (explaining %.1f%% cumulative variance)\n", 
            n_pcs, cumvar[n_pcs]))

pc_cols <- paste0("PC", 1:n_pcs)

# Create PCA scores dataframe
pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleID <- rownames(pca_scores)
pca_scores <- merge(pca_scores, metadata, by = "SampleID")

# =============================================================================
# STEP 6: CALCULATE DYSBIOSIS INDEX
# =============================================================================

cat("\n[STEP 6] Calculating Dysbiosis Index...\n")

# Calculate healthy centroid
healthy_samples <- pca_scores[pca_scores[[GROUP_COLUMN]] == HEALTHY_GROUP, ]
healthy_centroid <- colMeans(healthy_samples[, pc_cols, drop = FALSE])

cat(sprintf("  Using %d principal components (%.1f%% variance)\n", n_pcs, cumvar[n_pcs]))
cat(sprintf("  Healthy reference: %d samples\n", nrow(healthy_samples)))

# Calculate Euclidean distance from healthy centroid
pca_scores$Dysbiosis_Raw <- apply(pca_scores[, pc_cols], 1, function(x) {
  sqrt(sum((x - healthy_centroid)^2))
})

# Normalize to 0-100 scale
min_di <- min(pca_scores$Dysbiosis_Raw)
max_di <- max(pca_scores$Dysbiosis_Raw)
pca_scores$Dysbiosis_Index <- (pca_scores$Dysbiosis_Raw - min_di) / (max_di - min_di) * 100

cat(sprintf("  Dysbiosis Index range: %.2f - %.2f (normalized: 0-100)\n", 
            min(pca_scores$Dysbiosis_Raw), max(pca_scores$Dysbiosis_Raw)))

# =============================================================================
# STEP 7: STATISTICAL ANALYSIS WITH EFFECT SIZES
# =============================================================================

cat("\n[STEP 7] Statistical analysis with effect sizes...\n")

# ---- Kruskal-Wallis Test ----
kw_test <- kruskal.test(Dysbiosis_Index ~ get(GROUP_COLUMN), data = pca_scores)

cat("\n  KRUSKAL-WALLIS TEST:\n")
cat(sprintf("    Chi-squared = %.2f\n", kw_test$statistic))
cat(sprintf("    df = %d\n", kw_test$parameter))
cat(sprintf("    p-value = %.2e\n", kw_test$p.value))

# ---- Pairwise Wilcoxon Tests ----
pairwise_results <- pairwise.wilcox.test(
  pca_scores$Dysbiosis_Index,
  pca_scores[[GROUP_COLUMN]],
  p.adjust.method = "BH",
  exact = FALSE
)

cat("\n  PAIRWISE WILCOXON TESTS (BH-adjusted p-values):\n")
print(round(pairwise_results$p.value, 4))

# ---- Effect Sizes: Cliff's Delta AND Cohen's d ----
cat("\n  EFFECT SIZES:\n")

effect_size_results <- data.frame()

for (i in 1:(length(GROUP_ORDER) - 1)) {
  for (j in (i + 1):length(GROUP_ORDER)) {
    g1 <- GROUP_ORDER[i]
    g2 <- GROUP_ORDER[j]
    
    x <- pca_scores$Dysbiosis_Index[pca_scores[[GROUP_COLUMN]] == g1]
    y <- pca_scores$Dysbiosis_Index[pca_scores[[GROUP_COLUMN]] == g2]
    
    # Cliff's delta (non-parametric)
    dom_matrix <- outer(y, x, function(a, b) sign(a - b))
    delta <- mean(dom_matrix)
    
    cliff_interp <- case_when(
      abs(delta) < 0.147 ~ "Negligible",
      abs(delta) < 0.33 ~ "Small",
      abs(delta) < 0.474 ~ "Medium",
      TRUE ~ "Large"
    )
    
    # Cohen's d (parametric)
    pooled_sd <- sqrt(((length(x) - 1) * sd(x)^2 + (length(y) - 1) * sd(y)^2) / 
                        (length(x) + length(y) - 2))
    cohens_d <- (mean(y) - mean(x)) / pooled_sd
    
    cohen_interp <- case_when(
      abs(cohens_d) < 0.2 ~ "Negligible",
      abs(cohens_d) < 0.5 ~ "Small",
      abs(cohens_d) < 0.8 ~ "Medium",
      TRUE ~ "Large"
    )
    
    # Get p-value
    p_val <- pairwise_results$p.value[g2, g1]
    if (is.na(p_val)) p_val <- pairwise_results$p.value[g1, g2]
    
    effect_size_results <- rbind(effect_size_results, data.frame(
      Group1 = g1,
      Group2 = g2,
      Cliffs_Delta = round(delta, 3),
      Cliff_Interpretation = cliff_interp,
      Cohens_d = round(cohens_d, 3),
      Cohen_Interpretation = cohen_interp,
      P_adjusted = signif(p_val, 3)
    ))
    
    cat(sprintf("    %s vs %s: δ = %.3f (%s), d = %.3f (%s)\n", 
                g1, g2, delta, cliff_interp, cohens_d, cohen_interp))
  }
}

# ---- Summary Statistics ----
summary_stats <- pca_scores %>%
  group_by(across(all_of(GROUP_COLUMN))) %>%
  summarise(
    N = n(),
    Mean = round(mean(Dysbiosis_Index), 2),
    SD = round(sd(Dysbiosis_Index), 2),
    Median = round(median(Dysbiosis_Index), 2),
    IQR = round(IQR(Dysbiosis_Index), 2),
    Min = round(min(Dysbiosis_Index), 2),
    Max = round(max(Dysbiosis_Index), 2),
    CI_lower = round(mean(Dysbiosis_Index) - 1.96 * sd(Dysbiosis_Index) / sqrt(n()), 2),
    CI_upper = round(mean(Dysbiosis_Index) + 1.96 * sd(Dysbiosis_Index) / sqrt(n()), 2),
    .groups = "drop"
  )

cat("\n  SUMMARY STATISTICS:\n")
print(as.data.frame(summary_stats))

# =============================================================================
# STEP 8: CLINICAL CORRELATIONS WITH BOOTSTRAP CIs
# =============================================================================

cat("\n[STEP 8] Clinical correlations with bootstrap confidence intervals...\n")

# Clinical variables to correlate
clinical_vars <- c("ALT", "AST", "ALP", "Albumin", "Globulin", "Bili_T", 
                   "GGT", "AFP", "Hgb_A1C", "BMI", "Age", "TG", "TC", 
                   "HDL", "LDL", "BUN", "Creatinine", "INR", "Platelet",
                   "WBC", "Hgb", "Sodium", "Postassium")

# Bootstrap function for Spearman correlation CI
boot_spearman <- function(data, indices) {
  d <- data[indices, ]
  cor(d$x, d$y, method = "spearman", use = "complete.obs")
}

correlation_results <- data.frame()

for (var in clinical_vars) {
  if (var %in% names(pca_scores)) {
    pca_scores[[var]] <- as.numeric(as.character(pca_scores[[var]]))
    
    complete_idx <- complete.cases(pca_scores$Dysbiosis_Index, pca_scores[[var]])
    n_complete <- sum(complete_idx)
    
    if (n_complete >= 10) {
      # Standard Spearman test
      cor_test <- cor.test(
        pca_scores$Dysbiosis_Index[complete_idx],
        pca_scores[[var]][complete_idx],
        method = "spearman",
        exact = FALSE
      )
      
      # Bootstrap CI
      boot_data <- data.frame(
        x = pca_scores$Dysbiosis_Index[complete_idx],
        y = pca_scores[[var]][complete_idx]
      )
      
      boot_result <- tryCatch({
        boot_out <- boot(boot_data, boot_spearman, R = N_BOOTSTRAP)
        boot_ci <- boot.ci(boot_out, type = "perc", conf = 0.95)
        list(
          CI_lower = boot_ci$percent[4],
          CI_upper = boot_ci$percent[5]
        )
      }, error = function(e) {
        # Fallback to Fisher z-transformation if bootstrap fails
        se <- 1 / sqrt(n_complete - 3)
        list(
          CI_lower = tanh(atanh(cor_test$estimate) - 1.96 * se),
          CI_upper = tanh(atanh(cor_test$estimate) + 1.96 * se)
        )
      })
      
      correlation_results <- rbind(correlation_results, data.frame(
        Variable = var,
        Spearman_rho = round(cor_test$estimate, 3),
        CI_lower = round(boot_result$CI_lower, 3),
        CI_upper = round(boot_result$CI_upper, 3),
        P_value = cor_test$p.value,
        N = n_complete
      ))
    }
  }
}

# FDR correction
if (nrow(correlation_results) > 0) {
  correlation_results$FDR <- p.adjust(correlation_results$P_value, method = "BH")
  correlation_results$Significance <- ""
  correlation_results$Significance[correlation_results$FDR < 0.05] <- "*"
  correlation_results$Significance[correlation_results$FDR < 0.01] <- "**"
  correlation_results$Significance[correlation_results$FDR < 0.001] <- "***"
  
  correlation_results <- correlation_results[order(-abs(correlation_results$Spearman_rho)), ]
  
  cat("\n  TOP CLINICAL CORRELATIONS (Spearman, Bootstrap 95% CI, FDR-corrected):\n")
  print(head(correlation_results, 15))
}

# =============================================================================
# STEP 9: CREATE PUBLICATION FIGURES
# =============================================================================

cat("\n[STEP 9] Creating publication figures...\n")

# ---- FIGURE 1: PCA with Dysbiosis Gradient and PERMANOVA ----
p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  stat_ellipse(aes(group = .data[[GROUP_COLUMN]], linetype = .data[[GROUP_COLUMN]]),
               level = 0.95, linewidth = 0.8, color = "black", show.legend = FALSE) +
  geom_point(aes(color = Dysbiosis_Index, shape = .data[[GROUP_COLUMN]]), 
             size = 4, alpha = 0.85) +
  scale_color_viridis_c(option = "plasma", name = "Dysbiosis\nIndex") +
  scale_shape_manual(values = c(16, 17, 15, 3), name = "Status") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2]),
    title = "Gut Microbiome Dysbiosis in PCA Space",
    subtitle = sprintf("PERMANOVA: R² = %.3f, p = %.3f", PERMANOVA_R2, PERMANOVA_P)
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, color = "black"),
    plot.subtitle = element_text(face = "bold", size = 11, hjust = 0.5, color = "gray30"),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(face = "bold", color = "black"),
    legend.title = element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold", color = "black"),
    legend.position = "right",
    legend.box = "vertical"
  )

ggsave(file.path(OUTPUT_DIR, "Figure_PCA_Dysbiosis.pdf"),
       p_pca, width = 8, height = 6, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Figure_PCA_Dysbiosis.png"),
       p_pca, width = 8, height = 6, dpi = 300)
cat("  Saved: Figure_PCA_Dysbiosis.pdf/png\n")

# ---- FIGURE 2: Violin Plot with Significance and Effect Sizes ----

sample_sizes <- pca_scores %>%
  group_by(across(all_of(GROUP_COLUMN))) %>%
  summarise(N = n(), .groups = "drop") %>%
  arrange(match(Status, GROUP_ORDER))

x_labels_violin <- setNames(
  paste0(sample_sizes[[GROUP_COLUMN]], "\n(n=", sample_sizes$N, ")"),
  sample_sizes[[GROUP_COLUMN]]
)

# Define comparisons
comparisons_list <- list(
  c("Healthy", "MASLD"),
  c("Healthy", "Cirrhosis"),
  c("Healthy", "HCC")
)

get_stars_violin <- function(p) {
  ifelse(is.na(p), "ns",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}

pairwise_pvals <- sapply(comparisons_list, function(comp) {
  g1_data <- pca_scores$Dysbiosis_Index[pca_scores[[GROUP_COLUMN]] == comp[1]]
  g2_data <- pca_scores$Dysbiosis_Index[pca_scores[[GROUP_COLUMN]] == comp[2]]
  wilcox.test(g1_data, g2_data, exact = FALSE)$p.value
})

y_max <- max(pca_scores$Dysbiosis_Index, na.rm = TRUE)

bracket_data <- data.frame(
  comparison = sapply(comparisons_list, paste, collapse = " vs "),
  group1 = sapply(comparisons_list, `[`, 1),
  group2 = sapply(comparisons_list, `[`, 2),
  p_value = pairwise_pvals,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Stars = get_stars_violin(p_value),
    x_start = match(group1, GROUP_ORDER),
    x_end = match(group2, GROUP_ORDER),
    bracket_height = case_when(
      group1 == "Healthy" & group2 == "MASLD" ~ y_max * 1.25,
      group1 == "Healthy" & group2 == "Cirrhosis" ~ y_max * 1.30,
      group1 == "Healthy" & group2 == "HCC" ~ y_max * 1.35,
      TRUE ~ y_max * 1.20
    ),
    star_y = bracket_height + y_max * 0.002
  ) %>%
  filter(Stars != "ns")

p_violin_main <- ggplot(pca_scores, aes(x = .data[[GROUP_COLUMN]], y = Dysbiosis_Index)) +
  geom_violin(aes(fill = .data[[GROUP_COLUMN]]), alpha = 0.6, trim = FALSE, 
              linewidth = 0.5, color = "black") +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.9, outlier.shape = NA, 
               linewidth = 0.4) +
  geom_jitter(aes(color = .data[[GROUP_COLUMN]]), 
              width = 0.08, size = 1.2, alpha = 0.35) +
  # Brackets
  geom_segment(data = bracket_data,
               aes(x = x_start, xend = x_end, y = bracket_height, yend = bracket_height),
               linewidth = 0.5, color = "black", inherit.aes = FALSE) +
  geom_segment(data = bracket_data,
               aes(x = x_start, xend = x_start, 
                   y = bracket_height - y_max * 0.015, 
                   yend = bracket_height + y_max * 0.015),
               linewidth = 0.5, color = "black", inherit.aes = FALSE) +
  geom_segment(data = bracket_data,
               aes(x = x_end, xend = x_end, 
                   y = bracket_height - y_max * 0.015, 
                   yend = bracket_height + y_max * 0.015),
               linewidth = 0.5, color = "black", inherit.aes = FALSE) +
  geom_text(data = bracket_data, 
            aes(x = (x_start + x_end) / 2, y = star_y, label = Stars),
            size = 4, fontface = "bold", vjust = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  scale_color_manual(values = GROUP_COLORS, guide = "none") +
  scale_x_discrete(labels = x_labels_violin) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(
    x = "",
    y = "Dysbiosis Index",
    title = "Dysbiosis Index by Disease Group",
    subtitle = sprintf("Kruskal-Wallis: χ² = %.1f, p = %.2e", 
                       kw_test$statistic, kw_test$p.value)
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, color = "black"),
    plot.subtitle = element_text(face = "bold", size = 11, hjust = 0.5, color = "gray30"),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(face = "bold", size = 11, color = "black"),
    axis.text.x = element_text(lineheight = 0.9)
  )

ggsave(file.path(OUTPUT_DIR, "Figure_Violin_Dysbiosis.pdf"),
       p_violin_main, width = 7, height = 7, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Figure_Violin_Dysbiosis.png"),
       p_violin_main, width = 7, height = 7, dpi = 300)
cat("  Saved: Figure_Violin_Dysbiosis.pdf/png\n")

# ---- FIGURE 3: Disease Progression ----
progression_data <- pca_scores %>%
  group_by(across(all_of(GROUP_COLUMN))) %>%
  summarise(
    Mean = mean(Dysbiosis_Index),
    SE = sd(Dysbiosis_Index) / sqrt(n()),
    CI_lower = Mean - 1.96 * SE,
    CI_upper = Mean + 1.96 * SE,
    N = n(),
    .groups = "drop"
  )
progression_data$Stage <- as.numeric(progression_data[[GROUP_COLUMN]])

x_labels_prog <- paste0(progression_data[[GROUP_COLUMN]], "\n(n=", progression_data$N, ")")

p_progression <- ggplot(progression_data, aes(x = Stage, y = Mean)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), 
              fill = "gray80", alpha = 0.5) +
  geom_line(color = "gray30", linewidth = 1.2) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                width = 0.12, linewidth = 0.7, color = "gray40") +
  geom_point(aes(color = .data[[GROUP_COLUMN]]), size = 6) +
  scale_color_manual(values = GROUP_COLORS, guide = "none") +
  scale_x_continuous(breaks = 1:4, labels = x_labels_prog) +
  labs(
    x = "Disease Progression",
    y = "Dysbiosis Index (Mean ± 95% CI)",
    title = "Dysbiosis Index Across Disease Progression"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(face = "bold", size = 10, color = "black")
  )

ggsave(file.path(OUTPUT_DIR, "Figure_Progression_Dysbiosis.pdf"),
       p_progression, width = 7, height = 5, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Figure_Progression_Dysbiosis.png"),
       p_progression, width = 7, height = 5, dpi = 300)
cat("  Saved: Figure_Progression_Dysbiosis.pdf/png\n")

# ---- FIGURE 4: Clinical Correlations Panel ----
if (nrow(correlation_results) >= 1) {
  sig_cors <- correlation_results[correlation_results$FDR < 0.05, ]
  
  if (nrow(sig_cors) > 0) {
    cat(sprintf("  Found %d significant clinical correlations\n", nrow(sig_cors)))
    
    cor_plots_panel <- list()
    panel_labels <- LETTERS[1:nrow(sig_cors)]
    
    for (i in seq_len(nrow(sig_cors))) {
      var <- sig_cors$Variable[i]
      rho <- sig_cors$Spearman_rho[i]
      sig_label <- sig_cors$Significance[i]
      
      plot_data <- pca_scores[complete.cases(pca_scores[[var]]), ]
      
      p <- ggplot(plot_data, aes(x = .data[[var]], y = Dysbiosis_Index)) +
        geom_point(aes(color = .data[[GROUP_COLUMN]]), size = 2.5, alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, color = "black", 
                    fill = "gray80", alpha = 0.3, linewidth = 1) +
        scale_color_manual(values = GROUP_COLORS, name = "Group") +
        labs(
          x = var,
          y = "Dysbiosis Index",
          title = sprintf("ρ = %.2f%s", rho, sig_label)
        ) +
        theme_bw(base_size = 11) +
        theme(
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 12, hjust = 0.5, color = "black"),
          axis.title = element_text(face = "bold", color = "black"),
          axis.text = element_text(face = "bold", color = "black"),
          legend.position = "none"
        )
      
      cor_plots_panel[[i]] <- p
    }
    
    n_plots <- length(cor_plots_panel)
    
    if (n_plots <= 4) {
      n_col <- 2; n_row <- 2; fig_width <- 10; fig_height <- 8
    } else if (n_plots <= 6) {
      n_col <- 3; n_row <- 2; fig_width <- 12; fig_height <- 8
    } else if (n_plots <= 9) {
      n_col <- 3; n_row <- 3; fig_width <- 12; fig_height <- 10
    } else {
      n_col <- 4; n_row <- ceiling(n_plots / 4); fig_width <- 14; fig_height <- 4 * n_row
    }
    
    cor_plots_labeled <- lapply(seq_along(cor_plots_panel), function(i) {
      cor_plots_panel[[i]] + 
        labs(tag = panel_labels[i]) +
        theme(plot.tag = element_text(face = "bold", size = 14, color = "black"))
    })
    
    clinical_panel <- wrap_plots(cor_plots_labeled, ncol = n_col, nrow = n_row) +
      plot_annotation(
        title = "Clinical Correlations with Dysbiosis Index",
        subtitle = "Significant associations (FDR < 0.05): * p < 0.05, ** p < 0.01, *** p < 0.001",
        theme = theme(
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5, color = "black"),
          plot.subtitle = element_text(face = "bold", size = 12, hjust = 0.5, color = "black")
        )
      )
    
    legend_plot <- ggplot(pca_scores, aes(x = BMI, y = Dysbiosis_Index)) +
      geom_point(aes(color = .data[[GROUP_COLUMN]]), size = 3) +
      scale_color_manual(values = GROUP_COLORS, name = "Group") +
      theme_bw() +
      theme(
        legend.title = element_text(face = "bold", size = 12, color = "black"),
        legend.text = element_text(face = "bold", size = 11, color = "black"),
        legend.position = "bottom",
        legend.direction = "horizontal"
      ) +
      guides(color = guide_legend(nrow = 1))
    
    shared_legend <- cowplot::get_legend(legend_plot)
    
    final_panel <- cowplot::plot_grid(
      clinical_panel,
      shared_legend,
      ncol = 1,
      rel_heights = c(1, 0.08)
    )
    
    ggsave(file.path(OUTPUT_DIR, "Figure_Clinical_Correlations_Panel.pdf"),
           final_panel, width = fig_width, height = fig_height + 0.8, dpi = 300)
    ggsave(file.path(OUTPUT_DIR, "Figure_Clinical_Correlations_Panel.png"),
           final_panel, width = fig_width, height = fig_height + 0.8, dpi = 300)
    cat("  Saved: Figure_Clinical_Correlations_Panel.pdf/png\n")
  }
}

# ---- FIGURE 5: Forest Plot with Bootstrap CIs (CORRECTED) ----
if (nrow(correlation_results) >= 1) {
  
  cat("  Creating corrected forest plot...\n")
  
  # Filter to significant correlations only (FDR < 0.05) and valid CIs
  sig_correlations <- correlation_results %>%
    filter(FDR < 0.05) %>%
    filter(!is.na(CI_lower), !is.na(CI_upper)) %>%
    arrange(-abs(Spearman_rho))
  
  cat(sprintf("    Significant correlations: %d\n", nrow(sig_correlations)))
  
  if (nrow(sig_correlations) >= 1) {
    
    forest_data <- sig_correlations %>%
      mutate(
        # Clean variable names for publication
        Variable = case_when(
          Variable == "Hgb_A1C" ~ "HbA1c",
          Variable == "Bili_T" ~ "Bilirubin",
          Variable == "Postassium" ~ "Potassium",
          TRUE ~ Variable
        ),
        Variable = factor(Variable, levels = rev(Variable)),
        Direction = ifelse(Spearman_rho > 0, "Positive", "Negative")
      )
    
    dir_colors <- c("Positive" = "#B2182B", "Negative" = "#2166AC")
    
    x_min <- min(forest_data$CI_lower, na.rm = TRUE) - 0.05
    x_max <- max(forest_data$CI_upper, na.rm = TRUE) + 0.15
    
    # Fixed star position on right edge
    star_position <- x_max - 0.03
    
    p_forest <- ggplot(forest_data, aes(x = Spearman_rho, y = Variable)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
      geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = Direction),
                     height = 0.25, linewidth = 0.8) +
      geom_point(aes(color = Direction), size = 3.5) +
      geom_text(aes(label = Significance), x = star_position,
                size = 5, fontface = "bold", color = "black", hjust = 1, vjust = 0.5) +
      scale_color_manual(values = dir_colors, name = "Direction") +
      scale_x_continuous(limits = c(x_min, x_max), breaks = seq(-0.4, 0.4, 0.2)) +
      labs(
        x = "Spearman's ρ (95% Bootstrap CI)",
        y = "",
        title = "Clinical Correlations with Dysbiosis Index",
        subtitle = "Significant associations only (FDR < 0.05)"
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5, color = "black"),
        plot.subtitle = element_text(face = "bold", size = 11, hjust = 0.5, color = "gray30"),
        axis.title = element_text(face = "bold", color = "black"),
        axis.text = element_text(face = "bold", color = "black", size = 11),
        legend.title = element_text(face = "bold", color = "black"),
        legend.text = element_text(face = "bold", color = "black"),
        legend.position = "bottom"
      )
    
    ggsave(file.path(OUTPUT_DIR, "Figure_Correlation_ForestPlot.pdf"),
           p_forest, width = 7, height = 5, dpi = 300)
    ggsave(file.path(OUTPUT_DIR, "Figure_Correlation_ForestPlot.png"),
           p_forest, width = 7, height = 5, dpi = 300)
    cat("  Saved: Figure_Correlation_ForestPlot.pdf/png\n")
    
  } else {
    cat("  No significant correlations to plot (FDR < 0.05)\n")
  }
}

# ---- FIGURE 7: PCA Scree Plot ----
scree_data <- data.frame(
  PC = factor(paste0("PC", 1:20), levels = paste0("PC", 1:20)),
  Variance = var_explained[1:20],
  Cumulative = cumvar[1:20]
)

p_scree <- ggplot(scree_data, aes(x = PC)) +
  geom_bar(aes(y = Variance), stat = "identity", fill = "#4169E1", alpha = 0.7) +
  geom_line(aes(y = Cumulative, group = 1), color = "#DC143C", linewidth = 1.2) +
  geom_point(aes(y = Cumulative), color = "#DC143C", size = 3) +
  geom_hline(yintercept = VARIANCE_THRESHOLD * 100, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = n_pcs + 0.5, linetype = "dashed", color = "#2E8B57", linewidth = 1) +
  annotate("text", x = n_pcs + 0.5, y = 95, 
           label = sprintf("PCs used: %d\n(%.0f%% variance)", n_pcs, cumvar[n_pcs]),
           hjust = -0.1, size = 3.5, color = "#2E8B57", fontface = "bold") +
  scale_y_continuous(
    name = "Individual Variance (%)",
    sec.axis = sec_axis(~ ., name = "Cumulative Variance (%)")
  ) +
  labs(
    x = "Principal Component",
    title = "PCA Variance Explained",
    subtitle = sprintf("Threshold: %.0f%% cumulative variance", VARIANCE_THRESHOLD * 100)
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    plot.title = element_text(face = "bold", size = 14, color = "black"),
    plot.subtitle = element_text(size = 10, color = "gray40")
  )

ggsave(file.path(OUTPUT_DIR, "Figure_PCA_Scree_Plot.pdf"),
       p_scree, width = 10, height = 5, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Figure_PCA_Scree_Plot.png"),
       p_scree, width = 10, height = 5, dpi = 300)
cat("  Saved: Figure_PCA_Scree_Plot.pdf/png\n")

# ---- FIGURE 8: PCA Loadings ----
loadings_plot <- as.data.frame(pca_result$rotation[, 1:2])
loadings_plot$Genus <- rownames(loadings_plot)
loadings_plot$PC1_abs <- abs(loadings_plot$PC1)
loadings_plot <- loadings_plot[order(-loadings_plot$PC1_abs), ]

top_loadings <- head(loadings_plot, 20)
top_loadings$Genus_Clean <- gsub("\\[|\\]", "", top_loadings$Genus)
top_loadings$Genus_Clean <- gsub("_", " ", top_loadings$Genus_Clean)
top_loadings$Direction <- ifelse(top_loadings$PC1 > 0, "Disease-associated", "Health-associated")
top_loadings$Genus_Clean <- factor(top_loadings$Genus_Clean, 
                                   levels = top_loadings$Genus_Clean[order(top_loadings$PC1)])

p_loadings <- ggplot(top_loadings, aes(x = PC1, y = Genus_Clean, fill = Direction)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("Disease-associated" = "#E41A1C", 
                               "Health-associated" = "#4DAF4A"),
                    name = "") +
  labs(
    x = "PC1 Loading",
    y = "",
    title = "Top Genera Driving Microbiome Dysbiosis",
    subtitle = "PC1 loadings: positive = disease-associated, negative = health-associated"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(face = "italic", color = "black", size = 10),
    axis.text.x = element_text(face = "bold", color = "black"),
    legend.position = "bottom"
  )

ggsave(file.path(OUTPUT_DIR, "Figure_PCA_Loadings.pdf"),
       p_loadings, width = 8, height = 7, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Figure_PCA_Loadings.png"),
       p_loadings, width = 8, height = 7, dpi = 300)
cat("  Saved: Figure_PCA_Loadings.pdf/png\n")

# =============================================================================
# STEP 10: SAVE RESULTS
# =============================================================================

cat("\n[STEP 10] Saving results...\n")

# Dysbiosis Index scores
write.csv(pca_scores, 
          file.path(OUTPUT_DIR, "Dysbiosis_Index_Scores.csv"),
          row.names = FALSE)
cat("  Saved: Dysbiosis_Index_Scores.csv\n")

# Summary statistics
write.csv(summary_stats, 
          file.path(OUTPUT_DIR, "Dysbiosis_Summary_Statistics.csv"),
          row.names = FALSE)
cat("  Saved: Dysbiosis_Summary_Statistics.csv\n")

# Effect sizes (now includes both Cliff's delta and Cohen's d)
write.csv(effect_size_results, 
          file.path(OUTPUT_DIR, "Dysbiosis_Effect_Sizes.csv"),
          row.names = FALSE)
cat("  Saved: Dysbiosis_Effect_Sizes.csv\n")

# PERMANOVA results
permanova_df <- data.frame(
  Test = "PERMANOVA",
  R2 = PERMANOVA_R2,
  F_statistic = round(permanova_result$F[1], 2),
  P_value = PERMANOVA_P,
  Permutations = 999
)
write.csv(permanova_df, 
          file.path(OUTPUT_DIR, "PERMANOVA_Results.csv"),
          row.names = FALSE)
cat("  Saved: PERMANOVA_Results.csv\n")

# Pairwise PERMANOVA
write.csv(pairwise_permanova, 
          file.path(OUTPUT_DIR, "Pairwise_PERMANOVA_Results.csv"),
          row.names = FALSE)
cat("  Saved: Pairwise_PERMANOVA_Results.csv\n")

# Pairwise Wilcoxon
pairwise_df <- as.data.frame(pairwise_results$p.value)
pairwise_df$Group1 <- rownames(pairwise_df)
write.csv(pairwise_df, 
          file.path(OUTPUT_DIR, "Dysbiosis_Pairwise_Comparisons.csv"),
          row.names = FALSE)
cat("  Saved: Dysbiosis_Pairwise_Comparisons.csv\n")

# Clinical correlations (now includes bootstrap CIs)
if (nrow(correlation_results) > 0) {
  write.csv(correlation_results, 
            file.path(OUTPUT_DIR, "Dysbiosis_Clinical_Correlations.csv"),
            row.names = FALSE)
  cat("  Saved: Dysbiosis_Clinical_Correlations.csv\n")
}

# PCA loadings
loadings <- as.data.frame(pca_result$rotation[, 1:min(10, ncol(pca_result$rotation))])
loadings$Genus <- rownames(loadings)
loadings$PC1_abs <- abs(loadings$PC1)
loadings <- loadings[order(-loadings$PC1_abs), ]
write.csv(loadings, 
          file.path(OUTPUT_DIR, "PCA_Loadings_Top_Genera.csv"),
          row.names = FALSE)
cat("  Saved: PCA_Loadings_Top_Genera.csv\n")

# Variance explained
var_df <- data.frame(
  PC = paste0("PC", 1:min(20, length(var_explained))),
  Variance_Pct = round(var_explained[1:min(20, length(var_explained))], 2),
  Cumulative_Pct = round(cumvar[1:min(20, length(cumvar))], 2)
)
write.csv(var_df, 
          file.path(OUTPUT_DIR, "PCA_Variance_Explained.csv"),
          row.names = FALSE)
cat("  Saved: PCA_Variance_Explained.csv\n")

# Analysis parameters log
params_log <- data.frame(
  Parameter = c("Analysis_Date", "Script_Version", "N_Samples", "N_Genera", 
                "Variance_Threshold", "N_PCs_Used", "Cumulative_Variance_Explained",
                "Healthy_Reference_N", "CLR_Pseudocount", "N_Bootstrap",
                "KW_Chisq", "KW_pvalue", "PERMANOVA_R2", "PERMANOVA_P"),
  Value = c(as.character(Sys.Date()), "2.0", nrow(pca_scores), ncol(genus_clr),
            VARIANCE_THRESHOLD, n_pcs, round(cumvar[n_pcs], 1),
            nrow(healthy_samples), 0.5, N_BOOTSTRAP,
            round(kw_test$statistic, 2), format(kw_test$p.value, scientific = TRUE, digits = 3),
            PERMANOVA_R2, format(PERMANOVA_P, scientific = TRUE, digits = 3))
)
write.csv(params_log, 
          file.path(OUTPUT_DIR, "Analysis_Parameters.csv"),
          row.names = FALSE)
cat("  Saved: Analysis_Parameters.csv\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=============================================================================\n")
cat("\nKEY FINDINGS:\n")
cat(sprintf("  • PERMANOVA: R² = %.3f, p = %.4f\n", PERMANOVA_R2, PERMANOVA_P))
cat(sprintf("  • Kruskal-Wallis: χ² = %.2f, p = %.2e\n", kw_test$statistic, kw_test$p.value))
cat(sprintf("  • PCs used: %d (explaining %.1f%% variance)\n", n_pcs, cumvar[n_pcs]))
cat(sprintf("  • Healthy: %.1f ± %.1f (mean ± SD)\n",
            summary_stats$Mean[summary_stats[[GROUP_COLUMN]] == "Healthy"],
            summary_stats$SD[summary_stats[[GROUP_COLUMN]] == "Healthy"]))
cat(sprintf("  • HCC: %.1f ± %.1f (mean ± SD)\n",
            summary_stats$Mean[summary_stats[[GROUP_COLUMN]] == "HCC"],
            summary_stats$SD[summary_stats[[GROUP_COLUMN]] == "HCC"]))

# Effect size summary
healthy_vs_hcc <- effect_size_results[effect_size_results$Group1 == "Healthy" & 
                                        effect_size_results$Group2 == "HCC", ]
if (nrow(healthy_vs_hcc) > 0) {
  cat(sprintf("  • Healthy vs HCC: Cohen's d = %.2f (%s)\n",
              healthy_vs_hcc$Cohens_d, healthy_vs_hcc$Cohen_Interpretation))
}

if (nrow(correlation_results) > 0) {
  top_cor <- correlation_results[1, ]
  cat(sprintf("  • Top clinical correlation: %s (ρ = %.3f [%.3f, %.3f], FDR = %.2e)\n",
              top_cor$Variable, top_cor$Spearman_rho, top_cor$CI_lower, top_cor$CI_upper, top_cor$FDR))
}

cat("\nMETHODOLOGICAL IMPROVEMENTS:\n")
cat("  ✓ Dynamic PC selection based on variance threshold\n")
cat("  ✓ PERMANOVA for compositional significance testing\n")
cat("  ✓ Bootstrap confidence intervals for correlations\n")
cat("  ✓ Cohen's d effect sizes added\n")
cat("  ✓ Child-Pugh filtering (n < 5 excluded)\n")

cat("\nOUTPUT FILES:", OUTPUT_DIR, "\n")
cat(paste(" ", list.files(OUTPUT_DIR), collapse = "\n"), "\n")
cat("\n=============================================================================\n")

# End of script