#!/usr/bin/env Rscript
# ============================================================================
# FIGURE 8: COMPOSITE DYSBIOSIS INDEX
# Paper 1 — 16S rRNA Gut Microbiome | Saudi Arabian Liver Disease Cohort
# ----------------------------------------------------------------------------
# Composite microbiome dysbiosis index derived from CLR-transformed genus
# abundance matrix → PCA → Euclidean distance from Healthy centroid → 0-100 scale.
#
# REVISION (Reviewer 2, point 12 — circularity): healthy controls are scored
# against a LEAVE-ONE-OUT centroid (each control excluded from the centroid it is
# measured from), computed here from the stored PC scores of the published DI
# (same PCA, same number of PCs, same rescaling). Clinical correlations (Panel C)
# are recomputed on the LOO DI; partial Spearman (adjusted for group), split-half
# and disease-only tests are written to CSV.
#
# Layout (4 panels, 2x2 grid):
#   A. Violin plot — dysbiosis distribution per disease group + KW + pairwise
#   B. Trajectory — mean ± 95% CI across disease progression
#   C. Forest plot — clinical correlations (Spearman + bootstrap CI)
#   D. PCA loadings — top 20 genera driving PC1 (positive vs negative direction)
#
# Generates: Figure8_Dysbiosis_Index.pdf/.png/.svg
#
# JOURNAL SPECS — Gut Microbes (Taylor & Francis):
#   - Width: 174 mm (full page)
#   - Height: 180 mm (proportional for 2x2 layout)
#   - Resolution: 600 dpi (PNG); vector PDF; SVG editable
#   - Theme matches Figures 3, 4, 5, 6, 7 for cross-figure consistency
# ============================================================================

rm(list = ls())
cat("\014")

suppressMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(boot)
  library(svglite)
  library(ppcor)
})
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

# ============================================================================
# PATHS
# ============================================================================

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
dysbiosis_dir <- file.path(data_dir, "dysbiosis_analysis")   # written by Dysbiosis_Index_PCA.R
output_dir    <- getwd()
N_BOOT        <- 1000
SEED          <- 2024

# ============================================================================
# DEFENSIVE FILE-EXISTENCE CHECK
# ============================================================================

required_files <- c(
  scores       = file.path(dysbiosis_dir, "Dysbiosis_Index_Scores.csv"),
  loadings     = file.path(dysbiosis_dir, "PCA_Loadings_Top_Genera.csv"),
  correlations = file.path(dysbiosis_dir, "Dysbiosis_Clinical_Correlations.csv")
)
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop("Missing input files:\n  ",
       paste(names(missing), missing, sep = ": ", collapse = "\n  "),
       "\n\nRun Dysbiosis_Index_PCA.R first (it writes these files).")
}
cat("All input files located.\n\n")

# ============================================================
# LOAD DATA
# ============================================================

cat("Loading dysbiosis data...\n")

# Dysbiosis scores (from PCA analysis)
pca_scores <- read.csv(file.path(dysbiosis_dir, "Dysbiosis_Index_Scores.csv"))
pca_scores$Status <- factor(gsub("^NAFLD$", "MASLD", pca_scores$Status), levels = c("Healthy", "MASLD", "Cirrhosis", "HCC"))

# PCA loadings
loadings_data <- read.csv(file.path(dysbiosis_dir, "PCA_Loadings_Top_Genera.csv"))

# Correlation results
correlation_results <- read.csv(file.path(dysbiosis_dir, "Dysbiosis_Clinical_Correlations.csv"))

cat("  Loaded", nrow(pca_scores), "samples\n")

# ---- Replace DI with the leave-one-out version (Reviewer 2, point 12) ----
# The scores file stores every PC score. Find the number of PCs (k) that reproduces
# the published DI exactly, then rescore each control against a centroid of the
# other controls only. Same PCA, same k, same 0-100 rescaling as the original DI.
pc_cols <- grep("^PC[0-9]+$", names(pca_scores), value = TRUE)
is_h    <- pca_scores$Status == "Healthy"
rescale01 <- function(x) 100 * (x - min(x)) / (max(x) - min(x))
di_from <- function(k) {
  S <- as.matrix(pca_scores[, pc_cols[1:k], drop = FALSE])
  rescale01(sqrt(rowSums(sweep(S, 2, colMeans(S[is_h, , drop = FALSE]))^2)))
}
err   <- sapply(seq_along(pc_cols), function(k) max(abs(di_from(k) - pca_scores$Dysbiosis_Index)))
K_PC  <- which.min(err)
if (err[K_PC] > 1e-6) stop("Could not reproduce the published DI from the stored PC scores (max diff ", signif(err[K_PC], 3), ")")
cat(sprintf("  Published DI reproduced exactly with %d PCs\n", K_PC))
S_k   <- as.matrix(pca_scores[, pc_cols[1:K_PC]])
raw   <- sqrt(rowSums(sweep(S_k, 2, colMeans(S_k[is_h, , drop = FALSE]))^2))
hS    <- S_k[is_h, , drop = FALSE]; nh <- nrow(hS)
raw[is_h] <- sapply(seq_len(nh), function(i) sqrt(sum((hS[i, ] - colMeans(hS[-i, , drop = FALSE]))^2)))
pca_scores$Dysbiosis_Index_original <- pca_scores$Dysbiosis_Index
pca_scores$Dysbiosis_Index          <- rescale01(raw)
cat("  Dysbiosis_Index now = leave-one-out DI (controls scored against centroid without themselves)\n")
print(pca_scores %>% group_by(Status) %>% summarise(n = n(), median_LOO = median(Dysbiosis_Index),
                                                    median_original = median(Dysbiosis_Index_original)))

# ---- Robustness checks that avoid circularity ----
# (a) Split-half: half the controls define the centroid, the other half are an
#     independent comparison group (200 random splits)
set.seed(SEED)
h_ids <- which(is_h)
split_p <- t(replicate(200, {
  ref <- sample(h_ids, floor(length(h_ids) / 2)); held <- setdiff(h_ids, ref)
  d   <- sqrt(rowSums(sweep(S_k, 2, colMeans(S_k[ref, , drop = FALSE]))^2))
  sapply(c("MASLD", "Cirrhosis", "HCC"), function(g)
    wilcox.test(d[pca_scores$Status == g], d[held], exact = FALSE)$p.value)
}))
# (b) Disease groups only (no healthy-centroid contrast)
dis <- droplevels(pca_scores[!is_h, ])
kw_dis <- kruskal.test(Dysbiosis_Index ~ Status, data = dis)
robust_tab <- data.frame(
  Test  = c(paste("Split-half:", colnames(split_p), "vs held-out Healthy"),
            "Kruskal-Wallis among MASLD, cirrhosis, HCC (LOO DI)"),
  Median_P        = c(apply(split_p, 2, median), kw_dis$p.value),
  Pct_P_lt_0.05   = c(100 * colMeans(split_p < 0.05), NA),
  Iterations      = c(rep(200, 3), NA), row.names = NULL)
write.csv(robust_tab, file.path(output_dir, "Figure8_DI_LOO_robustness.csv"), row.names = FALSE)
cat("\n  Non-circular checks:\n"); print(robust_tab)

# ---- Recompute clinical correlations on the LOO DI (pooled + partial) ----
grp_dummies <- model.matrix(~ Status, data = pca_scores)[, -1, drop = FALSE]
set.seed(SEED)
cor_rows <- lapply(as.character(correlation_results$Variable), function(v) {
  if (!v %in% names(pca_scores)) { warning("Variable ", v, " not in scores file; skipped"); return(NULL) }
  y <- suppressWarnings(as.numeric(pca_scores[[v]])); ok <- !is.na(y)
  x <- pca_scores$Dysbiosis_Index[ok]; yy <- y[ok]
  ct <- suppressWarnings(cor.test(x, yy, method = "spearman", exact = FALSE))
  bs <- replicate(N_BOOT, { i <- sample(length(x), replace = TRUE); suppressWarnings(cor(x[i], yy[i], method = "spearman")) })
  g_ok <- droplevels(pca_scores$Status[ok])
  pc <- if (nlevels(g_ok) < 2) list(estimate = NA_real_, p.value = NA_real_) else {
    Z <- model.matrix(~ g_ok)[, -1, drop = FALSE]
    suppressWarnings(pcor.test(x, yy, Z, method = "spearman"))
  }
  data.frame(Variable = v, N = sum(ok), Spearman_rho = unname(ct$estimate), P_value = ct$p.value,
             CI_lower = quantile(bs, 0.025, na.rm = TRUE), CI_upper = quantile(bs, 0.975, na.rm = TRUE),
             Partial_rho = pc$estimate, Partial_P = pc$p.value, row.names = NULL)
})
correlation_results <- bind_rows(cor_rows) %>%
  mutate(FDR = p.adjust(P_value, "BH"), Partial_FDR = p.adjust(Partial_P, "BH"),
         Significance = ifelse(FDR < 0.001, "***", ifelse(FDR < 0.01, "**", ifelse(FDR < 0.05, "*", ""))))
write.csv(correlation_results, file.path(output_dir, "Figure8_DI_LOO_clinical_correlations.csv"), row.names = FALSE)
cat("\n  LOO DI clinical correlations (pooled vs adjusted for group):\n")
print(correlation_results %>% transmute(Variable, rho = round(Spearman_rho, 3), q = signif(FDR, 2),
                                        partial_rho = round(Partial_rho, 3), partial_q = signif(Partial_FDR, 2)))

# ============================================================
# STYLING
# ============================================================

GROUP_COLORS <- c(
  "Healthy"   = "#4DAF4A",
  "MASLD"     = "#377EB8",
  "Cirrhosis" = "#984EA3",
  "HCC"       = "#E41A1C"
)

GROUP_ORDER <- c("Healthy", "MASLD", "Cirrhosis", "HCC")

pub_theme <- theme_minimal(base_size = 8) +
  theme(
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.grid         = element_blank(),
    panel.border       = element_blank(),
    axis.line          = element_line(color = "black", linewidth = 0.4),
    axis.ticks         = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length  = unit(0.08, "cm"),
    axis.text          = element_text(size = 7, color = "black"),
    axis.title         = element_text(size = 8.5, face = "bold"),
    plot.tag           = element_text(size = 10, face = "bold"),
    plot.margin        = margin(5, 5, 5, 5)
  )

# ============================================================
# PANEL A: VIOLIN PLOT
# ============================================================

# ---- Within-group correlations (DI vs BMI, age, GGT, ALT; 16 tests, BH) ----
wg <- do.call(rbind, lapply(levels(pca_scores$Status), function(g) do.call(rbind, lapply(c("BMI", "Age", "GGT", "ALT"), function(v) {
  x <- pca_scores[pca_scores$Status == g, ]
  ct <- suppressWarnings(cor.test(x$Dysbiosis_Index, as.numeric(x[[v]]), method = "spearman", exact = FALSE))
  data.frame(Group = g, Variable = v, rho = unname(ct$estimate), P = ct$p.value)
}))))
wg$Q <- p.adjust(wg$P, "BH")
write.csv(wg, file.path(output_dir, "Figure8_DI_within_group_correlations.csv"), row.names = FALSE)
cat(sprintf("\n  Within-group correlations: %d tests, min q = %.2f\n", nrow(wg), min(wg$Q)))

# ---- Cirrhosis medication check (LOO DI) ----
cir <- pca_scores[pca_scores$Status == "Cirrhosis", ]
is_yes <- function(v) as.character(v) %in% c("Yes", "yes", "1", "TRUE")
treated <- is_yes(cir$Lactulose_Use) | is_yes(cir$Rifaximin_Use) | is_yes(cir$Diuretic)
med_tab <- data.frame(
  Test = c("Lactulose, rifaximin or diuretic (any)", "Lactulose", "Rifaximin", "Diuretic", "Non-selective beta-blocker"),
  n_treated = c(sum(treated), sum(is_yes(cir$Lactulose_Use)), sum(is_yes(cir$Rifaximin_Use)),
                sum(is_yes(cir$Diuretic)), sum(is_yes(cir$BetaBlocker_Use))),
  P_Wilcoxon = sapply(list(treated, is_yes(cir$Lactulose_Use), is_yes(cir$Rifaximin_Use),
                           is_yes(cir$Diuretic), is_yes(cir$BetaBlocker_Use)), function(t)
    if (length(unique(t)) < 2) NA else wilcox.test(cir$Dysbiosis_Index[t], cir$Dysbiosis_Index[!t], exact = FALSE)$p.value))
write.csv(med_tab, file.path(output_dir, "Figure8_DI_cirrhosis_medications.csv"), row.names = FALSE)
cat("\n  Cirrhosis medication check (LOO DI):\n"); print(med_tab)

cat("Creating Panel A (violin)...\n")

# Sample sizes (computed for use in figure caption, not labels)
sample_sizes <- pca_scores %>%
  group_by(Status) %>%
  summarise(N = n(), .groups = "drop")

# X-axis labels: just disease group names (sample sizes go in figure caption)
x_labels <- setNames(as.character(sample_sizes$Status), sample_sizes$Status)

# Kruskal-Wallis
kw_test <- kruskal.test(Dysbiosis_Index ~ Status, data = pca_scores)

# Pairwise Wilcoxon rank-sum tests with Benjamini-Hochberg (BH-FDR) correction.
# Brackets use BH-adjusted q-values (NOT raw p-values), matching the manuscript
# Methods/Results. This is the single source of truth for the DI group statistics.
comparisons_list <- list(c("Healthy", "MASLD"), c("Healthy", "Cirrhosis"), c("Healthy", "HCC"))

pw <- pairwise.wilcox.test(
  pca_scores$Dysbiosis_Index,
  pca_scores$Status,
  p.adjust.method = "BH",
  exact = FALSE
)

# Print the full BH-adjusted q-value matrix to the console (all pairwise comparisons,
# including Healthy vs Cirrhosis). Report these q-values in the manuscript text.
cat(sprintf("\n  Kruskal-Wallis (LOO DI): chi2 = %.2f, p = %.2g\n", kw_test$statistic, kw_test$p.value))
cat("\n  PAIRWISE WILCOXON RANK-SUM, BH-ADJUSTED q-VALUES (LOO DI):\n")
print(round(pw$p.value, 4))

# Helper: retrieve the BH-adjusted q-value for a given pair (matrix is lower-triangular)
get_q <- function(a, b) {
  q <- suppressWarnings(tryCatch(pw$p.value[b, a], error = function(e) NA))
  if (is.na(q)) q <- suppressWarnings(tryCatch(pw$p.value[a, b], error = function(e) NA))
  q
}
pairwise_qvals <- sapply(comparisons_list, function(comp) get_q(comp[1], comp[2]))

get_stars <- function(q) ifelse(is.na(q), "ns",
                          ifelse(q < 0.001, "***",
                          ifelse(q < 0.01,  "**",
                          ifelse(q < 0.05,  "*", "ns"))))

y_max <- max(pca_scores$Dysbiosis_Index, na.rm = TRUE)

bracket_data <- data.frame(
  group1 = sapply(comparisons_list, `[`, 1),
  group2 = sapply(comparisons_list, `[`, 2),
  q_value = pairwise_qvals,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Stars = get_stars(q_value),
    x_start = match(group1, GROUP_ORDER),
    x_end = match(group2, GROUP_ORDER),
    bracket_height = y_max + c(0.06, 0.13, 0.20) * y_max,
    star_y = bracket_height + y_max * 0.012
  ) %>%
  filter(Stars != "ns")

panel_A <- ggplot(pca_scores, aes(x = Status, y = Dysbiosis_Index)) +
  geom_violin(aes(fill = Status), alpha = 0.6, trim = TRUE,
              linewidth = 0.4, color = "black") +
  geom_boxplot(width = 0.12, fill = "white", alpha = 0.9,
               outlier.shape = NA, linewidth = 0.3) +
  geom_jitter(aes(color = Status), width = 0.06, size = 0.8, alpha = 0.3) +
  # Brackets
  geom_segment(data = bracket_data,
               aes(x = x_start, xend = x_end, y = bracket_height, yend = bracket_height),
               linewidth = 0.4, color = "black", inherit.aes = FALSE) +
  geom_segment(data = bracket_data,
               aes(x = x_start, xend = x_start,
                   y = bracket_height - y_max * 0.015, yend = bracket_height + y_max * 0.015),
               linewidth = 0.4, color = "black", inherit.aes = FALSE) +
  geom_segment(data = bracket_data,
               aes(x = x_end, xend = x_end,
                   y = bracket_height - y_max * 0.015, yend = bracket_height + y_max * 0.015),
               linewidth = 0.4, color = "black", inherit.aes = FALSE) +
  geom_text(data = bracket_data,
            aes(x = (x_start + x_end) / 2, y = star_y, label = Stars),
            size = 3.2, fontface = "bold", vjust = 0.5, inherit.aes = FALSE) +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  scale_color_manual(values = GROUP_COLORS, guide = "none") +
  scale_x_discrete(labels = x_labels) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.04)), breaks = seq(0, 100, 25)) +
  labs(x = NULL, y = "Dysbiosis Index (leave-one-out)", tag = "A") +
  pub_theme +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(color = "black", size = 7)
  )

# ============================================================
# PANEL B: GROUP MEANS ± 95% CI (unconnected; groups are cross-sectional, Reviewer 2 point 23)
# ============================================================

cat("Creating Panel B (progression)...\n")

progression_data <- pca_scores %>%
  group_by(Status) %>%
  summarise(
    Mean = mean(Dysbiosis_Index),
    SE = sd(Dysbiosis_Index) / sqrt(n()),
    CI_lower = Mean - 1.96 * SE,
    CI_upper = Mean + 1.96 * SE,
    N = n(),
    .groups = "drop"
  ) %>%
  mutate(Stage = as.numeric(Status))

x_labels_prog <- as.character(progression_data$Status)

panel_B <- ggplot(progression_data, aes(x = Stage, y = Mean)) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.1, linewidth = 0.5, color = "gray40") +
  geom_point(aes(color = Status), size = 4) +
  scale_color_manual(values = GROUP_COLORS, guide = "none") +
  scale_x_continuous(breaks = 1:4, labels = x_labels_prog) +
  labs(x = NULL, y = "Dysbiosis Index (Mean ± 95% CI)", tag = "B") +
  pub_theme +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(color = "black", size = 7)
  )

# ============================================================
# PANEL C: FOREST PLOT (clinical correlations)
# ============================================================

cat("Creating Panel C (forest plot)...\n")

sig_cors <- correlation_results %>%
  filter(FDR < 0.05) %>%
  filter(!is.na(CI_lower), !is.na(CI_upper)) %>%
  arrange(-abs(Spearman_rho)) %>%
  mutate(
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

x_min <- min(sig_cors$CI_lower, na.rm = TRUE) - 0.05
x_max <- max(sig_cors$CI_upper, na.rm = TRUE) + 0.12
star_pos <- x_max - 0.02

panel_C <- ggplot(sig_cors, aes(x = Spearman_rho, y = Variable)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = Direction),
                 height = 0.25, linewidth = 0.6) +
  geom_point(aes(color = Direction), size = 2.5) +
  geom_text(aes(label = Significance), x = star_pos,
            size = 3, fontface = "bold", color = "black", hjust = 1, vjust = 0.5) +
  scale_color_manual(values = dir_colors) +
  scale_x_continuous(limits = c(x_min, x_max), breaks = seq(-0.4, 0.4, 0.2)) +
  labs(x = expression("Spearman's " * rho * " (95% Bootstrap CI)"), y = NULL, color = NULL, tag = "C") +
  pub_theme +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position    = "bottom",
    legend.text        = element_text(size = 7),
    legend.key.size    = unit(0.35, "cm")
  )
# ============================================================
# PANEL D: PCA LOADINGS
# ============================================================

cat("Creating Panel D (PCA loadings)...\n")

top_loadings <- head(loadings_data, 21)  # grab extra in case Unassigned is removed

# Remove Unassigned
top_loadings <- top_loadings[!grepl("^Unassigned$", top_loadings$Genus, ignore.case = TRUE), ]
top_loadings <- head(top_loadings, 20)

# Clean genus names
top_loadings$Genus_Clean <- top_loadings$Genus
top_loadings$Genus_Clean <- gsub("\\[|\\]", "", top_loadings$Genus_Clean)
top_loadings$Genus_Clean <- gsub("_", " ", top_loadings$Genus_Clean)

# Add family context to ambiguous names
top_loadings$Genus_Clean <- gsub("^UCG-002$", "Oscillospiraceae UCG-002", top_loadings$Genus_Clean)
top_loadings$Genus_Clean <- gsub("^UCG-005$", "Oscillospiraceae UCG-005", top_loadings$Genus_Clean)
top_loadings$Genus_Clean <- gsub("^UCG-010$", "Oscillospirales UCG-010", top_loadings$Genus_Clean)
top_loadings$Genus_Clean <- gsub("^NK4A214 group$", "Oscillospiraceae NK4A214 group", top_loadings$Genus_Clean)
top_loadings$Genus_Clean <- gsub("^UCG-003$", "Oscillospiraceae UCG-003", top_loadings$Genus_Clean)
top_loadings$Genus_Clean <- gsub("^RF39$", "Bacilli RF39", top_loadings$Genus_Clean)
top_loadings$Genus_Clean <- gsub("^Prevotella 9$", "Prevotella", top_loadings$Genus_Clean)
top_loadings$Genus_Clean <- gsub("^Prevotella 7$", "Prevotella", top_loadings$Genus_Clean)
# Remove "f " prefix from family-level assignments
top_loadings$Genus_Clean <- gsub("^f ", "", top_loadings$Genus_Clean)
# Clean trailing numbers from family-level (e.g., "Ruminococcaceae 1" -> "Ruminococcaceae")
top_loadings$Genus_Clean <- gsub(" 1$", "", top_loadings$Genus_Clean)
# Family-level assignments (applied after the "f " prefix is removed)
top_loadings$Genus_Clean <- gsub("^Ruminococcaceae$", "Ruminococcaceae (unclassified)", top_loadings$Genus_Clean)
top_loadings$Genus_Clean <- gsub("^Christensenellaceae$", "Christensenellaceae (unclassified)", top_loadings$Genus_Clean)

top_loadings$Direction <- ifelse(top_loadings$PC1 > 0, "Positive", "Negative")
top_loadings$Genus_Clean <- factor(top_loadings$Genus_Clean,
                                    levels = top_loadings$Genus_Clean[order(top_loadings$PC1)])

panel_D <- ggplot(top_loadings, aes(x = PC1, y = Genus_Clean, fill = Direction)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
  scale_fill_manual(values = c("Positive" = "#E41A1C", "Negative" = "#4DAF4A"),
                    labels = c("Positive" = "PC1 loading (+)", "Negative" = "PC1 loading (\u2212)"),
                    breaks = c("Negative", "Positive")) +
  labs(x = "PC1 Loading", y = NULL, fill = NULL, tag = "D") +
  pub_theme +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 7, face = "italic", color = "black"),
    legend.position    = "bottom",
    legend.text        = element_text(size = 7),
    legend.key.size    = unit(0.35, "cm")
  )

# ============================================================
# COMBINE 4 PANELS (2x2)
# ============================================================

cat("Assembling Figure 8...\n")

figure8 <- (panel_A | panel_B) / (panel_C | panel_D) +
  plot_layout(heights = c(0.5, 0.5))

# ============================================================================
# SAVE — Gut Microbes specifications
# Width: 174 mm | Height: 180 mm | 600 dpi PNG | vector PDF/SVG
# ============================================================================

cat("Saving Figure 8...\n")

ggsave("Figure8_Dysbiosis_Index.pdf",
       plot = figure8, width = 174, height = 180, units = "mm")

ggsave("Figure8_Dysbiosis_Index.png",
       plot = figure8, width = 174, height = 180, units = "mm",
       dpi = 600, bg = "white")

ggsave("Figure8_Dysbiosis_Index.svg",
       plot = figure8, width = 174, height = 180, units = "mm",
       device = svglite)

cat("\nFigure 8 saved (Gut Microbes format):\n")
cat("  - Figure8_Dysbiosis_Index.pdf  (vector, submission-ready)\n")
cat("  - Figure8_Dysbiosis_Index.png  (600 dpi raster)\n")
cat("  - Figure8_Dysbiosis_Index.svg  (editable vector)\n")
cat("Location:", getwd(), "\n")
