#!/usr/bin/env Rscript
# ============================================================================
# FIGURE 6: ANCOM-BC2 DIFFERENTIAL ABUNDANCE HEATMAP (unadjusted + adjusted)
# Paper 1 — 16S rRNA Gut Microbiome | Saudi Arabian Liver Disease Cohort
# ----------------------------------------------------------------------------
# REVISION (Microbiology Spectrum, Reviewer 2 comments 6, 8, 24, 25):
#   * Differential abundance is now computed IN THIS SCRIPT with
#     ANCOMBC::ancombc2 (true ANCOM-BC2). The published version read QIIME2
#     `composition ancombc` output (= ANCOM-BC v1, Holm correction).
#   * BH-FDR within each contrast; pseudo-count sensitivity analysis;
#     structural-zero detection. All parameters printed + saved.
#   * Two models, one heatmap:
#       Left  block = unadjusted (Status)
#       Right block = adjusted for age, sex, BMI, PPI, metformin, statin
#     Rows = genera significant (q < 0.05 AND passing the pseudo-count
#     sensitivity 16S_liver_disease) in >= 1 unadjusted contrast.
#
# Generates:
#   Figure6_ANCOMBC2_Heatmap.pdf/.png/.svg
#   ANCOMBC2_R/ANCOMBC2_all_results.csv     (all tested genera, both models)
#   ANCOMBC2_R/ANCOMBC2_DA_taxa.csv         (significant genera -> Fig 7, S2, Table S3)
#   ANCOMBC2_R/ANCOMBC2_structural_zeros.csv
#   ANCOMBC2_R/ANCOMBC2_parameters.csv      (for Methods)
#
# JOURNAL SPECS — Gut Microbes (Taylor & Francis):
#   - Width: 174 mm | Height: scales with number of genera
#   - Resolution: 600 dpi (PNG); vector PDF; SVG editable
#   - Theme matches Figures 3, 4, 5, 7 for cross-figure consistency
# ============================================================================

rm(list = ls())
cat("\014")

suppressMessages({
  library(ANCOMBC)      # >= 2.8
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(svglite)
})
# Guard against masking (MASS, ANCOMBC deps, plyr, AnnotationDbi also export select/filter)
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename; mutate <- dplyr::mutate

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

# Paths RELATIVE to this script's folder (figure files are saved here)
#   script : 16S_liver_disease/Figure_6/
#   data   : Desktop/16S_reprocess/qiime2_analysis_v2/
proj_dir    <- file.path("..", "..", "16S_reprocess")
data_dir    <- file.path(proj_dir, "qiime2_analysis_v2")
counts_file <- file.path(data_dir, "Desq", "DeSeq_l6", "feature-table_L6.tsv")   # unrarefied genus counts
meta_file   <- file.path(proj_dir, "sample_metadata.tsv")
excl_file   <- file.path(proj_dir, "sample_exclusions.txt")   # optional: one SampleID per line
results_dir <- file.path(data_dir, "ANCOMBC2_R")               # shared with Fig 7, S2, Tables
order_file  <- file.path(results_dir, "Figure7_genus_order.csv")  # written by Figure 7
output_dir  <- getwd()

if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# ============================================================================
# SETTINGS (reported verbatim in Methods)
# ============================================================================

SEED      <- 2024
GROUPS    <- c("Healthy", "MASLD", "Cirrhosis", "HCC")
ANC <- list(prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
            struc_zero = TRUE,
            neg_lb = FALSE,          # HCC n = 17: package advises neg_lb = TRUE only for large groups
            alpha = 0.05, p_adj_method = "BH", pseudo_sens = TRUE)
ADJ_COVARS <- c("Age_z", "Gender", "BMI_z", "PPI_Use", "Metformin", "Statin_Use")
# Rifaximin (5 users) and lactulose (0 users in Healthy/MASLD) excluded: too sparse for
# taxon-level models. The 11-covariate model is reported in Supplementary Table S10.

# ============================================================================
# DEFENSIVE FILE-EXISTENCE CHECK
# ============================================================================

required_files <- c(counts = counts_file, metadata = meta_file)
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop("Missing input files:\n  ",
       paste(names(missing), missing, sep = ": ", collapse = "\n  "),
       "\n\nWorking directory: ", getwd(), "\nSave this script in 16S_liver_disease/Figure_6/ and click Source.")
}
cat("All input files located.\n\n")

# ============================================================================
# CONSISTENT COLOR PALETTE / THEME (same across all figures)
# ============================================================================

group_colors <- c("Healthy" = "#4DAF4A", "MASLD" = "#377EB8",
                  "Cirrhosis" = "#984EA3", "HCC" = "#E41A1C")

# ============================================================
# 1. LOAD DATA
# ============================================================

cat("Loading data...\n")

counts <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t", skip = 1,
                     check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
counts <- as.matrix(counts); storage.mode(counts) <- "numeric"

meta <- read.delim(meta_file, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
colnames(meta)[1] <- "SampleID"
meta <- meta[!grepl("^#", meta$SampleID), ]
meta$Status <- factor(gsub("^NAFLD$", "MASLD", meta$Status), levels = GROUPS)
rownames(meta) <- meta$SampleID

if (file.exists(excl_file)) {
  excl <- trimws(readLines(excl_file)); excl <- excl[nzchar(excl)]
  meta <- meta[!meta$SampleID %in% excl, ]
  cat("  Excluded", length(excl), "samples listed in sample_exclusions.txt\n")
}

common <- intersect(colnames(counts), meta$SampleID)
counts <- counts[, common]; meta <- meta[common, ]
counts <- counts[rowSums(counts) > 0, ]

# Covariates: same imputation rule as Figure 4 (median / mode); continuous z-scored
get_mode <- function(x) { x <- x[!is.na(x) & x != ""]; names(sort(table(x), decreasing = TRUE))[1] }
for (v in c("Age", "BMI")) {
  meta[[v]] <- suppressWarnings(as.numeric(meta[[v]]))
  meta[[v]][is.na(meta[[v]])] <- median(meta[[v]], na.rm = TRUE)
}
for (v in c("Gender", "PPI_Use", "Metformin", "Statin_Use")) {
  x <- as.character(meta[[v]]); x[is.na(x) | x == ""] <- get_mode(x); meta[[v]] <- factor(x)
}
meta$Age_z <- as.numeric(scale(meta$Age))
meta$BMI_z <- as.numeric(scale(meta$BMI))

cat(sprintf("  %d samples x %d genera | %s\n", ncol(counts), nrow(counts),
            paste(names(table(meta$Status)), table(meta$Status), sep = "=", collapse = ", ")))

# ============================================================
# 2. GENUS LABELS (same rule as previous Figure 6 / Supp S2)
# ============================================================

clean_genus <- function(taxon) {
  parts <- strsplit(taxon, ";")[[1]]
  for (offset in 1:7) {
    if (length(parts) >= offset) {
      name <- parts[length(parts) - offset + 1]
      name <- trimws(gsub("^[a-z]__", "", name))
      name <- gsub("_+", "_", name)
      if (!is.na(name) && nchar(name) >= 2 &&
          !name %in% c("unclassified", "uncultured", "metagenome", "NA", "Unknown") &&
          grepl("[a-zA-Z]", name)) {
        return(if (offset > 1) paste0("[", name, "]") else name)
      }
    }
  }
  paste0("Taxon_", substr(gsub("[^a-zA-Z0-9]", "", taxon), 1, 8))
}
clean_display <- function(name) gsub("_", " ", gsub("\\[|\\]", "", name))

# ============================================================
# 3. RUN ANCOM-BC2 (unadjusted and adjusted)
# ============================================================

run_ancombc2 <- function(formula_rhs, label) {
  cat(sprintf("Running ANCOM-BC2 [%s]: ~ %s ...\n", label, formula_rhs))
  set.seed(SEED)
  out <- ancombc2(data = counts, taxa_are_rows = TRUE, meta_data = meta,
                  fix_formula = formula_rhs, rand_formula = NULL,
                  p_adj_method = ANC$p_adj_method, pseudo_sens = ANC$pseudo_sens,
                  prv_cut = ANC$prv_cut, lib_cut = ANC$lib_cut, s0_perc = ANC$s0_perc,
                  group = "Status", struc_zero = ANC$struc_zero, neg_lb = ANC$neg_lb,
                  alpha = ANC$alpha, n_cl = 1, verbose = FALSE,
                  global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE)
  res <- out$res
  long <- bind_rows(lapply(GROUPS[-1], function(g) {
    col <- function(p) { k <- paste0(p, "Status", g); if (k %in% names(res)) res[[k]] else NA }
    data.frame(Model = label, Contrast = paste(g, "vs Healthy"), Taxon = res$taxon,
               LFC = col("lfc_"), SE = col("se_"), W = col("W_"), P = col("p_"), Q = col("q_"),
               Passed_SS = col("passed_ss_"), stringsAsFactors = FALSE)
  }))
  long$Significant <- !is.na(long$Q) & long$Q < ANC$alpha &
                      (is.na(long$Passed_SS) | long$Passed_SS)
  list(long = long, zero_ind = out$zero_ind, n_tested = nrow(res))
}

meta$Status <- relevel(meta$Status, ref = "Healthy")
fit_un  <- run_ancombc2("Status", "Unadjusted")
fit_adj <- run_ancombc2(paste(c("Status", ADJ_COVARS), collapse = " + "), "Adjusted")

all_res <- bind_rows(fit_un$long, fit_adj$long)
all_res$Genus <- vapply(all_res$Taxon, clean_genus, "")
all_res$Contrast <- factor(all_res$Contrast, levels = paste(GROUPS[-1], "vs Healthy"))

# ============================================================
# 4. SIGNIFICANT GENERA (unadjusted model defines the rows)
# ============================================================

sig_taxa <- unique(all_res$Taxon[all_res$Model == "Unadjusted" & all_res$Significant])
cat(sprintf("\n  Genera significant in >= 1 unadjusted contrast: %d\n", length(sig_taxa)))

# Diagnostic: version + counts with / without the pseudo-count sensitivity filter
cat(sprintf("  ANCOMBC %s | R %s | samples %d | taxa tested %d\n",
            as.character(packageVersion("ANCOMBC")), getRversion(), ncol(counts), fit_un$n_tested))
print(all_res %>% group_by(Model, Contrast) %>%
        dplyr::summarise(q_lt_0.05 = sum(Q < ANC$alpha, na.rm = TRUE),
                  q_and_passed_ss = sum(Significant), .groups = "drop") %>% as.data.frame())

da_tab <- all_res %>%
  dplyr::filter(Taxon %in% sig_taxa) %>%
  dplyr::select(Taxon, Genus, Model, Contrast, LFC, Q, Significant) %>%
  pivot_wider(names_from = c(Model, Contrast), values_from = c(LFC, Q, Significant),
              names_sep = " | ")

# ============================================================
# 5. ROW ORDER (import from Figure 7 if it matches; else Ward.D2 on LFC)
# ============================================================

lfc_un <- all_res %>% dplyr::filter(Model == "Unadjusted", Taxon %in% sig_taxa) %>%
  dplyr::select(Taxon, Contrast, LFC) %>% pivot_wider(names_from = Contrast, values_from = LFC)
lfc_mat <- as.matrix(lfc_un[, -1]); rownames(lfc_mat) <- lfc_un$Taxon

row_taxa <- NULL
if (file.exists(order_file)) {
  ord <- read.csv(order_file, stringsAsFactors = FALSE)
  if (setequal(ord$Taxon, sig_taxa)) {
    cat("  Row order imported from Figure 7 (same genera)\n")
    row_taxa <- ord$Taxon[order(ord$Order)]
  }
}
if (is.null(row_taxa)) {
  cat("  Row order: Ward.D2 clustering of unadjusted LFC\n")
  m <- lfc_mat; m[is.na(m)] <- 0
  row_taxa <- if (nrow(m) > 2) rownames(m)[hclust(dist(m), method = "ward.D2")$order] else rownames(m)
}

# ============================================================
# 6. PLOT DATA
# ============================================================

# * / ** / *** = q < 0.05 AND passed the pseudo-count sensitivity 16S_liver_disease
# \u2020 (dagger) = q < 0.05 but NOT robust to pseudo-count choice (failed passed_ss)
stars <- function(q, sig) ifelse(is.na(q), "",
                          ifelse(sig, ifelse(q < 0.001, "***", ifelse(q < 0.01, "**", "*")),
                                 ifelse(q < 0.05, "\u2020", "")))

plot_df <- all_res %>%
  dplyr::filter(Taxon %in% sig_taxa) %>%
  dplyr::select(Taxon, Model, Contrast, LFC, Q, Significant) %>%
  complete(Taxon = sig_taxa, Model = c("Unadjusted", "Adjusted"), Contrast) %>%   # untested cells -> grey
  dplyr::mutate(Significant = ifelse(is.na(Significant), FALSE, Significant),
         Sig    = stars(Q, Significant),
         Column = factor(sub(" vs Healthy", "", Contrast), levels = GROUPS[-1]),
         Panel  = factor(ifelse(Model == "Unadjusted", "Unadjusted",
                                "Adjusted (age, sex, BMI,\nPPI, metformin, statin)"),
                         levels = c("Unadjusted", "Adjusted (age, sex, BMI,\nPPI, metformin, statin)")),
         Genus_lab = factor(Taxon, levels = rev(row_taxa),
                            labels = clean_display(vapply(rev(row_taxa), clean_genus, ""))))

max_abs <- max(abs(plot_df$LFC), na.rm = TRUE)
lim <- ceiling(max_abs * 2) / 2

# ============================================================================
# 7. FIGURE 6: HEATMAP
# ============================================================================

cat("Creating Figure 6...\n")

figure6 <- ggplot(plot_df, aes(x = Column, y = Genus_lab, fill = LFC)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = Sig), size = 3.2, color = "black", fontface = "bold") +
  facet_grid(. ~ Panel) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-lim, lim), oob = scales::squish, na.value = "grey90",
                       name = "LFC (ln)\nvs Healthy",
                       breaks = seq(-floor(lim), floor(lim), by = 1)) +
  scale_x_discrete(position = "top") +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
        
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x.top   = element_text(size = 7.5, color = "black"),
    strip.text        = element_text(size = 8, face = "bold", lineheight = 0.9),
    strip.placement   = "outside",
    panel.spacing.x   = unit(4, "mm"),
    plot.caption.position = "plot",
    axis.text.y       = element_text(size = 7, face = "italic", color = "black"),
    axis.ticks        = element_blank(),
    panel.grid        = element_blank(),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    legend.position   = "right",
    legend.justification = c(0, 0.5),
    legend.title      = element_text(size = 8, face = "bold"),
    legend.text       = element_text(size = 7),
    legend.key.height = unit(1.5, "cm"),
    legend.key.width  = unit(0.35, "cm"),
    plot.caption      = element_text(size = 6, hjust = 0, color = "grey30"),
    plot.margin       = margin(10, 10, 10, 10)
  )

# ============================================================================
# 8. SAVE — Gut Microbes specifications
# ============================================================================

cat("Saving Figure 6...\n")
fig_height_mm <- max(100, length(sig_taxa) * 7 + 45)

ggsave("Figure6_ANCOMBC2_Heatmap.pdf", plot = figure6, width = 174, height = fig_height_mm, units = "mm")
ggsave("Figure6_ANCOMBC2_Heatmap.png", plot = figure6, width = 174, height = fig_height_mm, units = "mm",
       dpi = 600, bg = "white")
ggsave("Figure6_ANCOMBC2_Heatmap.svg", plot = figure6, width = 174, height = fig_height_mm, units = "mm",
       device = svglite)

# ============================================================
# 9. EXPORT RESULTS (consumed by Figure 7, Supp Fig S2, Supp Tables)
# ============================================================

write.csv(all_res, file.path(results_dir, "ANCOMBC2_all_results.csv"), row.names = FALSE)
write.csv(da_tab,  file.path(results_dir, "ANCOMBC2_DA_taxa.csv"), row.names = FALSE)

z <- fit_un$zero_ind
if (!is.null(z) && nrow(z)) {
  names(z) <- sub("structural_zero \\(Status = (.*)\\)", "\\1", names(z))
  z$Genus <- vapply(z$taxon, clean_genus, "")
  z <- z[rowSums(z[, intersect(GROUPS, names(z)), drop = FALSE]) > 0, ]
  write.csv(z, file.path(results_dir, "ANCOMBC2_structural_zeros.csv"), row.names = FALSE)
}

params <- data.frame(
  Parameter = c("Package", "Unadjusted formula", "Adjusted formula", "Reference group", "p_adj_method",
                "FDR scope", "prv_cut", "lib_cut", "s0_perc", "struc_zero", "neg_lb", "pseudo_sens",
                "alpha", "Significance rule", "Genera tested (unadjusted)", "Samples", "Seed"),
  Value = c(paste("ANCOMBC", packageVersion("ANCOMBC")), "Status",
            paste(c("Status", ADJ_COVARS), collapse = " + "), "Healthy", ANC$p_adj_method,
            "Across genera, within each contrast", ANC$prv_cut, ANC$lib_cut, ANC$s0_perc,
            ANC$struc_zero, ANC$neg_lb, ANC$pseudo_sens, ANC$alpha,
            "q < 0.05 AND passed pseudo-count sensitivity analysis",
            fit_un$n_tested, ncol(counts), SEED))
write.csv(params, file.path(results_dir, "ANCOMBC2_parameters.csv"), row.names = FALSE)

# ============================================================
# 10. SUMMARY
# ============================================================

cat("\n========================================\n")
cat("  FIGURE 6 — ANCOM-BC2 SUMMARY\n")
cat("========================================\n")
summ <- all_res %>% group_by(Model, Contrast) %>%
  dplyr::summarise(DA = sum(Significant), Enriched = sum(Significant & LFC > 0),
            Depleted = sum(Significant & LFC < 0), .groups = "drop")
print(as.data.frame(summ))
cat("\nSignificant genera (unadjusted) with adjusted estimates:\n")
print(all_res %>% dplyr::filter(Taxon %in% sig_taxa) %>%
        dplyr::select(Genus, Model, Contrast, LFC, Q) %>%
        dplyr::mutate(LFC = round(LFC, 2), Q = signif(Q, 2)) %>%
        pivot_wider(names_from = Model, values_from = c(LFC, Q)) %>%
        dplyr::arrange(Contrast, Genus) %>% as.data.frame())
cat("\nFiles:\n  ", output_dir, "/Figure6_ANCOMBC2_Heatmap.pdf/.png/.svg\n  ", results_dir, "/ANCOMBC2_*.csv\n", sep = "")
cat("\nRun order: Figure 6 -> Figure 7 -> Figure 6 again (to adopt Figure 7 row order) -> Supp Fig S2 -> Supp Tables\n")
