# ============================================================================
# FIGURE 7: SPEARMAN CORRELATION HEATMAP  (+ Supplementary Figure S4: partial)
# Paper 1 — 16S rRNA Gut Microbiome | Saudi Arabian Liver Disease Cohort
# ----------------------------------------------------------------------------
# REVISION (Microbiology Spectrum, Reviewer 2 comments 13-14):
#   * Genera are read from Figure 6 output (ANCOMBC2_DA_taxa.csv) instead of a
#     hard-coded list, so Figures 6 and 7 always describe the same genera.
#   * Pooled Spearman (Figure 7) AND partial Spearman adjusted for disease
#     group (Supplementary Figure S4) are computed in the same script.
#
# Cell color = Spearman rho; asterisks = BH-FDR (across all genus x parameter tests).
# Row order (Ward.D2 on pooled rho) is exported for Figure 6.
#
# Generates:
#   Figure7_Spearman_Heatmap.pdf/.png/.svg
#   SuppFigS4_Partial_Spearman_Heatmap.pdf/.png/.svg
#   Supplementary_Table_Spearman_Correlations.csv          (pooled  -> Table S5)
#   Supplementary_Table_Partial_Spearman_Correlations.csv  (partial -> Table S5)
#   ANCOMBC2_R/Figure7_genus_order.csv                     (row order -> Figure 6)
#
# JOURNAL SPECS — Gut Microbes (Taylor & Francis):
#   - Width: 174 mm | Height: scales with number of genera
#   - Resolution: 600 dpi (PNG); vector PDF; SVG editable
# ============================================================================

rm(list = ls())
cat("\014")

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(svglite)
  library(ppcor)
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

# Paths RELATIVE to this script's folder (figure files and tables are saved here)
#   script : 16S_liver_disease/Figure_7/
#   data   : Desktop/16S_reprocess/qiime2_analysis_v2/
proj_dir    <- file.path("..", "..", "16S_reprocess")
data_dir    <- file.path(proj_dir, "qiime2_analysis_v2")
counts_file <- file.path(data_dir, "Desq", "DeSeq_l6", "feature-table_L6.tsv")
meta_file   <- file.path(proj_dir, "sample_metadata.tsv")
excl_file   <- file.path(proj_dir, "sample_exclusions.txt")
results_dir <- file.path(data_dir, "ANCOMBC2_R")          # written by Figure 6
output_dir  <- getwd()
da_file <- file.path(results_dir, "ANCOMBC2_DA_taxa.csv")

# ============================================================================
# DEFENSIVE FILE-EXISTENCE CHECK
# ============================================================================

required_files <- c(counts = counts_file, metadata = meta_file, da_taxa = da_file)
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop("Missing input files:\n  ",
       paste(names(missing), missing, sep = ": ", collapse = "\n  "),
       "\n\nWorking directory: ", getwd(),
       "\nSave this script in 16S_liver_disease/Figure_7/ and click Source.",
       "\nRun Figure6_ANCOMBC2_Heatmap.R first (it writes ANCOMBC2_DA_taxa.csv).")
}
cat("All input files located.\n\n")

GROUPS <- c("Healthy", "MASLD", "Cirrhosis", "HCC")

# 11 clinical parameters: label = metadata column
clinical <- c(ALT = "ALT", AST = "AST", GGT = "GGT", ALP = "ALP",
              Albumin = "Albumin", Bilirubin = "Bili_T", Platelet = "Platelet",
              BMI = "BMI", HbA1c = "Hgb_A1C", TG = "TG", HDL = "HDL")

clean_genus <- function(taxon) {
  parts <- strsplit(taxon, ";")[[1]]
  for (offset in 1:7) {
    if (length(parts) >= offset) {
      name <- parts[length(parts) - offset + 1]
      name <- trimws(gsub("^[a-z]__", "", name)); name <- gsub("_+", "_", name)
      if (!is.na(name) && nchar(name) >= 2 &&
          !name %in% c("unclassified", "uncultured", "metagenome", "NA", "Unknown") &&
          grepl("[a-zA-Z]", name)) return(if (offset > 1) paste0("[", name, "]") else name)
    }
  }
  paste0("Taxon_", substr(gsub("[^a-zA-Z0-9]", "", taxon), 1, 8))
}
clean_display <- function(name) gsub("_", " ", gsub("\\[|\\]", "", name))

# ========================================================================
# 1. LOAD DATA
# ========================================================================

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
  excl <- trimws(readLines(excl_file)); meta <- meta[!meta$SampleID %in% excl[nzchar(excl)], ]
}
for (v in clinical) meta[[v]] <- suppressWarnings(as.numeric(meta[[v]]))

common <- intersect(colnames(counts), meta$SampleID)
counts <- counts[, common]; meta <- meta[common, ]

da <- read.csv(da_file, stringsAsFactors = FALSE, check.names = FALSE)
da_taxa <- da$Taxon
stopifnot(all(da_taxa %in% rownames(counts)))
cat(sprintf("  %d samples | %d DA genera from Figure 6\n", ncol(counts), length(da_taxa)))

# ========================================================================
# 2. RELATIVE ABUNDANCE (TSS on full table, then subset)
# ========================================================================

relab <- t(t(counts) / colSums(counts)) * 100
X <- t(relab[da_taxa, , drop = FALSE])            # samples x genera
labels <- vapply(da_taxa, clean_genus, "")

# ========================================================================
# 3. POOLED AND PARTIAL SPEARMAN
# ========================================================================

cat("Computing pooled and partial (group-adjusted) Spearman correlations...\n")
Z <- model.matrix(~ Status, data = meta)[, -1, drop = FALSE]

res <- expand.grid(Taxon = da_taxa, Parameter = names(clinical), stringsAsFactors = FALSE)
res$Genus <- labels[res$Taxon]
res[, c("Rho", "P", "Partial_Rho", "Partial_P", "N")] <- NA_real_
for (k in seq_len(nrow(res))) {
  x <- X[, res$Taxon[k]]; y <- meta[[clinical[[res$Parameter[k]]]]]
  ok <- complete.cases(x, y); if (sum(ok) < 10) next
  ct <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  pc <- pcor.test(x[ok], y[ok], Z[ok, , drop = FALSE], method = "spearman")
  res[k, c("Rho", "P", "Partial_Rho", "Partial_P", "N")] <- c(ct$estimate, ct$p.value, pc$estimate, pc$p.value, sum(ok))
}
res$Q         <- p.adjust(res$P, "BH")
res$Partial_Q <- p.adjust(res$Partial_P, "BH")
star <- function(q) ifelse(is.na(q), "", ifelse(q < 0.001, "***", ifelse(q < 0.01, "**", ifelse(q < 0.05, "*", ""))))
res$Sig <- star(res$Q); res$Partial_Sig <- star(res$Partial_Q)

cat(sprintf("  Pooled:  %d / %d significant (FDR < 0.05)\n", sum(res$Q < 0.05, na.rm = TRUE), sum(!is.na(res$Q))))
cat(sprintf("  Partial: %d / %d significant (FDR < 0.05); min q = %.3f\n",
            sum(res$Partial_Q < 0.05, na.rm = TRUE), sum(!is.na(res$Partial_Q)), min(res$Partial_Q, na.rm = TRUE)))

# ========================================================================
# 4. ROW ORDER (Ward.D2 on pooled rho) -> exported for Figure 6
# ========================================================================

rho_mat <- res %>% dplyr::select(Taxon, Parameter, Rho) %>% pivot_wider(names_from = Parameter, values_from = Rho)
m <- as.matrix(rho_mat[, -1]); rownames(m) <- rho_mat$Taxon; m[is.na(m)] <- 0
row_taxa <- if (nrow(m) > 2) rownames(m)[hclust(dist(m), method = "ward.D2")$order] else rownames(m)
write.csv(data.frame(Taxon = row_taxa, Genus = labels[row_taxa], Order = seq_along(row_taxa)),
          file.path(results_dir, "Figure7_genus_order.csv"), row.names = FALSE)

param_order <- c("ALT", "AST", "GGT", "ALP", "Albumin", "Bilirubin", "Platelet", "BMI", "HbA1c", "TG", "HDL")
res$Genus_lab <- factor(res$Taxon, levels = rev(row_taxa), labels = clean_display(labels[rev(row_taxa)]))
res$Parameter <- factor(res$Parameter, levels = param_order)

# ============================================================================
# 5. PLOTS (shared style)
# ============================================================================

heat <- function(df, fill_col, lab_col, legend_title, caption = NULL) {
  ggplot(df, aes(x = Parameter, y = Genus_lab, fill = .data[[fill_col]])) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = .data[[lab_col]]), size = 3.2, color = "black", fontface = "bold") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-0.4, 0.4), oob = scales::squish, breaks = seq(-0.4, 0.4, 0.2),
                         na.value = "grey90", name = legend_title) +
    labs(x = NULL, y = NULL, caption = caption) +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 7, face = "italic", color = "black"),
          axis.ticks = element_blank(), panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          legend.position = "right", legend.justification = c(0, 0.5),
          legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7),
          legend.key.height = unit(1.5, "cm"), legend.key.width = unit(0.35, "cm"),
          plot.margin = margin(10, 10, 10, 10))
}

# No in-figure captions: test details go in the manuscript figure legends
n_tests <- sum(!is.na(res$Q))
figure7 <- heat(res, "Rho",         "Sig",         expression(Spearman~rho))
suppS4  <- heat(res, "Partial_Rho", "Partial_Sig", expression(Partial~Spearman~rho))

# ============================================================================
# 6. SAVE — Gut Microbes specifications
# ============================================================================

h <- max(110, length(da_taxa) * 7 + 50)
for (ext in c("pdf", "png", "svg")) {
  args <- list(width = 174, height = h, units = "mm")
  if (ext == "png") args <- c(args, dpi = 600, bg = "white")
  if (ext == "svg") args$device <- svglite
  do.call(ggsave, c(list(filename = paste0("Figure7_Spearman_Heatmap.", ext), plot = figure7), args))
  do.call(ggsave, c(list(filename = paste0("SuppFigS4_Partial_Spearman_Heatmap.", ext), plot = suppS4), args))
}

out <- res %>% dplyr::arrange(match(Taxon, row_taxa), Parameter) %>%
  transmute(Genus = clean_display(Genus), Clinical_Parameter = as.character(Parameter), N,
            Spearman_Rho = round(Rho, 4), P_value = signif(P, 4), FDR_Q = signif(Q, 4), Significance = Sig,
            Partial_Spearman_Rho = round(Partial_Rho, 4), Partial_P_value = signif(Partial_P, 4),
            Partial_FDR_Q = signif(Partial_Q, 4), Partial_Significance = Partial_Sig)
write.csv(out[, c("Genus", "Clinical_Parameter", "Spearman_Rho", "P_value", "FDR_Q", "Significance")],
          "Supplementary_Table_Spearman_Correlations.csv", row.names = FALSE)
write.csv(out[, c("Genus", "Clinical_Parameter", "Partial_Spearman_Rho", "Partial_P_value", "Partial_FDR_Q", "Partial_Significance")],
          "Supplementary_Table_Partial_Spearman_Correlations.csv", row.names = FALSE)
write.csv(out, file.path(results_dir, "Spearman_pooled_and_partial.csv"), row.names = FALSE)

# ========================================================================
# 7. SUMMARY
# ========================================================================

cat("\n", strrep("=", 70), "\nFIGURE 7 — SPEARMAN SUMMARY\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("  Genera: %d | parameters: %d | tests: %d\n", length(da_taxa), length(clinical), n_tests))
cat(sprintf("  Pooled significant (q<0.05): %d | partial significant: %d\n",
            sum(res$Q < 0.05, na.rm = TRUE), sum(res$Partial_Q < 0.05, na.rm = TRUE)))
cat("  Strongest pooled correlations:\n")
top <- res %>% dplyr::filter(Q < 0.05) %>% dplyr::arrange(desc(abs(Rho))) %>% head(10)
for (i in seq_len(nrow(top))) cat(sprintf("    %-30s x %-9s rho = %+.2f  q = %.3g | partial rho = %+.2f\n",
                                          top$Genus[i], top$Parameter[i], top$Rho[i], top$Q[i], top$Partial_Rho[i]))
cat("\nFiles saved in:", output_dir, "\n")
