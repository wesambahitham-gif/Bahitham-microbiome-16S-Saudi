#!/usr/bin/env Rscript
# ============================================================================
# SUPPLEMENTARY FIGURE S2: CROSS-METHOD CONCORDANCE (ANCOM-BC2 vs DESeq2)
# Paper 1 — 16S rRNA Gut Microbiome | Saudi Arabian Liver Disease Cohort
# ----------------------------------------------------------------------------
# REVISION (Microbiology Spectrum, Reviewer 2 comments 7, 8):
#   * DESeq2 is a complementary method on the same data, not "independent
#     validation": wording changed to "concordance" throughout.
#   * ONE DESeq2 model containing all four groups (Healthy reference), matching
#     the ANCOM-BC2 design of Figure 6 (previously three pairwise subsets).
#   * ANCOM-BC2 results are read from Figure 6 output (true ANCOM-BC2).
#   * 90% "high-confidence" reference line removed; colours = main figures.
#
# Layout (6 panels):
#   A-C. Venn diagrams per contrast (ANCOM-BC2 vs DESeq2, q < 0.05)
#   D.   Concordance rate (% of ANCOM-BC2 genera also significant in DESeq2,
#        same direction)
#   E.   Summary table
#   F.   log2FC heatmap of concordant genera (both methods)
#
# Generates: Supplementary_Figure_S2_Concordance.pdf/.png/.svg
#            Supplementary_Tables/Table_S4_Concordance_Details.xlsx  (-> Table S4)
#
# JOURNAL SPECS — Gut Microbes (Taylor & Francis): 174 mm wide; 600 dpi PNG; PDF/SVG
# ============================================================================

rm(list = ls())
cat("\014")

suppressPackageStartupMessages({
  library(DESeq2); library(openxlsx); library(ggplot2); library(dplyr); library(tidyr)
  library(ggVennDiagram); library(gridExtra); library(cowplot); library(scales); library(svglite)
})
# Guard against masking (MASS, ANCOMBC deps, plyr, AnnotationDbi also export select/filter)
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename; mutate <- dplyr::mutate

# =============================================================================
# CONFIGURATION
# =============================================================================

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

# Paths RELATIVE to this script's folder (figure + Table S4 are saved here)
#   script : any sub-folder of Desktop (e.g. Desktop/Reviewers/Figure_7/Supp_Fig_S2/)
#   data   : Desktop/16S_reprocess/qiime2_analysis_v2/
# Find Desktop/16S_reprocess by walking up from this script's folder,
# so the script works from any sub-folder of Desktop (e.g. Reviewers/Figure_7/Supp_Fig_S2/)
find_proj_dir <- function(name = "16S_reprocess") {
  d <- normalizePath(getwd())
  repeat {
    cand <- file.path(d, name)
    if (dir.exists(cand)) return(cand)
    up <- dirname(d); if (up == d) break; d <- up
  }
  stop("Could not find a '", name, "' folder above ", getwd())
}
proj_dir    <- find_proj_dir()
data_dir    <- file.path(proj_dir, "qiime2_analysis_v2")
counts_file <- file.path(data_dir, "Desq", "DeSeq_l6", "feature-table_L6.tsv")
meta_file   <- file.path(proj_dir, "sample_metadata.tsv")
excl_file   <- file.path(proj_dir, "sample_exclusions.txt")
results_dir <- file.path(data_dir, "ANCOMBC2_R")          # written by Figure 6
output_dir  <- getwd()

anc_file <- file.path(results_dir, "ANCOMBC2_all_results.csv")
dir.create(file.path(output_dir, "Supplementary_Tables"), showWarnings = FALSE, recursive = TRUE)

FDR_CUTOFF <- 0.05
MIN_COUNT <- 3          # DESeq2 filter: >= 3 reads in >= 5% of samples (unchanged)
MIN_SAMPLE_PCT <- 0.05
GROUPS <- c("Healthy", "MASLD", "Cirrhosis", "HCC")
DISEASE <- GROUPS[-1]

required_files <- c(counts = counts_file, metadata = meta_file, ancombc2 = anc_file)
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) stop("Missing input files:\n  ", paste(names(missing), missing, sep = ": ", collapse = "\n  "),
                              "\n\nWorking directory: ", getwd(),
                              "\nOpen the script from its folder and click Source.",
                              "\nRun Figure6_ANCOMBC2_Heatmap.R first.")
cat("All input files located.\n\n")

# Same palette as all main figures
group_colors <- c("Healthy" = "#4DAF4A", "MASLD" = "#377EB8", "Cirrhosis" = "#984EA3", "HCC" = "#E41A1C")

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
clean_display <- function(x) gsub("_", " ", gsub("\\[|\\]", "", x))

# =============================================================================
# 1. LOAD DATA
# =============================================================================

counts <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t", skip = 1,
                     check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
counts <- as.matrix(counts); storage.mode(counts) <- "numeric"
meta <- read.delim(meta_file, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
colnames(meta)[1] <- "SampleID"; meta <- meta[!grepl("^#", meta$SampleID), ]
meta$Status <- factor(gsub("^NAFLD$", "MASLD", meta$Status), levels = GROUPS)
rownames(meta) <- meta$SampleID
if (file.exists(excl_file)) { excl <- trimws(readLines(excl_file)); meta <- meta[!meta$SampleID %in% excl[nzchar(excl)], ] }
common <- intersect(colnames(counts), meta$SampleID)
counts <- counts[, common]; meta <- meta[common, ]
n_text <- paste(sprintf("n=%d %s", as.integer(table(meta$Status)), GROUPS), collapse = ", ")

anc <- read.csv(anc_file, stringsAsFactors = FALSE) %>% dplyr::filter(Model == "Unadjusted")

# =============================================================================
# 2. DESeq2 — one model, all groups, Healthy reference
# =============================================================================

cat("Running DESeq2 (single model, ~ Status)...\n")
keep <- rowSums(counts >= MIN_COUNT) >= ceiling(MIN_SAMPLE_PCT * ncol(counts))
dds <- DESeqDataSetFromMatrix(round(counts[keep, ]), meta, ~ Status)
dds$Status <- relevel(dds$Status, ref = "Healthy")
dds <- estimateSizeFactors(dds, type = "poscounts")    # geometric mean of positive counts (as before)
dds <- DESeq(dds, fitType = "local", quiet = TRUE)
des <- bind_rows(lapply(DISEASE, function(g) {
  r <- results(dds, contrast = c("Status", g, "Healthy"), cooksCutoff = FALSE, independentFiltering = FALSE)
  data.frame(Contrast = paste(g, "vs Healthy"), Taxon = rownames(r), DES_LFC = r$log2FoldChange,
             DES_SE = r$lfcSE, DES_P = r$pvalue, DES_Q = r$padj, stringsAsFactors = FALSE)
}))
cat(sprintf("  %d genera tested after filter\n", sum(keep)))

# =============================================================================
# 3. CONCORDANCE
# =============================================================================

comb <- anc %>% dplyr::select(Contrast, Taxon, ANC_LFC = LFC, ANC_Q = Q, ANC_Sig = Significant) %>%
  full_join(des, by = c("Contrast", "Taxon")) %>%
  dplyr::mutate(Genus = vapply(Taxon, clean_genus, ""),
         ANC_Sig = ANC_Sig %in% TRUE,
         DES_Sig = !is.na(DES_Q) & DES_Q < FDR_CUTOFF,
         Concordant = ANC_Sig & DES_Sig & sign(ANC_LFC) == sign(DES_LFC),
         Discordant_direction = ANC_Sig & DES_Sig & sign(ANC_LFC) != sign(DES_LFC))

summary_df <- comb %>% dplyr::mutate(Comparison = sub(" vs Healthy", "", Contrast)) %>%
  group_by(Comparison) %>%
  dplyr::summarise(`ANCOM-BC2` = sum(ANC_Sig), DESeq2 = sum(DES_Sig), Concordant = sum(Concordant),
            .groups = "drop") %>%
  dplyr::mutate(`Concordance (%)` = round(100 * Concordant / pmax(1, `ANCOM-BC2`), 1),
         Comparison = factor(Comparison, levels = DISEASE)) %>% dplyr::arrange(Comparison)
print(as.data.frame(summary_df))
if (any(comb$Discordant_direction)) cat("  Discordant direction:", paste(comb$Genus[comb$Discordant_direction], collapse = ", "), "\n")

conc_taxa <- unique(comb$Taxon[comb$Concordant])
cat(sprintf("  Concordant genera (any contrast): %d -> %s\n", length(conc_taxa),
            paste(vapply(conc_taxa, clean_genus, ""), collapse = ", ")))

# =============================================================================
# 4. TABLE S4 (Excel)
# =============================================================================

stars <- function(q) ifelse(is.na(q), "", ifelse(q < 0.001, "***", ifelse(q < 0.01, "**", ifelse(q < 0.05, "*", ""))))
anc_sig_taxa <- unique(comb$Taxon[comb$ANC_Sig])
tab_s4 <- comb %>% dplyr::filter(Taxon %in% anc_sig_taxa) %>%
  dplyr::mutate(Comparison = sub(" vs Healthy", "", Contrast)) %>%
  dplyr::select(Taxon, Genus, Comparison, ANC_LFC, ANC_Q, ANC_Sig, DES_LFC, DES_Q, DES_Sig, Concordant) %>%
  pivot_wider(names_from = Comparison, values_from = c(ANC_LFC, ANC_Q, ANC_Sig, DES_LFC, DES_Q, DES_Sig, Concordant)) %>%
  dplyr::mutate(Contrasts_concordant = rowSums(across(starts_with("Concordant_")), na.rm = TRUE)) %>%
  dplyr::arrange(desc(Contrasts_concordant), Genus)

wb <- createWorkbook()
hs <- createStyle(fontSize = 11, fontColour = "#FFFFFF", halign = "center", fgFill = "#2C3E50", textDecoration = "bold")
for (sh in list(list("Summary", summary_df), list("Concordance_Details", tab_s4), list("All_genera", comb),
                list("Legend", data.frame(Item = c("ANC_*", "DES_*", "Concordant", "Q"),
                                          Meaning = c("ANCOM-BC2 (Figure 6, unadjusted)", "DESeq2 (single model, ~Status)",
                                                      "q < 0.05 in both methods, same direction", "BH-FDR within contrast"))))) {
  addWorksheet(wb, sh[[1]]); writeData(wb, sh[[1]], sh[[2]])
  addStyle(wb, sh[[1]], hs, rows = 1, cols = seq_len(ncol(sh[[2]])), gridExpand = TRUE)
  setColWidths(wb, sh[[1]], cols = seq_len(ncol(sh[[2]])), widths = "auto")
}
saveWorkbook(wb, file.path(output_dir, "Supplementary_Tables", "Table_S4_Concordance_Details.xlsx"), overwrite = TRUE)
write.csv(comb, file.path(results_dir, "DESeq2_vs_ANCOMBC2_unadjusted.csv"), row.names = FALSE)

# Agreement beyond the q < 0.05 overlap (Reviewer 2, point 8):
#   (a) Spearman correlation of effect sizes across all genera tested by both methods
#   (b) direction agreement and DESeq2 nominal p for every ANCOM-BC2 hit
both_tested <- comb %>% dplyr::filter(!is.na(ANC_LFC), !is.na(DES_LFC))
agree_all <- both_tested %>% dplyr::group_by(Contrast) %>%
  dplyr::summarise(Genera_tested_both = dplyr::n(),
                   Spearman_rho_LFC   = cor(ANC_LFC, DES_LFC, method = "spearman"),
                   Spearman_P         = cor.test(ANC_LFC, DES_LFC, method = "spearman", exact = FALSE)$p.value,
                   .groups = "drop")
agree_hits <- both_tested %>% dplyr::filter(ANC_Sig) %>%
  dplyr::transmute(Contrast, Genus, ANC_LFC, ANC_Q, DES_LFC, DES_P, DES_Q,
                   Same_direction = sign(ANC_LFC) == sign(DES_LFC),
                   DESeq2_nominal_P_lt_0.05 = DES_P < 0.05,
                   DESeq2_q_lt_0.05 = DES_Q < 0.05)
write.csv(agree_all,  file.path(results_dir, "Concordance_LFC_correlation.csv"), row.names = FALSE)
write.csv(agree_hits, file.path(results_dir, "Concordance_ANCOMBC2_hits_in_DESeq2.csv"), row.names = FALSE)
cat("\nEffect-size agreement (all genera tested by both methods):\n"); print(as.data.frame(agree_all))
cat("\nANCOM-BC2 hits in DESeq2:\n"); print(as.data.frame(agree_hits))

# =============================================================================
# 5. PANELS
# =============================================================================

venn_theme <- theme_void() + theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
                                   plot.subtitle = element_text(size = 7.5, hjust = 0.5, color = "gray30"),
                                   legend.position = "none", plot.margin = margin(5, 5, 5, 5))
venn_plots <- lapply(DISEASE, function(g) {
  sub <- comb %>% dplyr::filter(Contrast == paste(g, "vs Healthy"))
  a <- sub$Taxon[sub$ANC_Sig]; d <- sub$Taxon[sub$DES_Sig]; n_c <- sum(sub$Concordant)
  ggVennDiagram(list("ANCOM-BC2" = a, "DESeq2" = d), label = "none", label_alpha = 0, edge_size = 1, set_size = 0) +
    scale_fill_gradient(low = "white", high = group_colors[[g]]) +
    annotate("text", x = 0.5, y = 6,    label = "DESeq2",    size = 3, fontface = "bold") +
    annotate("text", x = 0.5, y = 4.8,  label = length(d),   size = 3) +
    annotate("text", x = 0.5, y = -0.42, label = "ANCOM-BC2", size = 3, fontface = "bold") +
    annotate("text", x = 0.5, y = -1.5, label = length(a),   size = 3) +
    annotate("text", x = 0.5, y = 2.5,  label = n_c,         size = 3.6, fontface = "bold") +
    labs(title = sprintf("%s vs Healthy", g),
         subtitle = sprintf("Concordant: %d/%d (%.0f%%)", n_c, length(a), 100 * n_c / max(1, length(a)))) +
    venn_theme
})

conc_bar <- ggplot(summary_df, aes(x = Comparison, y = `Concordance (%)`, fill = Comparison)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%", `Concordance (%)`)), vjust = -0.5, size = 2.8, fontface = "bold") +
  scale_fill_manual(values = group_colors[DISEASE]) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "ANCOM-BC2 genera also\nsignificant in DESeq2 (%)", subtitle = n_text) +
  theme_minimal(base_size = 8) +
  theme(axis.text = element_text(size = 7, color = "black"), axis.title.y = element_text(size = 8, face = "bold"),
        axis.line = element_line(color = "black", linewidth = 0.4), axis.ticks = element_line(color = "black", linewidth = 0.3),
        panel.grid = element_blank(), legend.position = "none",
        plot.subtitle = element_text(size = 6.5, hjust = 0.5, color = "gray40"))

tab_show <- as.data.frame(summary_df); names(tab_show)[names(tab_show) == "Concordance (%)"] <- "Concord. (%)"
table_plot <- ggdraw() + draw_grob(tableGrob(tab_show, rows = NULL,
  theme = ttheme_minimal(core = list(fg_params = list(fontsize = 6.5), bg_params = list(fill = c("grey95", "white"))),
                         colhead = list(fg_params = list(fontsize = 6.5, fontface = "bold"), bg_params = list(fill = "grey80")))))

# Panel F: heatmap of concordant genera
hm <- comb %>% dplyr::filter(Taxon %in% conc_taxa) %>%
  dplyr::select(Taxon, Genus, Contrast, ANC_LFC, ANC_Q, ANC_Sig, DES_LFC, DES_Q) %>%
  pivot_longer(c(ANC_LFC, DES_LFC), names_to = "Method", values_to = "LFC") %>%
  dplyr::mutate(Method = ifelse(Method == "ANC_LFC", "ANCOM-BC2", "DESeq2"),
         Sig = ifelse(Method == "ANCOM-BC2", ifelse(ANC_Sig, stars(ANC_Q), ""), stars(DES_Q)),
         Column = factor(sub(" vs Healthy", "", Contrast), levels = DISEASE),
         Genus_lab = clean_display(Genus))
p_heat <- ggplot(hm, aes(x = Column, y = Genus_lab, fill = LFC)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = Sig), size = 3, fontface = "bold") +
  facet_grid(. ~ Method) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-4, 4),
                       oob = squish, na.value = "grey90",
                       name = "LFC\n(ANCOM-BC2: ln;\nDESeq2: log2)") +
  labs(x = NULL, y = NULL) +            # no in-figure caption: details go in the figure legend
  theme_minimal(base_size = 8) +
  theme(axis.text.y      = element_text(face = "italic", size = 7, color = "black"),
        axis.text.x      = element_text(size = 7, color = "black"),
        axis.ticks       = element_blank(),
        strip.text       = element_text(size = 8, face = "bold"),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.grid       = element_blank(),
        legend.title     = element_text(size = 8, face = "bold"),
        legend.text      = element_text(size = 7),
        plot.margin      = margin(5, 5, 5, 15))

# =============================================================================
# 6. ASSEMBLE + EXPORT
# =============================================================================

top_row    <- plot_grid(plotlist = venn_plots, ncol = 3, labels = c("A", "B", "C"), label_size = 10)
middle_row <- plot_grid(conc_bar, table_plot, ncol = 2, labels = c("D", "E"), label_size = 10, rel_widths = c(1, 1.4))
bottom_row <- plot_grid(p_heat, labels = "F", label_size = 10)
fig_s2 <- plot_grid(top_row, middle_row, bottom_row, ncol = 1,
                    rel_heights = c(0.8, 0.7, max(0.6, 0.12 * length(conc_taxa) + 0.3)))

h <- 150 + 12 * length(conc_taxa)
ggsave(file.path(output_dir, "Supplementary_Figure_S2_Concordance.pdf"), fig_s2, width = 174, height = h, units = "mm")
ggsave(file.path(output_dir, "Supplementary_Figure_S2_Concordance.png"), fig_s2, width = 174, height = h, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(output_dir, "Supplementary_Figure_S2_Concordance.svg"), fig_s2, width = 174, height = h, units = "mm", device = svglite)

cat("\nSupplementary Figure S2 saved:", output_dir, "\n")
cat("  Table_S4_Concordance_Details.xlsx (Supplementary_Tables)\n")
