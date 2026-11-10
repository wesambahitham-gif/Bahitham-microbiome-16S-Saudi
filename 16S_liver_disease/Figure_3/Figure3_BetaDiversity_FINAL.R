#!/usr/bin/env Rscript
# ============================================================================
# FIGURE 3: BETA DIVERSITY — PCoA (centroid separation) + BETADISPER (variance)
# Paper 1 — 16S rRNA Gut Microbiome | Saudi Arabian Liver Disease Cohort
# ----------------------------------------------------------------------------
# Layout (4 rows x 2 cols = 8 panels, A–H):
#   Row 1:  A. Bray-Curtis PCoA          | B. Bray-Curtis Betadisper
#   Row 2:  C. Jaccard PCoA              | D. Jaccard Betadisper
#   Row 3:  E. Unweighted UniFrac PCoA   | F. Unweighted UniFrac Betadisper
#   Row 4:  G. Weighted UniFrac PCoA     | H. Weighted UniFrac Betadisper
#
# Logic: each ROW = one distance metric; left column = centroid result (PCoA
# + PERMANOVA stats), right column = dispersion result (betadisper violin +
# Tukey HSD significance bars). This visually enacts the methodological
# distinction between PERMANOVA (centroid) and PERMDISP (variance homogeneity).
#
# JOURNAL SPECS — Gut Microbes (Taylor & Francis):
#   - Width: 174 mm (double-column / full page width)
#   - Height: 215 mm (proportional for 4-row layout)
#   - Resolution: 600 dpi (PNG); vector PDF
#   - Min line width: 0.4 pt; embedded fonts
# ============================================================================

rm(list = ls())
cat("\014")

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(vegan)
  library(ape)
  library(rstatix)
  library(patchwork)
  library(cowplot)
  library(svglite)
})

# ============================================================================
# PATHS — relative to this script's folder (outputs are saved here)
#   script : 16S_liver_disease/Figure_3/
#   data   : Desktop/16S_reprocess/qiime2_analysis_v2/
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

data_dir <- file.path("..", "..", "16S_reprocess", "qiime2_analysis_v2")
beta_dir <- file.path(data_dir, "Beta_diversity")

required_files <- c(
  bray     = file.path(beta_dir, "Bray-curtis",        "Bray-Curtis_raw_data.tsv"),
  jaccard  = file.path(beta_dir, "Jaccard",            "Jaccardraw_data.tsv"),
  unweight = file.path(beta_dir, "UnWeighted_UniFrac", "Unweighted.tsv"),
  weighted = file.path(beta_dir, "Weighted_UniFrac",   "Weighted_rowdata.tsv"))
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) stop("Missing input files:\n  ",
  paste(names(missing), missing, sep = ": ", collapse = "\n  "),
  "\n\nWorking directory: ", getwd(),
  "\nSave this script in 16S_liver_disease/Figure_3/ and click Source.")
cat("All input files located.\n\n")

# ============================================================================
# CONSISTENT COLOR PALETTE (must match Figures 2, 4–8)
# ============================================================================

group_colors <- c(
  "Healthy"   = "#4DAF4A",
  "MASLD"     = "#377EB8",
  "Cirrhosis" = "#984EA3",
  "HCC"       = "#E41A1C"
)

# ============================================================================
# UNIFIED PUBLICATION THEME (Gut Microbes specs)
# ============================================================================

base_theme <- theme_minimal(base_size = 8) +
  theme(
    axis.title.x      = element_text(size = 8.5, face = "bold", margin = margin(t = 4)),
    axis.title.y      = element_text(size = 7, face = "bold", margin = margin(r = 4)),
    axis.text.x       = element_text(size = 7,  color = "black",
    angle = 45, hjust = 1, vjust = 1),
    axis.text.y       = element_text(size = 7,  color = "black"),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "none",
    panel.border      = element_blank(),
    axis.line         = element_line(color = "black", linewidth = 0.4),
    axis.ticks        = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = unit(0.08, "cm"),
    plot.margin       = margin(5, 5, 5, 5),
    plot.tag          = element_text(size = 10, face = "bold"),
    plot.title        = element_text(size = 8.5, face = "bold", hjust = 0.5,
                                     margin = margin(b = 3))
  )

# Theme variant for betadisper (no x-axis title, group labels are self-explanatory)
betadisper_theme <- base_theme +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 7, color = "black"))

# ============================================================================
# SHARED I/O FUNCTIONS
# ============================================================================

load_beta_data <- function(file_path) {
  raw <- read.table(file_path, sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE, comment.char = "", check.names = FALSE)
  names(raw)[names(raw) == ""] <- "RowID"
  raw$Group1 <- gsub("NAFLD", "MASLD", raw$Group1)
  raw$Group2 <- gsub("NAFLD", "MASLD", raw$Group2)
  return(raw)
}

build_dist_matrix <- function(raw_data) {
  samples_info <- unique(rbind(
    data.frame(SampleID = raw_data$SubjectID1, Group = raw_data$Group1),
    data.frame(SampleID = raw_data$SubjectID2, Group = raw_data$Group2)
  ))
  samples_info <- samples_info[!duplicated(samples_info$SampleID), ]
  samples_info$Group <- factor(samples_info$Group,
                               levels = c("Healthy", "MASLD", "Cirrhosis", "HCC"))

  sample_ids <- sort(unique(c(raw_data$SubjectID1, raw_data$SubjectID2)))
  n <- length(sample_ids)
  mat <- matrix(0, nrow = n, ncol = n)
  rownames(mat) <- colnames(mat) <- sample_ids

  for (i in 1:nrow(raw_data)) {
    pos1 <- which(sample_ids == raw_data$SubjectID1[i])
    pos2 <- which(sample_ids == raw_data$SubjectID2[i])
    mat[pos1, pos2] <- raw_data$Distance[i]
    mat[pos2, pos1] <- raw_data$Distance[i]
  }

  dist_obj <- as.dist(mat)
  samples_info <- samples_info[match(sample_ids, samples_info$SampleID), ]

  return(list(dist = dist_obj, samples = samples_info, ids = sample_ids))
}

# ============================================================================
# PCoA PANEL FUNCTION
# ============================================================================

create_pcoa_plot <- function(dist_obj, samples_info, metric_name) {

  pcoa_result  <- pcoa(dist_obj)
  pcoa_coords  <- pcoa_result$vectors
  eigenvalues  <- pcoa_result$values$Eigenvalues
  var_exp      <- eigenvalues / sum(eigenvalues) * 100

  pcoa_df <- data.frame(
    SampleID = rownames(pcoa_coords),
    PC1      = pcoa_coords[, 1],
    PC2      = pcoa_coords[, 2],
    Group    = samples_info$Group[match(rownames(pcoa_coords), samples_info$SampleID)]
  )

  sample_counts <- pcoa_df %>% count(Group)
  label_levels  <- paste0(levels(pcoa_df$Group), " (n=",
                          sample_counts$n[match(levels(pcoa_df$Group), sample_counts$Group)], ")")
  pcoa_df$Group_label <- factor(
    paste0(pcoa_df$Group, " (n=",
           sample_counts$n[match(pcoa_df$Group, sample_counts$Group)], ")"),
    levels = label_levels
  )

  # PERMANOVA
  set.seed(123)
  perm   <- adonis2(dist_obj ~ Group, data = samples_info, permutations = 999)
  f_stat <- perm$F[1]
  r2     <- perm$R2[1]
  p_val  <- perm$`Pr(>F)`[1]

  perm_label <- sprintf("PERMANOVA: F=%.2f, R\u00B2=%.3f, p=%s",
                        f_stat, r2,
                        ifelse(p_val < 0.001, "<0.001", sprintf("%.3f", p_val)))

  y_max   <- max(pcoa_df$PC2)
  y_min   <- min(pcoa_df$PC2)
  y_range <- y_max - y_min
  x_min   <- min(pcoa_df$PC1)

  p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, fill = Group_label, color = Group_label)) +
    stat_ellipse(type = "norm", level = 0.68, linewidth = 0.5,
                 alpha = 0.9, linetype = "solid") +
    geom_point(size = 1.6, alpha = 0.85, shape = 21,
               color = "black", stroke = 0.3) +
    scale_fill_manual(values  = setNames(group_colors, label_levels)) +
    scale_color_manual(values = setNames(group_colors, label_levels)) +
    labs(
      x    = sprintf("PC1 (%.1f%%)", var_exp[1]),
      y    = sprintf("PC2 (%.1f%%)", var_exp[2]),
      fill = NULL, color = NULL
    ) +
    ggtitle(metric_name) +
    annotate("text", x = x_min, y = y_max + (y_range * 0.30),
             label = perm_label,
             hjust = 0, vjust = 1, size = 2.2, fontface = "italic") +
    coord_cartesian(ylim = c(y_min - (y_range * 0.05),
                             y_max + (y_range * 0.35)),
                    clip = "off") +
    base_theme

  return(list(plot = p, pcoa_df = pcoa_df, label_levels = label_levels))
}

# ============================================================================
# BETADISPER VIOLIN PANEL FUNCTION
# ============================================================================

create_betadisper_violin <- function(dist_info, y_label, title_text) {

  bd <- betadisper(dist_info$dist, dist_info$samples$Group)

  set.seed(123)
  perm_test  <- permutest(bd, permutations = 999)
  overall_p  <- perm_test$tab$`Pr(>F)`[1]
  overall_F  <- perm_test$tab$F[1]

  within_data <- data.frame(
    SampleID       = names(bd$distances),
    Distance       = bd$distances,
    Disease_Status = bd$group
  )
  within_data$Disease_Status <- factor(within_data$Disease_Status,
                                       levels = c("Healthy", "MASLD", "Cirrhosis", "HCC"))

  sample_counts <- within_data %>%
    group_by(Disease_Status) %>%
    summarise(n = n(), .groups = "drop")
  # X-axis labels: just disease group names (sample sizes shown in PCoA panel legends)
  status_labels <- setNames(
    as.character(sample_counts$Disease_Status),
    sample_counts$Disease_Status
  )

  # PERMDISP header label
  perm_label <- sprintf("PERMDISP: F=%.2f, p=%s",
                        overall_F,
                        ifelse(overall_p < 0.001, "<0.001",
                               sprintf("%.3f", overall_p)))

  p <- ggplot(within_data, aes(x = Disease_Status, y = Distance, fill = Disease_Status)) +
    geom_violin(alpha = 0.85, trim = TRUE, scale = "width",
                adjust = 1.1, width = 0.7, color = "black", linewidth = 0.3) +
    geom_boxplot(width = 0.13, alpha = 0.9, outlier.shape = NA,
                 color = "black", fill = "white", linewidth = 0.4) +
    geom_jitter(width = 0.07, alpha = 0.20, size = 0.25,
                color = "black", stroke = 0.05) +
    stat_summary(fun = mean, geom = "point", shape = 18,
                 size = 2.2, color = "red") +
    scale_fill_manual(values = group_colors) +
    scale_x_discrete(labels = status_labels) +
    labs(y = y_label) +
    ggtitle(title_text) +
    betadisper_theme

  # Tukey HSD pairwise — add significance bars
  tukey       <- TukeyHSD(bd)
  tukey_df    <- as.data.frame(tukey$group)
  tukey_df$comparison <- gsub("NAFLD", "MASLD", rownames(tukey_df))
  q_vals      <- tukey_df$`p adj`
  sig_idx     <- which(q_vals < 0.05)

  if (length(sig_idx) > 0) {
    max_y       <- max(within_data$Distance)
    y_range     <- max(within_data$Distance) - min(within_data$Distance)
    bar_start   <- max_y + (y_range * 0.35)
    bar_spacing <- y_range * 0.25

    group_pos <- setNames(1:4, c("Healthy", "MASLD", "Cirrhosis", "HCC"))

    for (k in seq_along(sig_idx)) {
      idx  <- sig_idx[k]
      grps <- strsplit(tukey_df$comparison[idx], "-")[[1]]
      x1   <- group_pos[trimws(grps[1])]
      x2   <- group_pos[trimws(grps[2])]
      if (x1 > x2) { tmp <- x1; x1 <- x2; x2 <- tmp }
      y_pos <- bar_start + (k - 1) * bar_spacing

      asterisks <- ifelse(q_vals[idx] < 0.001, "***",
                   ifelse(q_vals[idx] < 0.01,  "**", "*"))

      p <- p +
        annotate("segment", x = x1, xend = x2, y = y_pos, yend = y_pos,
                 color = "black", linewidth = 0.3) +
        annotate("segment", x = x1, xend = x1,
                 y = y_pos - (y_range * 0.025), yend = y_pos + (y_range * 0.025),
                 color = "black", linewidth = 0.3) +
        annotate("segment", x = x2, xend = x2,
                 y = y_pos - (y_range * 0.025), yend = y_pos + (y_range * 0.025),
                 color = "black", linewidth = 0.3) +
        annotate("text", x = (x1 + x2) / 2, y = y_pos + (y_range * 0.025),
                 label = asterisks, size = 2.5, fontface = "bold", color = "black")
    }

    new_max <- bar_start + length(sig_idx) * bar_spacing + (y_range * 0.06)
    p <- p + scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                                limits = c(NA, new_max))
  }

  # Add PERMDISP header annotation
  p <- p + labs(subtitle = perm_label) +
    theme(plot.subtitle = element_text(size = 6.8, face = "italic",
                                       hjust = 0.5, margin = margin(b = 2)))

  return(p)
}

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading beta diversity data...\n")

bc_raw  <- load_beta_data(file.path(beta_dir, "Bray-curtis",        "Bray-Curtis_raw_data.tsv"))
jac_raw <- load_beta_data(file.path(beta_dir, "Jaccard",            "Jaccardraw_data.tsv"))
uw_raw  <- load_beta_data(file.path(beta_dir, "UnWeighted_UniFrac", "Unweighted.tsv"))
wt_raw  <- load_beta_data(file.path(beta_dir, "Weighted_UniFrac",   "Weighted_rowdata.tsv"))

bc_info  <- build_dist_matrix(bc_raw)
jac_info <- build_dist_matrix(jac_raw)
uw_info  <- build_dist_matrix(uw_raw)
wt_info  <- build_dist_matrix(wt_raw)

cat("  All data loaded\n\n")

# ============================================================================
# CREATE 8 PANELS
# ============================================================================

cat("Creating PCoA panels (left column)...\n")
pcoa_bc  <- create_pcoa_plot(bc_info$dist,  bc_info$samples,  "Bray-Curtis")
pcoa_jac <- create_pcoa_plot(jac_info$dist, jac_info$samples, "Jaccard")
pcoa_uw  <- create_pcoa_plot(uw_info$dist,  uw_info$samples,  "Unweighted UniFrac")
pcoa_wt  <- create_pcoa_plot(wt_info$dist,  wt_info$samples,  "Weighted UniFrac")

cat("Creating betadisper panels (right column)...\n")
vio_bc  <- create_betadisper_violin(bc_info,  "Distance to centroid", "Bray-Curtis")
vio_jac <- create_betadisper_violin(jac_info, "Distance to centroid", "Jaccard")
vio_uw  <- create_betadisper_violin(uw_info,  "Distance to centroid", "Unweighted UniFrac")
vio_wt  <- create_betadisper_violin(wt_info,  "Distance to centroid", "Weighted UniFrac")

cat("  All 8 panels done\n\n")

# ============================================================================
# SHARED LEGEND (bottom of figure)
# ============================================================================

legend_plot <- ggplot(pcoa_bc$pcoa_df,
                      aes(x = PC1, y = PC2, fill = Group_label, color = Group_label)) +
  geom_point(size = 2.5, shape = 21, color = "black", stroke = 0.4) +
  scale_fill_manual(values  = setNames(group_colors, pcoa_bc$label_levels)) +
  scale_color_manual(values = setNames(group_colors, pcoa_bc$label_levels)) +
  labs(fill = NULL, color = NULL) +
  theme_minimal() +
  theme(
    legend.position  = "bottom",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 8, face = "bold", margin = margin(l = 2, r = 10)),
    legend.key.size  = unit(0.4, "cm"),
    legend.spacing.x = unit(0.25, "cm")
  ) +
  guides(fill = guide_legend(nrow = 1), color = guide_legend(nrow = 1))

shared_legend <- cowplot::get_legend(legend_plot)

# ============================================================================
# ASSEMBLE 4×2 GRID
# Layout: row = metric, col 1 = PCoA, col 2 = betadisper
# ============================================================================

cat("Assembling Figure 3 (4×2 grid)...\n")

panels <- (pcoa_bc$plot  | vio_bc)  /
          (pcoa_jac$plot | vio_jac) /
          (pcoa_uw$plot  | vio_uw)  /
          (pcoa_wt$plot  | vio_wt)  +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

figure3 <- cowplot::plot_grid(
  panels,
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.035)
)

# ============================================================================
# EXPORT — Gut Microbes specifications
# Width: 174 mm (double-column / full page)
# Height: 215 mm (4 rows ≈ 50 mm each + legend)
# Resolution: 600 dpi for raster; vector PDF/SVG for editable
# ============================================================================

cat("Saving Figure 3...\n")

ggsave("Figure3_BetaDiversity_Combined.pdf",
       plot = figure3,
       width = 174, height = 215, units = "mm")

ggsave("Figure3_BetaDiversity_Combined.png",
       plot = figure3,
       width = 174, height = 215, units = "mm",
       dpi = 600, bg = "white")

ggsave("Figure3_BetaDiversity_Combined.svg",
       plot = figure3,
       width = 174, height = 215, units = "mm",
       device = svglite)

cat("\nFigure 3 saved (Gut Microbes format):\n")
cat("  - Figure3_BetaDiversity_Combined.pdf  (vector, submission-ready)\n")
cat("  - Figure3_BetaDiversity_Combined.png  (600 dpi raster)\n")
cat("  - Figure3_BetaDiversity_Combined.svg  (editable vector)\n")
cat("Location:", getwd(), "\n\n")

# ============================================================================
# PRINT STATISTICAL SUMMARY (for Methods/Results/legend)
# ============================================================================

cat("=== STATISTICAL SUMMARY ===\n")
for (m in list(
  list(info = bc_info,  name = "Bray-Curtis"),
  list(info = jac_info, name = "Jaccard"),
  list(info = uw_info,  name = "Unweighted UniFrac"),
  list(info = wt_info,  name = "Weighted UniFrac")
)) {
  cat(sprintf("\n--- %s ---\n", m$name))

  # PERMANOVA
  set.seed(123)
  perm <- adonis2(m$info$dist ~ Group, data = m$info$samples, permutations = 999)
  cat(sprintf("  PERMANOVA:  F=%.2f, R\u00B2=%.3f, p=%.4f\n",
              perm$F[1], perm$R2[1], perm$`Pr(>F)`[1]))

  # PERMDISP
  bd <- betadisper(m$info$dist, m$info$samples$Group)
  set.seed(123)
  pd <- permutest(bd, permutations = 999)
  cat(sprintf("  PERMDISP:   F=%.2f, p=%.4f\n",
              pd$tab$F[1], pd$tab$`Pr(>F)`[1]))

  # Tukey HSD pairwise
  tk <- TukeyHSD(bd)
  tk_df <- as.data.frame(tk$group)
  for (i in 1:nrow(tk_df)) {
    sig <- ifelse(tk_df$`p adj`[i] < 0.05, " *", "")
    cat(sprintf("    %s: diff=%.4f, p.adj=%.4f%s\n",
                rownames(tk_df)[i], tk_df$diff[i], tk_df$`p adj`[i], sig))
  }
}

cat("\nDone.\n")
