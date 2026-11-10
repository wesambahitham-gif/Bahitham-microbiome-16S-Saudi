#!/usr/bin/env Rscript
# ============================================================================
# FIGURE 2: ALPHA DIVERSITY
# Paper 1 — 16S rRNA Gut Microbiome | Saudi Arabian Liver Disease Cohort
# ----------------------------------------------------------------------------
# Layout (5 panels, 1x5 horizontal row):
#   A. Chao1 Richness Index
#   B. Observed Features
#   C. Shannon Diversity Index
#   D. Simpson's Diversity Index
#   E. Faith's Phylogenetic Diversity
#
# Each panel: violin + boxplot + jitter + mean (red diamond)
# Statistics: Kruskal-Wallis per metric + Dunn's post-hoc (BH-adjusted)
# Significance bars only shown when cross-metric BH-adjusted KW q < 0.05 AND Dunn p.adj < 0.05.
# Cross-metric BH correction also reported in console.
#
# Generates: Figure2_Alpha_Diversity.pdf/.png/.svg
#
# JOURNAL SPECS — Gut Microbes (Taylor & Francis):
#   - Width: 174 mm (full page)
#   - Height: 90 mm (proportional for 1x5 horizontal layout)
#   - Resolution: 600 dpi (PNG); vector PDF; SVG editable
#   - Theme matches Figures 3, 4, 5, 6, 7, 8 for cross-figure consistency
# ============================================================================

rm(list = ls())
cat("\014")

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(rstatix)
  library(patchwork)
  library(svglite)
})

# ============================================================================
# PATHS
# ============================================================================

# Working directory = folder containing this script; outputs are saved here.
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

# Input data, RELATIVE to this script's folder
#   script : Desktop/Reviewers/Figure_2/
#   data   : Desktop/16S_reprocess/qiime2_analysis_v2/
data_dir  <- file.path("..", "..", "16S_reprocess", "qiime2_analysis_v2")
alpha_dir <- file.path(data_dir, "Alpha_diversity")

# ============================================================================
# DEFENSIVE FILE-EXISTENCE CHECK
# ============================================================================

required_files <- c(
  chao1    = file.path(alpha_dir, "Chao1", "Chao1.tsv"),
  observed = file.path(alpha_dir, "Observed feature", "Observed_feature_metadata.tsv"),
  shannon  = file.path(alpha_dir, "Shannon", "Shannon_Liver.tsv"),
  simpson  = file.path(alpha_dir, "Simpson", "Simpson_liver_micorbiome.tsv"),
  faith_pd = file.path(alpha_dir, "Faith-pd", "faith-pd.tsv")
)
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop("Missing input files:\n  ",
       paste(names(missing), missing, sep = ": ", collapse = "\n  "),
       "\n\nCheck data_dir (relative to this script's folder): ", normalizePath(data_dir, mustWork = FALSE))
}
cat("All input files located.\n\n")

# ============================================================
# CONSISTENT COLOR PALETTE (same across all figures)
# ============================================================

group_colors <- c(
  "Healthy"   = "#4DAF4A",
  "MASLD"     = "#377EB8",
  "Cirrhosis" = "#984EA3",
  "HCC"       = "#E41A1C"
)

# ============================================================================
# PUBLICATION THEME (matches Figures 3, 4, 5, 6, 7, 8)
# ============================================================================

pub_theme <- theme_minimal(base_size = 8) +
  theme(
    axis.title.x      = element_blank(),
    axis.title.y      = element_text(size = 8.5, face = "bold", margin = margin(r = 4)),
    axis.text.x       = element_text(size = 7, color = "black",
    angle = 45, hjust = 1, vjust = 1),    
    axis.text.y       = element_text(size = 7, color = "black"),
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
    plot.tag          = element_text(size = 10, face = "bold")
  )

# ============================================================
# FUNCTIONS
# ============================================================

load_alpha_data <- function(file_path, column_name) {
  data <- read.table(file_path, sep = "\t", header = TRUE,
                     stringsAsFactors = FALSE, comment.char = "", check.names = FALSE)
  data <- data[data$id != "#q2:types", ]
  data$Status <- gsub("NAFLD", "MASLD", data$Status)
  data$Value <- as.numeric(data[[column_name]])
  data$Status <- factor(data$Status, levels = c("Healthy", "MASLD", "Cirrhosis", "HCC"))
  data <- data[!is.na(data$Value), ]
  return(data)
}

create_alpha_plot <- function(data, y_label, y_breaks = waiver(), show_sig = TRUE) {

  # X-axis labels: just disease group names (sample sizes go in figure caption)
  status_labels <- setNames(
    as.character(levels(data$Status)),
    levels(data$Status)
  )

  kw_test <- kruskal.test(Value ~ Status, data = data)

  p <- ggplot(data, aes(x = Status, y = Value, fill = Status)) +
    geom_violin(alpha = 0.9, trim = TRUE, scale = "width",
                adjust = 1.1, width = 0.7, color = "black", linewidth = 0.4) +
    geom_boxplot(width = 0.14, alpha = 0.9, outlier.shape = NA,
                 color = "black", fill = "white", linewidth = 0.6) +
    geom_jitter(width = 0.1, alpha = 0.35, size = 0.8, color = "black", stroke = 0.1) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3.5, color = "red") +
    scale_fill_manual(values = group_colors) +
    scale_x_discrete(labels = status_labels) +
    labs(y = y_label) +
    pub_theme

  if (!identical(y_breaks, waiver())) {
    p <- p + scale_y_continuous(breaks = y_breaks,
                                expand = expansion(mult = c(0.02, 0.05)))
  }

  # Significance bars
  if (show_sig && kw_test$p.value < 0.05) {
    dunn_result <- data %>% dunn_test(Value ~ Status, p.adjust.method = "BH")
    sig_dunn <- dunn_result[dunn_result$p.adj < 0.05, ]

    if (nrow(sig_dunn) > 0) {
      max_y <- max(data$Value)
      y_range <- max(data$Value) - min(data$Value)
      bar_start <- max_y + (y_range * 0.30)
      bar_spacing <- y_range * 0.15

      group_pos <- setNames(1:4, c("Healthy", "MASLD", "Cirrhosis", "HCC"))

      for (i in 1:nrow(sig_dunn)) {
        x1 <- group_pos[sig_dunn$group1[i]]
        x2 <- group_pos[sig_dunn$group2[i]]
        y_pos <- bar_start + (i - 1) * bar_spacing

        asterisks <- ifelse(sig_dunn$p.adj[i] < 0.001, "***",
                            ifelse(sig_dunn$p.adj[i] < 0.01, "**", "*"))

        p <- p +
          annotate("segment", x = x1, xend = x2, y = y_pos, yend = y_pos,
                   color = "black", linewidth = 0.3) +
          annotate("segment", x = x1, xend = x1,
                   y = y_pos - (y_range * 0.015), yend = y_pos + (y_range * 0.015),
                   color = "black", linewidth = 0.3) +
          annotate("segment", x = x2, xend = x2,
                   y = y_pos - (y_range * 0.015), yend = y_pos + (y_range * 0.015),
                   color = "black", linewidth = 0.3) +
          annotate("text", x = (x1 + x2) / 2, y = y_pos + (y_range * 0.03),
                   label = asterisks, size = 3.2, fontface = "bold", color = "black")
      }

      new_max <- bar_start + nrow(sig_dunn) * bar_spacing + (y_range * 0.06)
      if (!identical(y_breaks, waiver())) {
        p <- p + scale_y_continuous(breaks = y_breaks,
                                    expand = expansion(mult = c(0.02, 0.02)),
                                    limits = c(NA, new_max))
      } else {
        p <- p + scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)),
                                    limits = c(NA, new_max))
      }
    }
  }

  return(p)
}

# ============================================================
# LOAD DATA
# ============================================================

cat("Loading alpha diversity data...\n")

chao1_data    <- load_alpha_data(file.path(alpha_dir, "Chao1", "Chao1.tsv"), "chao1")
observed_data <- load_alpha_data(file.path(alpha_dir, "Observed feature", "Observed_feature_metadata.tsv"), "observed_features")
shannon_data  <- load_alpha_data(file.path(alpha_dir, "Shannon", "Shannon_Liver.tsv"), "shannon_entropy")

# Simpson - special column handling
simpson_raw <- read.table(file.path(alpha_dir, "Simpson", "Simpson_liver_micorbiome.tsv"),
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                          comment.char = "", check.names = FALSE)
simpson_raw <- simpson_raw[simpson_raw$id != "#q2:types", ]
simpson_raw$Status <- gsub("NAFLD", "MASLD", simpson_raw$Status)
simp_col <- grep("^simpson$|simpson", names(simpson_raw), value = TRUE)
simp_col <- simp_col[length(simp_col)]   # QIIME 2 names it "simpson" or "functools.partial(<..._simpsons_dominance..>)"
simpson_raw$Value <- as.numeric(simpson_raw[[simp_col]])
simpson_raw$Status <- factor(simpson_raw$Status, levels = c("Healthy", "MASLD", "Cirrhosis", "HCC"))
simpson_data <- simpson_raw[!is.na(simpson_raw$Value), ]

# Faith PD - load using helper (column name is faith_pd, lowercase)
faith_pd_data <- load_alpha_data(file.path(alpha_dir, "Faith-pd", "faith-pd.tsv"), "faith_pd")

cat("  All data loaded\n\n")

# ============================================================
# CREATE PANELS
# ============================================================

cat("Creating panels...\n")

# Cross-metric BH correction (5 KW tests) gates the significance bars, so the
# figure matches the text: bars are drawn only if the metric's BH-adjusted KW q < 0.05.
kw_q <- p.adjust(c(
  Chao1    = kruskal.test(Value ~ Status, data = chao1_data)$p.value,
  Observed = kruskal.test(Value ~ Status, data = observed_data)$p.value,
  Shannon  = kruskal.test(Value ~ Status, data = shannon_data)$p.value,
  Simpson  = kruskal.test(Value ~ Status, data = simpson_data)$p.value,
  Faith    = kruskal.test(Value ~ Status, data = faith_pd_data)$p.value), method = "BH")

# y-axis breaks are data-driven (pretty breaks) so they fit the reprocessed data
p_chao1    <- create_alpha_plot(chao1_data,    "Chao1 Richness Index",           show_sig = kw_q["Chao1"]    < 0.05)
p_observed <- create_alpha_plot(observed_data, "Observed Features",              show_sig = kw_q["Observed"] < 0.05)
p_shannon  <- create_alpha_plot(shannon_data,  "Shannon Diversity Index",        show_sig = kw_q["Shannon"]  < 0.05)
p_simpson  <- create_alpha_plot(simpson_data,  "Simpson's Diversity Index",      show_sig = kw_q["Simpson"]  < 0.05)
p_faith    <- create_alpha_plot(faith_pd_data, "Faith's Phylogenetic Diversity", show_sig = kw_q["Faith"]    < 0.05)

# ============================================================
# COMBINE: 1x5 horizontal row
# ============================================================

cat("Assembling Figure 2...\n")

figure2 <- (p_chao1 + p_observed + p_shannon) /
            (p_simpson + p_faith + plot_spacer()) +
  plot_layout(nrow = 2) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

# ============================================================================
# SAVE — Gut Microbes specifications
# Width: 174 mm | Height: 150 mm | 600 dpi PNG | vector PDF/SVG
# 2x3 grid layout: top row [A B C], bottom row [D E spacer]
# Each panel ~55 mm wide; horizontal x-axis labels fit cleanly
# ============================================================================

cat("Saving Figure 2...\n")

ggsave("Figure2_Alpha_Diversity.pdf",
       plot = figure2, width = 174, height = 150, units = "mm")

ggsave("Figure2_Alpha_Diversity.png",
       plot = figure2, width = 174, height = 150, units = "mm",
       dpi = 600, bg = "white")

ggsave("Figure2_Alpha_Diversity.svg",
       plot = figure2, width = 174, height = 150, units = "mm",
       device = svglite)

cat("\nFigure 2 saved (Gut Microbes format):\n")
cat("  - Figure2_Alpha_Diversity.pdf  (vector, submission-ready)\n")
cat("  - Figure2_Alpha_Diversity.png  (600 dpi raster)\n")
cat("  - Figure2_Alpha_Diversity.svg  (editable vector)\n")
cat("Location:", getwd(), "\n")


cat("\n=== Kruskal-Wallis Results (per metric) ===\n")
datasets <- list(
  list(data = chao1_data,    name = "Chao1"),
  list(data = observed_data, name = "Observed Features"),
  list(data = shannon_data,  name = "Shannon"),
  list(data = simpson_data,  name = "Simpson"),
  list(data = faith_pd_data, name = "Faith PD")
)

# Collect raw p-values dynamically (no hardcoded values)
kw_pvals <- numeric(length(datasets))
names(kw_pvals) <- sapply(datasets, function(m) m$name)
for (i in seq_along(datasets)) {
  kw <- kruskal.test(Value ~ Status, data = datasets[[i]]$data)
  kw_pvals[i] <- kw$p.value
  cat(sprintf("  %s: H(3) = %.2f, p = %.4f\n",
              datasets[[i]]$name, kw$statistic, kw$p.value))
}

# Cross-metric BH correction (5 tests)
p_adj <- p.adjust(kw_pvals, method = "BH")
cat("\n=== BH-Adjusted Across 5 Metrics ===\n")
for (i in seq_along(kw_pvals)) {
  cat(sprintf("  %s: raw p = %.4f, p.adj = %.4f\n",
              names(kw_pvals)[i], kw_pvals[i], p_adj[i]))
}
cat("\nDone.\n")
