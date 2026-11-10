################################################################################
# Supplementary Figure S1: Rarefaction Curves
# Observed Features and Shannon Diversity vs Sequencing Depth
# Styled to match Main Figures 1-8
#
# INPUT FILES (update paths below):
#   - qiime2_analysis_v2/revision_exports/asv_table/asv_table.tsv (ASV table, filtered, unrarefied)
#   - 16S_reprocess/sample_metadata.tsv
#
# OUTPUT: FigS1_Rarefaction_Curves.pdf/.tiff/.png/.svg (7 x 4 inches)
################################################################################

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)

# ── Paths: relative to this script's folder (outputs saved here) ────────────
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
find_proj_dir <- function(name = "16S_reprocess") {
  d <- normalizePath(getwd())
  repeat {
    cand <- file.path(d, name); if (dir.exists(cand)) return(cand)
    up <- dirname(d); if (up == d) break; d <- up
  }
  stop("Could not find a '", name, "' folder above ", getwd())
}
proj_dir <- find_proj_dir()
data_dir <- file.path(proj_dir, "qiime2_analysis_v2")

# ASV-level table (filtered, UNRAREFIED). The previous version used the GENUS
# table (ANCOM-L6/exported_genus), so its curves counted genera, not ASVs.
feature_table_path <- file.path(data_dir, "revision_exports", "asv_table", "asv_table.tsv")
if (!file.exists(feature_table_path)) feature_table_path <- file.path(data_dir, "feature_table.tsv")
metadata_path      <- file.path(proj_dir, "sample_metadata.tsv")
output_dir         <- getwd()

RAREFACTION_DEPTH <- 13000   # depth used for all diversity analyses (min library = 13,241)

for (f in c(feature_table_path, metadata_path)) if (!file.exists(f))
  stop("Missing input file: ", f, "\nRun QIIME2 Step 6d (ASV table export) first.")
cat("All input files located.\n")

# ── 1. Load Data ─────────────────────────────────────────────────────────────

otu_raw <- read.delim(feature_table_path, skip = 1, row.names = 1, check.names = FALSE)
otu <- as.data.frame(t(otu_raw))  # samples x ASVs

meta <- read.delim(metadata_path, comment.char = "", check.names = FALSE)
names(meta)[1] <- "id"
meta <- meta[!grepl("^#", meta$id), ]
rownames(meta) <- meta$id

common <- intersect(rownames(otu), rownames(meta))
otu <- otu[common, ]
meta <- meta[common, ]

cat("Samples:", nrow(otu), "| ASVs:", ncol(otu), "\n")
cat("Depth range:", min(rowSums(otu)), "-", max(rowSums(otu)), "\n")

# ── 2. Disease Group Setup ───────────────────────────────────────────────────

# Rename NAFLD → MASLD to match current nomenclature
meta$Status[meta$Status == "NAFLD"] <- "MASLD"

meta$Group <- factor(meta$Status,
                     levels = c("Healthy", "MASLD", "Cirrhosis", "HCC"))

group_colors <- c("Healthy"   = "#4DAF4A",
                  "MASLD"     = "#377EB8",
                  "Cirrhosis" = "#984EA3",
                  "HCC"       = "#E41A1C")

group_n <- table(meta$Group)

# ── 3. Compute Rarefaction ──────────────────────────────────────────────────

min_depth <- min(rowSums(otu))
# Curves stop where every group still retains >= 75% of its samples
# (= smallest group-wise 25th percentile of library size). Beyond that, group means
# are computed on a shrinking, non-random subset of samples and the curves jump.
lib <- rowSums(otu)
max_depth <- min(tapply(lib, meta$Group, quantile, probs = 0.25))
cat(sprintf("Curves computed to %s reads (>= 75%% of every group retained; max library = %s)\n",
            format(round(max_depth), big.mark = ","), format(max(lib), big.mark = ",")))

depths <- sort(unique(round(c(
  seq(100, min_depth, length.out = 15),
  seq(min_depth, max_depth, length.out = 20)
))))

cat("Computing rarefaction at", length(depths), "depths...\n")

rarefy_at_depth <- function(otu_mat, depth) {
  keep <- rowSums(otu_mat) >= depth
  otu_sub <- otu_mat[keep, ]
  
  obs <- rarefy(otu_sub, sample = depth, se = FALSE)
  
  set.seed(42)
  rarefied <- rrarefy(otu_sub, sample = depth)
  shan <- diversity(rarefied, index = "shannon", base = 2)   # log2, as QIIME 2 (Figure 2)
  
  data.frame(
    sample_id = rownames(otu_sub),
    depth = depth,
    observed = as.numeric(obs),
    shannon = as.numeric(shan)
  )
}

rare_df <- bind_rows(lapply(depths, function(d) {
  cat("  Depth:", format(d, big.mark = ","), "\n")
  rarefy_at_depth(otu, d)
}))

rare_df <- rare_df %>%
  left_join(meta %>% select(id, Group), by = c("sample_id" = "id"))

# ── 4. Summarize by Group ───────────────────────────────────────────────────

rare_summary <- rare_df %>%
  group_by(depth, Group) %>%
  summarise(
    obs_mean  = mean(observed, na.rm = TRUE),
    obs_se    = sd(observed, na.rm = TRUE) / sqrt(n()),
    shan_mean = mean(shannon, na.rm = TRUE),
    shan_se   = sd(shannon, na.rm = TRUE) / sqrt(n()),
    n_samples = n(),
    .groups = "drop"
  )

# ── 5. Publication Theme ─────────────────────────────────────────────────────

pub_theme <- theme_minimal(base_size = 9) +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    panel.grid        = element_blank(),
    axis.line         = element_line(color = "black", linewidth = 0.4),
    axis.ticks        = element_line(color = "black", linewidth = 0.3),
    axis.text         = element_text(size = 8, color = "black"),
    axis.title        = element_text(size = 10, face = "bold"),
    plot.tag          = element_text(size = 12, face = "bold"),
    legend.position   = "none",
    plot.margin       = margin(10, 15, 10, 10)
  )

# ── 6. Panel A: Observed Features ────────────────────────────────────────────

pA <- ggplot(rare_summary, aes(x = depth, y = obs_mean, 
                                color = Group, fill = Group)) +
  geom_ribbon(aes(ymin = obs_mean - obs_se, ymax = obs_mean + obs_se),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = RAREFACTION_DEPTH, linetype = "dashed", 
             color = "gray40", linewidth = 0.4) +
  annotate("text", x = RAREFACTION_DEPTH * 1.08, 
           y = max(rare_summary$obs_mean) * 0.15,
           label = paste0("Rarefaction depth\n(", 
                          format(RAREFACTION_DEPTH, big.mark = ","), ")"),
           size = 2.5, color = "gray30", hjust = 0) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Sequencing Depth", y = "Observed Features", tag = "A") +
  pub_theme

# ── 7. Panel B: Shannon Diversity ────────────────────────────────────────────

pB <- ggplot(rare_summary, aes(x = depth, y = shan_mean, 
                                color = Group, fill = Group)) +
  geom_ribbon(aes(ymin = shan_mean - shan_se, ymax = shan_mean + shan_se),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = RAREFACTION_DEPTH, linetype = "dashed", 
             color = "gray40", linewidth = 0.4) +
  annotate("text", x = RAREFACTION_DEPTH * 1.08,
           y = min(rare_summary$shan_mean) + 
             diff(range(rare_summary$shan_mean)) * 0.1,
           label = paste0("Rarefaction depth\n(", 
                          format(RAREFACTION_DEPTH, big.mark = ","), ")"),
           size = 2.5, color = "gray30", hjust = 0) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Sequencing Depth", y = "Shannon Diversity Index", tag = "B") +
  pub_theme

# ── 8. Shared Legend ─────────────────────────────────────────────────────────

legend_plot <- ggplot(rare_summary, aes(x = depth, y = obs_mean, color = Group)) +
  geom_line() +
  scale_color_manual(
    values = group_colors,
    labels = paste0(names(group_n), " (n=", group_n, ")")
  ) +
  pub_theme +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 9, margin = margin(l = 2, r = 10)))

shared_legend <- get_legend(legend_plot)

# ── 9. Combine & Save ───────────────────────────────────────────────────────

combined <- (pA | pB) / 
  ggdraw(shared_legend) +
  plot_layout(heights = c(1, 0.08))

# Vector PDF with embedded fonts (Taylor & Francis requirement)
ggsave(file.path(output_dir, "FigS1_Rarefaction_Curves.pdf"), combined,
       width = 7, height = 4, device = cairo_pdf)

# High-resolution TIFF (600 dpi, LZW-compressed) for submission
ggsave(file.path(output_dir, "FigS1_Rarefaction_Curves.tiff"), combined,
       width = 7, height = 4, dpi = 600, compression = "lzw")

ggsave(file.path(output_dir, "FigS1_Rarefaction_Curves.png"), combined,
       width = 7, height = 4, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "FigS1_Rarefaction_Curves.svg"), combined,
       width = 7, height = 4)

cat("\n✅ Saved FigS1_Rarefaction_Curves.pdf (vector, fonts embedded) and .tiff (600 dpi)\n")
cat("Dimensions: 7 x 4 inches\n")
