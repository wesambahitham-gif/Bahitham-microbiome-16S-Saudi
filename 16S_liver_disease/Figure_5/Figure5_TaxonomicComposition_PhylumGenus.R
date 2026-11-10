#!/usr/bin/env Rscript
# ============================================================================
# FIGURE 5: TAXONOMIC COMPOSITION — PHYLUM + GENUS LEVELS  (REVISED v4: F/B median 95% CI)
# Paper 1 — 16S rRNA Gut Microbiome | Saudi Arabian Liver Disease Cohort
# ----------------------------------------------------------------------------
# Layout (4 panels, 2x2 grid + shared bottom legend):
#   Row 1 (Phylum, top 6 + Other):
#     A. Mean stacked bar per disease group
#     B. Horizontal bar chart with SE error bars + KW(BH) significance
#   Row 2 (Genus, top 15 + Other):
#     C. Mean stacked bar per disease group
#     D. Horizontal bar chart with SE error bars + KW(BH) significance
#
# REVISIONS in this version (v3):
#   1. X-axis labels in A & C: 45-degree rotation explicitly preserved
#   2. "Other" in legend: drop=FALSE forces visibility even at <1% abundance
#   3. Genus legend (C): 2-column layout to prevent over-tall legend
#   4. Spacing between rows: explicit spacer with proportional height
#   5. Bar chart x-axis ranges normalized: B uses 0-80%, D uses 0-30%
#   6. Significance asterisks: padded right margin to prevent edge clipping
#   7. Single shared "Disease Status" legend at figure bottom
#   8. Genus color #6: replaced #C71585 -> #8DD3C7 to avoid HCC red collision
#   9. Wider error bar caps for visibility at small print size
#
# JOURNAL SPECS — Microbiology Spectrum (ASM):
#   - Width: 174 mm (full page)
#   - Height: 220 mm
#   - Resolution: 600 dpi PNG; vector PDF; SVG editable
# ============================================================================

rm(list = ls())
cat("\014")

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(cowplot)
  library(svglite)
})

# ============================================================================
# 1. PATHS
# ============================================================================

# Paths are RELATIVE to this script's folder (outputs are saved here)
#   script : 16S_liver_disease/Figure_5/
#   data   : Desktop/16S_reprocess/qiime2_analysis_v2/
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

data_dir      <- file.path("..", "..", "16S_reprocess", "qiime2_analysis_v2")
taxa_file     <- file.path(data_dir, "Taxanomy", "level-6(1).csv")
metadata_file <- file.path("..", "..", "16S_reprocess", "sample_metadata.tsv")
output_dir    <- "."

# ============================================================================
# 2. DEFENSIVE FILE-EXISTENCE CHECK
# ============================================================================

required_files <- c(taxa = taxa_file, metadata = metadata_file)
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop("Missing input files:\n  ",
       paste(names(missing), missing, sep = ": ", collapse = "\n  "),
       "\n\nWorking directory: ", getwd(),
       "\nSave this script in 16S_liver_disease/Figure_5/ and click Source.")
}
cat("All input files located.\n\n")

# ============================================================================
# 3. COLOR PALETTES
# ============================================================================

group_colors <- c(
  "Healthy"   = "#4DAF4A",
  "MASLD"     = "#377EB8",
  "Cirrhosis" = "#984EA3",
  "HCC"       = "#E41A1C"
)

phyla_palette <- c(
  "#440154",  # Firmicutes
  "#3B528B",  # Bacteroidota
  "#21908C",  # Proteobacteria
  "#5DC863",  # Actinobacteriota
  "#FDE725",  # Verrucomicrobiota
  "#B12A90"   # Desulfobacterota
)

# Genera — color #6 = mint (#8DD3C7) NOT magenta (avoids HCC red collision)
genera_palette <- c(
  "#A6CEE3", "#1F78B4",
  "#B2DF8A", "#33A02C",
  "#FB9A99", "#8DD3C7",
  "#FDBF6F", "#FF7F00",
  "#CAB2D6", "#6A3D9A",
  "#FFFF99", "#B15928",
  "#80B1D3",
  "#BEBADA",
  "#FCCDE5"
)

# ============================================================================
# 4. THEME (matches Figures 3, 4)
# ============================================================================

base_theme <- theme_minimal(base_size = 8) +
  theme(
    axis.title.x       = element_text(size = 7.5, face = "bold",
                                      margin = margin(t = 0, b = 0),
                                      vjust = 1),
    axis.title.y       = element_text(size = 8.5, face = "bold", margin = margin(r = 4)),
    axis.text.x        = element_text(size = 7,  color = "black"),
    axis.text.y        = element_text(size = 7,  color = "black"),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.grid.minor   = element_blank(),
    panel.border       = element_blank(),
    axis.line          = element_line(color = "black", linewidth = 0.4),
    axis.ticks         = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length  = unit(0.08, "cm"),
    legend.title       = element_text(size = 7,   face = "bold"),
    legend.text        = element_text(size = 6.5, face = "italic"),
    legend.key.size    = unit(0.3, "cm"),
    legend.key.height  = unit(0.3, "cm"),
    plot.margin        = margin(5, 5, 5, 5),
    plot.tag           = element_text(size = 10, face = "bold")
  )

# ============================================================================
# 5. LOAD & PROCESS
# ============================================================================

cat("Loading taxonomic data...\n")
taxa_data <- read.csv(taxa_file, row.names = 1, check.names = FALSE)
metadata  <- read.delim(metadata_file, row.names = 1, check.names = FALSE)

extract_phylum <- function(taxa_name) {
  parts <- strsplit(taxa_name, ";")[[1]]
  phylum_part <- parts[grepl("p__", parts)]
  if (length(phylum_part) > 0) {
    phylum <- gsub("p__", "", trimws(phylum_part[1]))
    phylum <- gsub("_+$", "", phylum)
    if (phylum == "" || grepl("^[Uu]n(classified|cultured|known)", phylum)) return(NA)
    return(phylum)
  }
  return(NA)
}

extract_genus <- function(taxa_name) {
  parts <- strsplit(taxa_name, ";")[[1]]
  genus_part <- parts[grepl("g__", parts)]
  if (length(genus_part) > 0) {
    genus <- gsub("g__", "", genus_part[1])
    genus <- gsub("_+$", "", genus)
    if (genus == "" || grepl("^[Uu]n(classified|cultured|known)", genus) || genus == "metagenome") return(NA)
    return(genus)
  }
  return(NA)
}

# Phylum aggregation
phylum_names <- sapply(colnames(taxa_data), extract_phylum)
phy_classified <- !is.na(phylum_names)
phylum_names <- phylum_names[phy_classified]
taxa_phy <- taxa_data[, phy_classified]

phylum_data <- data.frame(row.names = rownames(taxa_phy))
for (phy in unique(phylum_names)) {
  cols <- names(phylum_names)[phylum_names == phy]
  phylum_data[[phy]] <- if (length(cols) == 1) taxa_phy[[cols]] else rowSums(taxa_phy[, cols, drop = FALSE])
}

# Genus aggregation
genus_names <- sapply(colnames(taxa_data), extract_genus)
gen_classified <- !is.na(genus_names)
genus_names <- genus_names[gen_classified]
taxa_gen <- taxa_data[, gen_classified]

genus_data <- data.frame(row.names = rownames(taxa_gen))
for (g in unique(genus_names)) {
  cols <- names(genus_names)[genus_names == g]
  genus_data[[g]] <- if (length(cols) == 1) taxa_gen[[cols]] else rowSums(taxa_gen[, cols, drop = FALSE])
}

top6  <- names(sort(colSums(phylum_data), decreasing = TRUE))[1:6]
top15 <- names(sort(colSums(genus_data),  decreasing = TRUE))[1:15]

cat("Top 6 phyla:\n")
for (i in seq_along(top6))  cat(sprintf("  %d. %s\n", i, top6[i]))
cat("Top 15 genera:\n")
for (i in seq_along(top15)) cat(sprintf("  %2d. %s\n", i, top15[i]))

common_phy <- intersect(rownames(phylum_data), rownames(metadata))
common_gen <- intersect(rownames(genus_data),  rownames(metadata))
phylum_matched <- phylum_data[common_phy, ]
genus_matched  <- genus_data[common_gen, ]
meta_phy       <- metadata[common_phy, ]
meta_gen       <- metadata[common_gen, ]

disease_phy <- factor(gsub("NAFLD", "MASLD", meta_phy$Status),
                      levels = c("Healthy", "MASLD", "Cirrhosis", "HCC"))
disease_gen <- factor(gsub("NAFLD", "MASLD", meta_gen$Status),
                      levels = c("Healthy", "MASLD", "Cirrhosis", "HCC"))

rel_phy <- phylum_matched / rowSums(phylum_matched) * 100
rel_gen <- genus_matched  / rowSums(genus_matched)  * 100

phy_df <- rel_phy[, top6]
phy_df$Other <- 100 - rowSums(phy_df)
phy_df$Disease_Status <- disease_phy

gen_df <- rel_gen[, top15]
gen_df$Other <- 100 - rowSums(gen_df)
gen_df$Disease_Status <- disease_gen

phy_cols <- c(top6, "Other")
gen_cols <- c(top15, "Other")

n_phy <- table(disease_phy)
n_gen <- table(disease_gen)
# X-axis labels: just disease group names (sample sizes appear in shared bottom legend)
xlab_phy <- names(n_phy)
xlab_gen <- names(n_gen)

cat("\n'Other' phylum abundance (mean % per group):\n")
phy_other_means <- tapply(phy_df$Other, phy_df$Disease_Status, mean)
print(round(phy_other_means, 3))

cat("\n'Other' genus abundance (mean % per group):\n")
gen_other_means <- tapply(gen_df$Other, gen_df$Disease_Status, mean)
print(round(gen_other_means, 3))

# DECISION: If 'Other' phylum is consistently <1% across all groups,
# drop it from Panel A entirely (the slice would be invisible anyway,
# and showing a phantom legend entry confuses readers). Threshold = 1%.
# Genus 'Other' typically >2% so always retained.
phy_other_max <- max(phy_other_means, na.rm = TRUE)
keep_phy_other <- phy_other_max >= 1.0

if (!keep_phy_other) {
  cat(sprintf("\nNote: Phylum 'Other' max = %.2f%% across groups (<1%% threshold).\n",
              phy_other_max))
  cat("  Dropping 'Other' from Panel A — top 6 phyla account for ~100%.\n")
  phy_df$Other <- NULL
  phy_cols <- top6  # remove "Other" from levels
} else {
  cat(sprintf("\nPhylum 'Other' max = %.2f%% — retained in Panel A.\n", phy_other_max))
}

gen_other_max <- max(gen_other_means, na.rm = TRUE)
cat(sprintf("Genus 'Other' max = %.2f%% — retained in Panel C.\n", gen_other_max))

# ============================================================================
# 6. CLEAN LABELS
# ============================================================================

make_clean_label <- function(name) {
  out <- name
  out <- gsub("^Prevotella_9$", "Prevotella", out)
  out <- gsub("^Prevotella_7$", "Prevotella", out)
  out <- gsub("^UCG-002$", "Oscillospiraceae UCG-002", out)   # SILVA 138: f__Oscillospiraceae
  out <- gsub("^UCG-005$", "Oscillospiraceae UCG-005", out)
  out <- gsub("^\\[Eubacterium\\]_hallii_group$", "[Eubacterium] hallii group", out)
  out <- gsub("^\\[Eubacterium\\]_coprostanoligenes_group$",
              "[Eubacterium] coprostanoligenes group", out)
  out <- gsub("_", " ", out)
  return(out)
}

phy_clean_labels <- setNames(sapply(phy_cols, make_clean_label), phy_cols)
if ("Other" %in% phy_cols) phy_clean_labels["Other"] <- "Other"
gen_clean_labels <- setNames(sapply(gen_cols, make_clean_label), gen_cols)
if ("Other" %in% gen_cols) gen_clean_labels["Other"] <- "Other"

cat("\nData processed.\n\n")

# ============================================================================
# 7. PANEL A — PHYLUM STACKED
# ============================================================================

cat("Building Panel A...\n")

phy_mean <- phy_df %>%
  group_by(Disease_Status) %>%
  summarise(across(all_of(phy_cols), mean), .groups = "drop") %>%
  pivot_longer(cols = all_of(phy_cols), names_to = "Phylum", values_to = "Mean_Abundance") %>%
  mutate(Phylum = factor(Phylum, levels = phy_cols))

phy_stacked_colors <- setNames(phyla_palette[1:6], top6)
if ("Other" %in% phy_cols) phy_stacked_colors["Other"] <- "#CCCCCC"

panel_A <- ggplot(phy_mean, aes(x = Disease_Status, y = Mean_Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2, width = 0.65) +
  scale_fill_manual(values = phy_stacked_colors, labels = phy_clean_labels, drop = FALSE) +
  scale_x_discrete(labels = xlab_phy) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)), breaks = seq(0, 100, 25)) +
  labs(x = NULL, y = "Relative Abundance (%)", fill = "Phylum", tag = "A") +
  base_theme +
  theme(
    axis.text.x        = element_text(size = 7, face = "bold", color = "black",
                                      angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#E5E5E5", linewidth = 0.3),
    legend.position    = "right"
  ) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE,
                             keywidth = unit(0.3, "cm"),
                             keyheight = unit(0.3, "cm")))

# ============================================================================
# 8. PANEL B — PHYLUM BAR CHART
# ============================================================================

cat("Building Panel B...\n")

phy_b <- phy_df[, c(top6, "Disease_Status")] %>%
  pivot_longer(cols = all_of(top6), names_to = "Phylum", values_to = "Abundance") %>%
  group_by(Disease_Status, Phylum) %>%
  summarise(Mean = mean(Abundance), SE = sd(Abundance) / sqrt(n()), .groups = "drop")

phy_order <- phy_b %>% group_by(Phylum) %>%
  summarise(Overall = mean(Mean), .groups = "drop") %>%
  arrange(Overall) %>% pull(Phylum)
phy_b$Phylum <- factor(phy_b$Phylum, levels = phy_order)
phy_y_clean <- setNames(sapply(phy_order, make_clean_label), phy_order)

phy_raw_p <- sapply(top6, function(p) kruskal.test(phy_df[[p]], phy_df$Disease_Status)$p.value)
phy_adj_p <- p.adjust(phy_raw_p, method = "BH")

cat("\nPhylum KW(BH) results:\n")
for (p in top6) cat(sprintf("  %-25s  raw p=%.4e   p.adj=%.4e\n", p, phy_raw_p[p], phy_adj_p[p]))

cat("\nFirmicutes/Bacteroidota ratio (median per group):\n")
fb_ratio <- phy_df$Firmicutes / phy_df$Bacteroidota
fb_med <- tapply(fb_ratio, phy_df$Disease_Status, function(x) median(x, na.rm = TRUE))
# Bootstrap 95% CI of the group median (Reviewer 2, point 17). Samples with zero
# Bacteroidota give F/B = Inf; they are kept (they rank highest), as in the median.
set.seed(2024)
fb_ci <- t(sapply(levels(phy_df$Disease_Status), function(g) {
  x <- fb_ratio[phy_df$Disease_Status == g]
  b <- replicate(2000, median(sample(x, replace = TRUE), na.rm = TRUE))
  c(n = length(x), zero_Bacteroidota = sum(phy_df$Bacteroidota[phy_df$Disease_Status == g] == 0),
    median = median(x, na.rm = TRUE), CI_lower = unname(quantile(b, 0.025, na.rm = TRUE)),
    CI_upper = unname(quantile(b, 0.975, na.rm = TRUE)))
}))
fb_ci <- data.frame(Group = rownames(fb_ci), fb_ci, row.names = NULL)
for (i in seq_len(nrow(fb_ci))) with(fb_ci[i, ], cat(sprintf("  %-10s median %.2f (95%% CI %.2f-%.2f); n = %d, zero Bacteroidota = %d\n",
                                                                Group, median, CI_lower, CI_upper, n, zero_Bacteroidota)))
write.csv(fb_ci, file.path(output_dir, "Figure5_FB_ratio_median_CI.csv"), row.names = FALSE)

phy_sig_df <- data.frame()
for (p in top6) {
  if (phy_adj_p[p] < 0.05) {
    ast <- ifelse(phy_adj_p[p] < 0.001, "***", ifelse(phy_adj_p[p] < 0.01, "**", "*"))
    pdat <- phy_b[phy_b$Phylum == p, ]
    phy_sig_df <- rbind(phy_sig_df,
                        data.frame(Phylum = p, label = ast, x_pos = max(pdat$Mean + pdat$SE)))
  }
}
if (nrow(phy_sig_df) > 0) phy_sig_df$Phylum <- factor(phy_sig_df$Phylum, levels = phy_order)

panel_B <- ggplot(phy_b, aes(x = Mean, y = Phylum, fill = Disease_Status)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           color = "black", linewidth = 0.2, alpha = 0.85) +
  geom_errorbar(aes(xmin = pmax(0, Mean - SE), xmax = Mean + SE),
                position = position_dodge(width = 0.8), width = 0.4,
                color = "#333333", linewidth = 0.3) +
  scale_fill_manual(values = group_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.20)), breaks = seq(0, 80, 20)) +
  scale_y_discrete(labels = phy_y_clean) +
  labs(x = "Mean Relative Abundance (%)", y = NULL, tag = "B") +
  base_theme +
  theme(
    axis.text.y        = element_text(size = 7, face = "italic", color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "#E5E5E5", linewidth = 0.3),
    legend.position    = "none"
  )

if (nrow(phy_sig_df) > 0) {
  panel_B <- panel_B +
    geom_text(data = phy_sig_df,
              aes(x = x_pos, y = Phylum, label = label),
              inherit.aes = FALSE,
              hjust = -0.3, size = 2.8, fontface = "bold", nudge_x = 0.5)
}

# ============================================================================
# 9. PANEL C — GENUS STACKED
# ============================================================================

cat("Building Panel C...\n")

gen_mean <- gen_df %>%
  group_by(Disease_Status) %>%
  summarise(across(all_of(gen_cols), mean), .groups = "drop") %>%
  pivot_longer(cols = all_of(gen_cols), names_to = "Genus", values_to = "Mean_Abundance") %>%
  mutate(Genus = factor(Genus, levels = gen_cols))

gen_stacked_colors <- setNames(genera_palette[1:15], top15)
gen_stacked_colors["Other"] <- "#CCCCCC"

panel_C <- ggplot(gen_mean, aes(x = Disease_Status, y = Mean_Abundance, fill = Genus)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2, width = 0.65) +
  scale_fill_manual(values = gen_stacked_colors, labels = gen_clean_labels, drop = FALSE) +
  scale_x_discrete(labels = xlab_gen) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)), breaks = seq(0, 100, 25)) +
  labs(x = NULL, y = "Relative Abundance (%)", fill = "Genus", tag = "C") +
  base_theme +
  theme(
    axis.text.x        = element_text(size = 7, face = "bold", color = "black",
                                      angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#E5E5E5", linewidth = 0.3),
    legend.position    = "right",
    legend.title       = element_text(size = 7, face = "bold"),
    legend.text        = element_text(size = 5.5, face = "italic"),
    legend.key.size    = unit(0.22, "cm"),
    legend.key.height  = unit(0.22, "cm"),
    legend.spacing.y   = unit(0.02, "cm"),
    legend.spacing.x   = unit(0.05, "cm")
  ) +
  guides(fill = guide_legend(ncol = 1, reverse = TRUE,
                             keywidth  = unit(0.25, "cm"),
                             keyheight = unit(0.22, "cm")))

# ============================================================================
# 10. PANEL D — GENUS BAR CHART
# ============================================================================

cat("Building Panel D...\n")

gen_b <- gen_df[, c(top15, "Disease_Status")] %>%
  pivot_longer(cols = all_of(top15), names_to = "Genus", values_to = "Abundance") %>%
  group_by(Disease_Status, Genus) %>%
  summarise(Mean = mean(Abundance), SE = sd(Abundance) / sqrt(n()), .groups = "drop")

gen_order <- gen_b %>% group_by(Genus) %>%
  summarise(Overall = mean(Mean), .groups = "drop") %>%
  arrange(Overall) %>% pull(Genus)
gen_b$Genus <- factor(gen_b$Genus, levels = gen_order)
gen_y_clean <- setNames(sapply(gen_order, make_clean_label), gen_order)

gen_raw_p <- sapply(top15, function(g) kruskal.test(gen_df[[g]], gen_df$Disease_Status)$p.value)
gen_adj_p <- p.adjust(gen_raw_p, method = "BH")

cat("\nGenus KW(BH) results:\n")
for (g in top15) cat(sprintf("  %-35s  raw p=%.4e   p.adj=%.4e\n",
                              make_clean_label(g), gen_raw_p[g], gen_adj_p[g]))

gen_sig_df <- data.frame()
for (g in top15) {
  if (gen_adj_p[g] < 0.05) {
    ast <- ifelse(gen_adj_p[g] < 0.001, "***", ifelse(gen_adj_p[g] < 0.01, "**", "*"))
    gdat <- gen_b[gen_b$Genus == g, ]
    gen_sig_df <- rbind(gen_sig_df,
                        data.frame(Genus = g, label = ast, x_pos = max(gdat$Mean + gdat$SE)))
  }
}
if (nrow(gen_sig_df) > 0) gen_sig_df$Genus <- factor(gen_sig_df$Genus, levels = gen_order)

panel_D <- ggplot(gen_b, aes(x = Mean, y = Genus, fill = Disease_Status)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           color = "black", linewidth = 0.2, alpha = 0.85) +
  geom_errorbar(aes(xmin = pmax(0, Mean - SE), xmax = Mean + SE),
                position = position_dodge(width = 0.8), width = 0.4,
                color = "#333333", linewidth = 0.3) +
  scale_fill_manual(values = group_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)), breaks = seq(0, 30, 10)) +
  scale_y_discrete(labels = gen_y_clean) +
  labs(x = "Mean Relative Abundance (%)", y = NULL, tag = "D") +
  base_theme +
  theme(
    axis.text.y        = element_text(size = 7, face = "italic", color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "#E5E5E5", linewidth = 0.3),
    legend.position    = "none",
    plot.margin        = margin(5, 5, 5, -10)   # Negative left margin shifts Panel D left
  )

if (nrow(gen_sig_df) > 0) {
  panel_D <- panel_D +
    geom_text(data = gen_sig_df,
              aes(x = x_pos, y = Genus, label = label),
              inherit.aes = FALSE,
              hjust = -0.3, size = 2.8, fontface = "bold", nudge_x = 0.3)
}

# ============================================================================
# 11. SHARED DISEASE LEGEND
# ============================================================================

cat("Building shared disease legend...\n")

shared_legend_plot <- ggplot(phy_b, aes(x = Mean, y = Phylum, fill = Disease_Status)) +
  geom_col() +
  scale_fill_manual(
    values = group_colors,
    labels = paste0(names(group_colors), " (n=", n_phy[names(group_colors)], ")")
  ) +
  labs(fill = "Disease Status   ") +
  theme_minimal() +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(size = 8,   face = "bold", margin = margin(r = 8)),
    legend.text        = element_text(size = 7.5, face = "bold", margin = margin(l = 2, r = 10)),
    legend.key.size    = unit(0.4, "cm"),
    legend.spacing.x   = unit(0.3, "cm"),
    legend.background  = element_rect(fill = "white", color = NA)
  ) +
  guides(fill = guide_legend(nrow = 1))

shared_disease_legend <- cowplot::get_legend(shared_legend_plot)

# ============================================================================
# 12. ASSEMBLE FIGURE
# ============================================================================

cat("Assembling Figure 5...\n")

# Strategy: Align Panel B and D vertically FIRST so their y-axis label widths
# match (forcing identical bar start positions). Same for A and C.
# This produces 4 panels with matching plot areas.

# Step 1: Vertically align column 1 (A above C) - forces y-axis label widths to match
left_col  <- cowplot::align_plots(panel_A, panel_C, align = "v", axis = "lr")

# Step 2: Vertically align column 2 (B above D) - same logic
right_col <- cowplot::align_plots(panel_B, panel_D, align = "v", axis = "lr")

# Step 3: Build rows from the now-aligned panels
row_phy <- cowplot::plot_grid(left_col[[1]], right_col[[1]],
                              ncol = 2, rel_widths = c(0.48, 0.52))
row_gen <- cowplot::plot_grid(left_col[[2]], right_col[[2]],
                              ncol = 2, rel_widths = c(0.48, 0.52))

# Spacer for visual separation between rows
spacer <- ggplot() + theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA))

figure5 <- cowplot::plot_grid(row_phy, spacer, row_gen, shared_disease_legend,
                              ncol = 1, rel_heights = c(0.42, 0.025, 0.50, 0.055))

# ============================================================================
# 13. EXPORT
# ============================================================================

cat("Saving Figure 5...\n")

ggsave("Figure5_TaxonomicComposition_PhylumGenus.pdf",
       plot = figure5, width = 174, height = 220, units = "mm")

ggsave("Figure5_TaxonomicComposition_PhylumGenus.png",
       plot = figure5, width = 174, height = 220, units = "mm",
       dpi = 600, bg = "white")

ggsave("Figure5_TaxonomicComposition_PhylumGenus.svg",
       plot = figure5, width = 174, height = 220, units = "mm",
       device = svglite)

cat("\nFigure 5 saved (Gut Microbes format):\n")
cat("  - Figure5_TaxonomicComposition_PhylumGenus.pdf  (vector, submission-ready)\n")
cat("  - Figure5_TaxonomicComposition_PhylumGenus.png  (600 dpi raster)\n")
cat("  - Figure5_TaxonomicComposition_PhylumGenus.svg  (editable vector)\n")
cat("Location:", getwd(), "\n\n")

# ============================================================================
# 14. STATISTICAL SUMMARY
# ============================================================================

cat("=", rep("=", 70), "\n", sep = "")
cat("FIGURE 5 — STATISTICAL SUMMARY\n")
cat("=", rep("=", 70), "\n", sep = "")

cat("\n--- PHYLUM (Kruskal-Wallis, BH-corrected) ---\n")
print(data.frame(
  Phylum = top6,
  Raw_p  = round(phy_raw_p, 6),
  Adj_p  = round(phy_adj_p, 6),
  Sig    = ifelse(phy_adj_p < 0.001, "***",
           ifelse(phy_adj_p < 0.01,  "**",
           ifelse(phy_adj_p < 0.05,  "*",   "ns")))
), row.names = FALSE)

cat("\n--- GENUS (Kruskal-Wallis, BH-corrected) ---\n")
print(data.frame(
  Genus = sapply(top15, make_clean_label),
  Raw_p = round(gen_raw_p, 6),
  Adj_p = round(gen_adj_p, 6),
  Sig   = ifelse(gen_adj_p < 0.001, "***",
          ifelse(gen_adj_p < 0.01,  "**",
          ifelse(gen_adj_p < 0.05,  "*",   "ns")))
), row.names = FALSE)

cat("\n--- F/B RATIO (median, bootstrap 95% CI) ---\n")
print(transform(fb_ci, median = round(median, 2), CI_lower = round(CI_lower, 2), CI_upper = round(CI_upper, 2)),
      row.names = FALSE)

cat("\n--- 'OTHER' ABUNDANCE (mean % per group) ---\n")
cat("Phylum 'Other':\n"); print(round(phy_other_means, 3))
cat("Genus 'Other':\n");  print(round(gen_other_means, 3))

cat("\nDone.\n")
