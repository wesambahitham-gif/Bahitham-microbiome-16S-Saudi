################################################################################
# SUPPLEMENTARY FIGURE S3: VARIANCE DECOMPOSITION & DF-ADJUSTED EFFECT SIZES
# Paper 1 — 16S rRNA Gut Microbiome | Saudi Arabian Liver Disease Cohort
# ------------------------------------------------------------------------------
# Shows that the disease effect persists after adjustment for 11 confounders and
# compares disease with each covariate on a degrees-of-freedom-adjusted scale.
# Companion to main Figure 4.
#
# REVISION (Reviewer 2, point 10): raw R² of a 3-df disease term is not directly
# comparable with 1-df covariates. Panel C now compares terms by omega-squared
# (ω², corrected for df and chance) instead of disease-R²/confounder-R² ratios,
# and Panel A no longer stacks the strongest confounder's R² onto the disease bar.
#
# 3-Panel Layout:
#   A: Disease variance decomposition (independent vs shared with confounders)
#   B: R² reduction (Unadjusted -> Adjusted)
#   C: ω² of disease status vs each covariate (adjusted marginal model)
#
# Generates: Supplementary_Figure_S3_Variance_Decomposition.pdf/.png/.svg
#            Table_S3_Variance_Decomposition.xlsx
#
# DEPENDENCIES:
#   1. Figure4_Confounder_PERMANOVA.R -> PERMANOVA_results_for_supplementary.RData
#   2. Supplementary_Analyses_Revision.R (section PERMANOVA_ROBUSTNESS)
#      -> Revision_MicrobiolSpectrum/Rev03_PERMANOVA_Robustness.xlsx (ω², R²/df)
#
# JOURNAL SPECS — Gut Microbes (Taylor & Francis):
#   - Width: 174 mm (full page)
#   - Resolution: 600 dpi (PNG); vector PDF; SVG editable
#   - Theme matches Figures 2-8 and Supp Fig S1 for cross-figure consistency
################################################################################

# ============================================================================
# 1. CONFIGURATION
# ============================================================================

# Paths are RELATIVE: outputs go to this script's folder; inputs are found automatically
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
BASE_DIR   <- find_proj_dir()                                   # Desktop/16S_reprocess
OUTPUT_DIR <- getwd()

# RData written by Figure4_Confounder_PERMANOVA.R (searched for under 16S_liver_disease/)
rdata_hits <- list.files(dirname(getwd()),
                         pattern = "^PERMANOVA_results_for_supplementary\\.RData$",
                         recursive = TRUE, full.names = TRUE)
RDATA_FILE <- if (length(rdata_hits)) rdata_hits[which.max(file.mtime(rdata_hits))] else
  file.path(dirname(getwd()), "Figure_4", "PERMANOVA_results_for_supplementary.RData")

# omega^2 and R^2/df for every term of the adjusted marginal model (same 11 covariates)
REV03_FILE <- file.path(dirname(getwd()), "Supp_Analyses", "Rev03_PERMANOVA_Robustness.xlsx")   # from Supplementary_Analyses_Revision.R
cat("RData :", RDATA_FILE, "\nRev03 :", REV03_FILE, "\n")

# ============================================================================
# 2. PACKAGES
# ============================================================================

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "scales", "patchwork",
                   "openxlsx", "RColorBrewer", "svglite", "readxl")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cran.r-project.org")
  library(pkg, character.only = TRUE)
}

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ============================================================================
# 3. DEFENSIVE FILE-EXISTENCE CHECK
# ============================================================================

if (!file.exists(REV03_FILE)) {
  stop("Missing input file:\n  REV03: ", REV03_FILE,
       "\n\nRun Supplementary_Analyses_Revision.R (RUN$PERMANOVA_ROBUSTNESS = TRUE) first.")
}
if (!file.exists(RDATA_FILE)) {
  stop("Missing input file:\n  RDATA: ", RDATA_FILE,
       "\n\nRun Figure4_Confounder_PERMANOVA.R first to generate ",
       "PERMANOVA_results_for_supplementary.RData, OR update RDATA_FILE path ",
       "if your Figure 4 outputs landed elsewhere.")
}
cat("Input files located.\n\n")
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate; rename <- dplyr::rename

# ============================================================================
# 4. LOAD PRECOMPUTED RESULTS
# ============================================================================

load(RDATA_FILE)
cat("Loaded precomputed PERMANOVA results.\n")
cat(sprintf("  Samples: %d | Metrics: %d | Confounders: %d\n",
            length(common_samples), length(dist_sub), length(confounders)))

# ============================================================================
# 4. THEME & COLORS (matching Figure 7)
# ============================================================================

metric_colors <- c(
  "Weighted UniFrac"   = "#4DAF4A",
  "Unweighted UniFrac" = "#984EA3",
  "Bray-Curtis"        = "#377EB8"
)

metric_order <- c("Weighted UniFrac", "Unweighted UniFrac", "Bray-Curtis")

metric_labels <- c("Weighted UniFrac"   = "Weighted\nUniFrac",
                    "Unweighted UniFrac" = "Unweighted\nUniFrac",
                    "Bray-Curtis"        = "Bray-Curtis")

theme_pub <- theme_minimal(base_size = 8) +
  theme(
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "#E5E5E5", linewidth = 0.3),
    panel.border       = element_blank(),
    axis.line          = element_line(color = "black", linewidth = 0.4),
    axis.ticks         = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length  = unit(0.08, "cm"),
    axis.text          = element_text(size = 7, color = "black"),
    axis.title         = element_text(size = 8.5, face = "bold"),
    legend.text        = element_text(size = 7),
    legend.title       = element_text(size = 7.5, face = "bold"),
    legend.key.size    = unit(0.35, "cm"),
    strip.text         = element_text(size = 8, face = "bold"),
    strip.background   = element_blank(),
    plot.title         = element_text(face = "bold", size = 9, hjust = 0),
    plot.subtitle      = element_text(size = 7.5, color = "grey40"),
    plot.tag           = element_text(size = 10, face = "bold"),
    plot.margin        = margin(5, 5, 5, 5)
  )

# ============================================================================
# 5. COMPUTE VARIANCE DECOMPOSITION
# ============================================================================

cat("\n--- Computing variance decomposition ---\n")

decomp <- merged %>%
  filter(Variable == "Disease Status") %>%
  select(Metric, R2_unadjusted, R2_adjusted) %>%
  mutate(
    R2_shared = R2_unadjusted - R2_adjusted,
    R2_independent = R2_adjusted,
    Reduction_pct = R2_shared / R2_unadjusted * 100,
    Retained_pct = R2_independent / R2_unadjusted * 100,
    # Percentage versions
    R2_unadj_pct = R2_unadjusted * 100,
    R2_adj_pct = R2_adjusted * 100,
    R2_shared_pct = R2_shared * 100,
    R2_indep_pct = R2_independent * 100
  )

# Max single confounder per metric (NOT sum)
conf_stats <- merged %>%
  filter(Type == "Confounder", !is.na(R2_adjusted)) %>%
  group_by(Metric) %>%
  summarise(
    Total_confounder_R2 = sum(R2_adjusted),
    Max_confounder_name = Variable[which.max(R2_adjusted)],
    Max_confounder_R2 = max(R2_adjusted),
    .groups = "drop"
  )

decomp <- decomp %>% left_join(conf_stats, by = "Metric")

decomp <- decomp %>%
  mutate(
    # Ratio vs strongest single confounder (meaningful comparison)
    Disease_vs_Max_Ratio = R2_independent / Max_confounder_R2,
    Max_conf_pct = Max_confounder_R2 * 100,
    Total_conf_pct = Total_confounder_R2 * 100
  )

cat("Variance Decomposition:\n")
for (i in seq_len(nrow(decomp))) {
  cat(sprintf("  %s: Unadj=%.2f%%, Adj=%.2f%% (-%.1f%%), Max confounder=%s (%.2f%%), Disease:Max=%.1fx\n",
              decomp$Metric[i], decomp$R2_unadj_pct[i], decomp$R2_adj_pct[i],
              decomp$Reduction_pct[i], decomp$Max_confounder_name[i],
              decomp$Max_conf_pct[i], decomp$Disease_vs_Max_Ratio[i]))
}

# ---------------------------------------------------------------------------
# df-adjusted effect sizes (ω², R²/df) from the adjusted marginal model
# ---------------------------------------------------------------------------
term_labels <- c(Status = "Disease Status", Gender = "Sex", Age_numeric = "Age", BMI_numeric = "BMI",
                 TG_numeric = "TG", LDL_numeric = "LDL", HDL_numeric = "HDL", PPI_Use = "PPI Use",
                 Metformin = "Metformin", Statin_Use = "Statin Use", Rifaximin_Use = "Rifaximin",
                 Lactulose_Use = "Lactulose")
omega <- readxl::read_excel(REV03_FILE, sheet = "A_Effect_sizes_all_terms") %>%
  filter(Model == "Adjusted (marginal)", Term %in% names(term_labels), Metric %in% metric_order) %>%
  mutate(Variable = unname(term_labels[Term]),
         Type = ifelse(Term == "Status", "Disease", "Confounder")) %>%
  select(Metric, Variable, Type, Df, R2, R2_per_df, Omega2, F, P)
if (!setequal(unique(omega$Metric), metric_order))
  stop("Rev03 workbook does not contain all metrics: ", paste(metric_order, collapse = ", "))

omega_summary <- omega %>% group_by(Metric) %>%
  summarise(Disease_omega2   = Omega2[Type == "Disease"],
            Disease_R2_per_df = R2_per_df[Type == "Disease"],
            Max_cov_omega2_name = Variable[Type == "Confounder"][which.max(Omega2[Type == "Confounder"])],
            Max_cov_omega2   = max(Omega2[Type == "Confounder"]),
            Max_cov_R2df_name = Variable[Type == "Confounder"][which.max(R2_per_df[Type == "Confounder"])],
            Max_cov_R2_per_df = max(R2_per_df[Type == "Confounder"]),
            Disease_rank_omega2 = rank(-Omega2)[Type == "Disease"],
            Disease_rank_R2_per_df = rank(-R2_per_df)[Type == "Disease"],
            .groups = "drop") %>%
  mutate(Omega2_ratio = Disease_omega2 / Max_cov_omega2,
         R2_per_df_ratio = Disease_R2_per_df / Max_cov_R2_per_df)
decomp <- decomp %>% left_join(omega_summary, by = "Metric")

cat("\ndf-adjusted comparison (adjusted marginal model):\n")
for (i in seq_len(nrow(omega_summary))) with(omega_summary[i, ], cat(sprintf(
  "  %s: disease ω²=%.4f (rank %d) vs %s %.4f -> %.1fx | R²/df %.4f vs %s %.4f -> %.2fx\n",
  Metric, Disease_omega2, Disease_rank_omega2, Max_cov_omega2_name, Max_cov_omega2, Omega2_ratio,
  Disease_R2_per_df, Max_cov_R2df_name, Max_cov_R2_per_df, R2_per_df_ratio)))

# Per-confounder ratios (raw R²; kept in the Excel output for reference only)
confounder_ratios <- merged %>%
  filter(Type == "Confounder", !is.na(R2_adjusted)) %>%
  group_by(Metric) %>%
  mutate(
    Disease_R2 = merged$R2_adjusted[merged$Variable == "Disease Status" &
                                      merged$Metric == Metric[1]],
    Ratio_vs_Disease = Disease_R2 / R2_adjusted
  ) %>%
  ungroup()

# ============================================================================
# 6. GROUPED VARIABLE ORDERING (matching Figure 7)
# ============================================================================

demo_clinical <- c("Sex", "Age", "BMI", "TG", "LDL", "HDL")
medications   <- c("PPI Use", "Metformin", "Statin Use", "Rifaximin",
                    "Lactulose")

mean_r2 <- merged %>%
  filter(Type == "Confounder") %>%
  group_by(Variable) %>%
  summarise(mean_r2 = mean(R2_adjusted, na.rm = TRUE), .groups = "drop")

demo_order <- mean_r2 %>%
  filter(Variable %in% demo_clinical) %>%
  arrange(desc(mean_r2)) %>%
  pull(Variable)

med_order <- mean_r2 %>%
  filter(Variable %in% medications) %>%
  arrange(desc(mean_r2)) %>%
  pull(Variable)

var_order_confounders <- c(demo_order, med_order)

# ============================================================================
# 7. PANEL A: Variance Decomposition (Disease vs Max Single Confounder)
# ============================================================================

stack_data <- decomp %>%
  select(Metric,
         `Independent Disease` = R2_indep_pct,
         `Shared (Disease + Confounders)` = R2_shared_pct) %>%
  pivot_longer(-Metric, names_to = "Component", values_to = "R2_pct") %>%
  mutate(
    Component = factor(Component, levels = c(
      "Shared (Disease + Confounders)",
      "Independent Disease"
    )),
    Metric = factor(Metric, levels = metric_order)
  )

panelA <- ggplot(stack_data, aes(x = Metric, y = R2_pct, fill = Component)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.6) +
  geom_text(aes(label = sprintf("%.2f%%", R2_pct)),
            position = position_stack(vjust = 0.5), size = 2.8, fontface = "bold") +
  scale_fill_manual(
    values = c(
      "Independent Disease" = "#D73027",
      "Shared (Disease + Confounders)" = "#FEE08B"
    ),
    name = "Variance Component"
  ) +
  scale_x_discrete(labels = metric_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Variance Explained (%)",
       title = "Disease variance: independent vs. shared with covariates") +
  theme_pub +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 8.5, hjust = 0))

# ============================================================================
# 8. PANEL B: R-squared Reduction Waterfall (percentages)
# ============================================================================

waterfall_data <- decomp %>%
  select(Metric, Unadjusted = R2_unadj_pct, Adjusted = R2_adj_pct,
         Reduction_pct) %>%
  pivot_longer(cols = c(Unadjusted, Adjusted),
               names_to = "Model", values_to = "R2_pct") %>%
  mutate(
    Model = factor(Model, levels = c("Unadjusted", "Adjusted")),
    Metric = factor(Metric, levels = metric_order)
  )

reduction_labels <- decomp %>%
  mutate(
    y_pos = (R2_unadj_pct + R2_adj_pct) / 2,
    label = sprintf("-%.1f%%", Reduction_pct),
    Metric = factor(Metric, levels = metric_order)
  )

panelB <- ggplot(waterfall_data, aes(x = Model, y = R2_pct, fill = Model)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.5) +
  geom_text(aes(label = sprintf("%.2f%%", R2_pct)), vjust = -0.5, 
            size = 2.8, fontface = "bold") +
  geom_segment(data = reduction_labels,
               aes(x = 1.3, xend = 1.7, y = R2_unadj_pct, yend = R2_adj_pct),
               inherit.aes = FALSE, color = "red", linewidth = 0.6,
               arrow = arrow(length = unit(0.12, "cm"), type = "closed")) +
  geom_label(data = reduction_labels,
             aes(x = 1.5, y = y_pos, label = label),
             inherit.aes = FALSE, size = 2.8, fontface = "bold",
             color = "red", fill = "white", label.size = 0.25) +
  scale_fill_manual(values = c("Unadjusted" = "#2C3E50", "Adjusted" = "#E74C3C")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  facet_wrap(~ Metric, ncol = 3) +
  labs(x = NULL, y = "Disease Variance Explained (%)",
       title = expression(bold(paste("Disease ", R^2, " before and after adjustment for 11 covariates")))) +
  theme_pub +
  theme(legend.position = "none",
        strip.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(face = "bold", size = 8.5, hjust = 0),
        plot.subtitle = element_text(size = 7, color = "grey40"))

# ============================================================================
# 9. PANEL C: df-adjusted effect size (ω²) of every term in the adjusted model
#    Replaces the R² ratio panel (Reviewer 2, point 10)
# ============================================================================

omega_plot <- omega %>%
  mutate(Variable = factor(Variable, levels = rev(c("Disease Status", var_order_confounders))),
         Metric   = factor(Metric, levels = metric_order),
         Fill     = ifelse(Type == "Disease", "Disease status (3 df)", "Covariate (1 df)"))
facet_lab <- setNames(sprintf("%s\nDisease %.1fx next\n(%s)",
                              gsub("\n", " ", metric_labels[omega_summary$Metric]),
                              omega_summary$Omega2_ratio, omega_summary$Max_cov_omega2_name),
                      omega_summary$Metric)

panelC <- ggplot(omega_plot, aes(x = Variable, y = Omega2, fill = Fill)) +
  geom_col(color = "black", linewidth = 0.25, width = 0.7) +
  geom_vline(xintercept = length(var_order_confounders) + 0.5, linetype = "dotted",
             color = "grey50", linewidth = 0.4) +
  facet_wrap(~ Metric, nrow = 1, scales = "free_x", labeller = labeller(Metric = facet_lab)) +
  coord_flip() +
  scale_fill_manual(values = c("Disease status (3 df)" = "#D73027", "Covariate (1 df)" = "grey70"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                     breaks = function(l) { b <- pretty(c(0, l[2]), n = 2); b[b <= 0.85 * l[2]] },
                     labels = function(x) sprintf("%.3f", x)) +
  labs(x = NULL,
       y = expression(paste(omega^2, " (adjusted marginal model)")),
       title = expression(bold(paste("Degrees-of-freedom-adjusted effect size (", omega^2, ") of each term")))) +
  theme_pub +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_line(color = "#E5E5E5", linewidth = 0.3),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(size = 7, face = "bold", lineheight = 0.9),
        panel.spacing.x = unit(0.8, "cm"),
        plot.title = element_text(face = "bold", size = 8.5, hjust = 0))

# ============================================================================
# 10. ASSEMBLE 3-PANEL SUPPLEMENTARY FIGURE S3
# ============================================================================

cat("\n--- Assembling Supplementary Figure S3 (3 panels) ---\n")

suppfig_S3 <- panelA / panelB / panelC +
  plot_layout(heights = c(1, 0.9, 1.1)) +
  plot_annotation(   # no figure title/subtitle: details go in the figure legend
    tag_levels = "A",
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  )

# ============================================================================
# EXPORT — Gut Microbes specifications
# Width: 174 mm (full page) | Height: 280 mm | 600 dpi PNG | vector PDF/SVG
# 3 panels stacked vertically (each gets full page width for readability)
# ============================================================================

cat("Saving Supplementary Figure S3...\n")

ggsave(file.path(OUTPUT_DIR, "Supplementary_Figure_S3_Variance_Decomposition.pdf"),
       suppfig_S3, width = 174, height = 280, units = "mm")

ggsave(file.path(OUTPUT_DIR, "Supplementary_Figure_S3_Variance_Decomposition.png"),
       suppfig_S3, width = 174, height = 280, units = "mm",
       dpi = 600, bg = "white")

ggsave(file.path(OUTPUT_DIR, "Supplementary_Figure_S3_Variance_Decomposition.svg"),
       suppfig_S3, width = 174, height = 280, units = "mm",
       device = svglite)

cat("\nSupplementary Figure S3 saved:\n")
cat("  - Supplementary_Figure_S3_Variance_Decomposition.pdf  (vector, submission-ready)\n")
cat("  - Supplementary_Figure_S3_Variance_Decomposition.png  (600 dpi raster)\n")
cat("  - Supplementary_Figure_S3_Variance_Decomposition.svg  (editable vector)\n")
cat("Location:", OUTPUT_DIR, "\n")

# ============================================================================
# 11. SUPPLEMENTARY EXCEL TABLE
# ============================================================================

cat("\n--- Saving Supplementary Excel ---\n")

wb2 <- createWorkbook()

# Sheet 1: Variance decomposition summary
addWorksheet(wb2, "Variance_Decomposition")
decomp_export <- decomp %>%
  mutate(Metric = factor(Metric, levels = metric_order)) %>%
  arrange(Metric) %>%
  select(Metric,
         R2_unadjusted, R2_adjusted, R2_shared, R2_independent,
         Reduction_pct, Retained_pct,
         Disease_omega2, Max_cov_omega2_name, Max_cov_omega2, Omega2_ratio, Disease_rank_omega2,
         Disease_R2_per_df, Max_cov_R2df_name, Max_cov_R2_per_df, R2_per_df_ratio, Disease_rank_R2_per_df,
         Max_confounder_name, Max_confounder_R2,
         Raw_R2_ratio_not_df_adjusted = Disease_vs_Max_Ratio,
         Total_confounder_R2)
writeData(wb2, "Variance_Decomposition", decomp_export)

# Sheet 2: df-adjusted effect sizes for every term
addWorksheet(wb2, "Omega2_All_Terms")
writeData(wb2, "Omega2_All_Terms", omega %>%
            mutate(Metric = factor(Metric, levels = metric_order)) %>% arrange(Metric, desc(Omega2)))

# Sheet 3: Per-confounder raw R² ratios (reference only; not df-adjusted)
addWorksheet(wb2, "Confounder_Ratios")
ratio_export <- confounder_ratios %>%
  mutate(
    Group = ifelse(Variable %in% demo_clinical, "Demographics/Clinical", "Medications"),
    Variable = factor(Variable, levels = var_order_confounders),
    Metric = factor(Metric, levels = metric_order)
  ) %>%
  arrange(Metric, Group, Variable) %>%
  select(Metric, Group, Variable, R2_adjusted, FDR_adjusted,
         Significance_adjusted = Sig_adjusted,
         Disease_R2, Ratio_vs_Disease)
writeData(wb2, "Confounder_Ratios", ratio_export)

# Sheet 3: Interpretation guide
interp_omega <- data.frame(
  Component = c("Omega2", "R2_per_df", "Omega2_ratio", "R2_per_df_ratio"),
  Description = c("(SS_term - df_term x MS_residual) / (SS_total + MS_residual); truncated at 0",
                  "R2 / df of the term",
                  "Disease omega2 / largest covariate omega2",
                  "Disease R2 per df / largest covariate R2 per df"),
  Interpretation = c("Effect size corrected for degrees of freedom and chance; comparable across terms",
                     "Average variance per degree of freedom",
                     ">1: disease has the largest df-adjusted effect",
                     "~1: per df, disease is comparable to the strongest single covariate"),
  Confounders = "Same 11 confounders")
interp <- data.frame(
  Component = c("R2_unadjusted", "R2_adjusted (independent)", "R2_shared",
                "Max_confounder_R2", "Raw_R2_ratio_not_df_adjusted", "Total_confounder_R2"),
  Description = c(
    "Disease R2 from univariate PERMANOVA (disease only)",
    "Disease R2 from marginal PERMANOVA (adjusted for all 11 confounders)",
    "Variance shared between disease and confounders (unadjusted - adjusted)",
    "R2 of the single strongest confounder (marginal, adjusted for other variables)",
    "Disease adjusted R2 / max single confounder R2. NOT df-adjusted (disease has 3 df, covariates 1 df); use Omega2_ratio instead",
    "Sum of all 11 confounder marginal R2 values (provided for reference; NOT used for ratios)"
  ),
  Interpretation = c(
    "Total variance attributable to disease (may include confounder effects)",
    "True independent disease effect after removing confounder contributions",
    "Variance confounded by clinical covariates",
    "The most influential individual confounder after adjustment",
    "Reference only; inflated by the extra degrees of freedom of the disease term",
    "Note: summing marginal R2 is not directly comparable to single-variable R2"
  ),
  Confounders = c(
    "N/A",
    "Sex, Age, BMI, TG, LDL, HDL, PPI, Metformin, Statin, Rifaximin, Lactulose",
    "Same 11 confounders",
    "Same 11 confounders",
    "Same 11 confounders",
    "Same 11 confounders"
  )
)
addWorksheet(wb2, "Interpretation_Guide")
writeData(wb2, "Interpretation_Guide", rbind(interp, interp_omega))

saveWorkbook(wb2, file.path(OUTPUT_DIR, "Table_S3_Variance_Decomposition.xlsx"),
             overwrite = TRUE)
cat("  Supplementary Excel saved (4 sheets).\n")

# ============================================================================
# 12. CONSOLE SUMMARY
# ============================================================================

cat("\n")
cat("=", rep("=", 70), "\n", sep = "")
cat("SUPPLEMENTARY FIGURE S3 - VARIANCE DECOMPOSITION COMPLETE\n")
cat("=", rep("=", 70), "\n", sep = "")

cat("\n11 Confounders:\n")
cat("  Demographics/Clinical: Sex, Age, BMI, TG, LDL, HDL\n")
cat("  Medications: PPI, Metformin, Statin, Rifaximin, Lactulose\n")

cat("\nVariance Decomposition:\n")
for (i in seq_len(nrow(decomp))) {
  cat(sprintf(
    "  %s:\n    Unadjusted  = %.2f%%\n    Adjusted    = %.2f%% (reduction: %.1f%%)\n    Shared      = %.2f%%\n    omega2      = disease %.4f vs %s %.4f (%.1fx)\n    R2 per df   = disease %.4f vs %s %.4f (%.2fx)\n",
    decomp$Metric[i], decomp$R2_unadj_pct[i], decomp$R2_adj_pct[i],
    decomp$Reduction_pct[i], decomp$R2_shared_pct[i],
    decomp$Disease_omega2[i], decomp$Max_cov_omega2_name[i], decomp$Max_cov_omega2[i], decomp$Omega2_ratio[i],
    decomp$Disease_R2_per_df[i], decomp$Max_cov_R2df_name[i], decomp$Max_cov_R2_per_df[i], decomp$R2_per_df_ratio[i]
  ))
}

cat(sprintf("\nFiles generated in: %s\n", OUTPUT_DIR))
cat("  Supplementary_Figure_S3_Variance_Decomposition.pdf / .png / .svg (3 panels)\n")
cat("  Table_S3_Variance_Decomposition.xlsx (4 sheets)\n")
cat("=", rep("=", 70), "\n", sep = "")
