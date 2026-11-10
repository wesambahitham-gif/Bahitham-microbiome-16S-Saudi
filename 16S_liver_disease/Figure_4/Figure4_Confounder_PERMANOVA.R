################################################################################
# FIGURE 4: CONFOUNDER-ADJUSTED PERMANOVA
# Paper 1 — 16S rRNA Gut Microbiome | Saudi Arabian Liver Disease Cohort
# ------------------------------------------------------------------------------
# Demonstrates that disease status remains a significant explanatory variable
# for microbiome composition after adjusting for 11 potential confounders
# (demographics, clinical chemistry, medications).
#
# Layout:
#   Panel A: Disease-explained R² across 3 distance metrics (Unadjusted vs Adjusted)
#   Panel B: Individual variable R² per metric, ranked
#
# Generates: Figure4_Confounder_PERMANOVA.pdf/.png/.svg
#            PERMANOVA_Confounder_Results.xlsx (5 sheets)
#            Manuscript_PERMANOVA_Text.txt
#            PERMANOVA_results_for_supplementary.RData
#
# Input: Pairwise distance files (long format) + metadata
# Pipeline: QIIME2 → exported pairwise distances → this script
#
# JOURNAL SPECS — Gut Microbes (Taylor & Francis):
#   - Width: 174 mm (double-column / full page width)
#   - Height: 130 mm (proportional for 2-row layout)
#   - Resolution: 600 dpi (PNG); vector PDF
################################################################################

# ============================================================================
# 1. CONFIGURATION - UPDATE PATHS HERE
# ============================================================================

# Paths are RELATIVE to this script's folder (outputs are saved here)
#   script : 16S_liver_disease/Figure_4/
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

DATA_DIR <- file.path("..", "..", "16S_reprocess", "qiime2_analysis_v2")

DIST_FILES <- list(
  "Bray-Curtis"         = file.path(DATA_DIR, "Beta_diversity/Bray-curtis/Bray-Curtis_raw_data.tsv"),
  "Unweighted UniFrac"  = file.path(DATA_DIR, "Beta_diversity/UnWeighted_UniFrac/Unweighted.tsv"),
  "Weighted UniFrac"    = file.path(DATA_DIR, "Beta_diversity/Weighted_UniFrac/Weighted_rowdata.tsv")
)

METADATA_FILE <- file.path("..", "..", "16S_reprocess", "sample_metadata.tsv")

OUTPUT_DIR <- "."

# Reproducibility
set.seed(42)
N_PERM <- 999

# ============================================================================
# DEFENSIVE FILE-EXISTENCE CHECK (fail early with clear message)
# ============================================================================

required_files <- c(unlist(DIST_FILES), METADATA_FILE)
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop("Missing input files:\n  ",
       paste(missing, collapse = "\n  "),
       "\n\nWorking directory: ", getwd(),
       "\nSave this script in 16S_liver_disease/Figure_4/ and click Source.")
}
cat("All input files located.\n\n")

# ============================================================================
# 2. PACKAGES
# ============================================================================

required_pkgs <- c("vegan", "dplyr", "tidyr", "ggplot2", "scales",
                   "patchwork", "openxlsx", "RColorBrewer", "grid", "gridExtra",
                   "svglite")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("vegan")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg, repos = "https://cran.r-project.org")
    }
  }
  library(pkg, character.only = TRUE)
}

# ============================================================================
# 3. COLOR SCHEME & THEME
# ============================================================================

group_colors <- c(
  "Healthy"   = "#4DAF4A",
  "MASLD"     = "#377EB8",
  "Cirrhosis" = "#984EA3",
  "HCC"       = "#E41A1C"
)

colors_model <- c("Unadjusted" = "#2C3E50", "Adjusted" = "#E74C3C")

colors_vars <- c(
  "Disease Status"     = "#E41A1C",
  "Sex"                = "#377EB8",
  "PPI Use"            = "#4DAF4A",
  "BMI"                = "#984EA3",
  "Age"                = "#FF7F00",
  "Metformin"          = "#A65628",
  "Rifaximin"          = "#F781BF",
  "Lactulose"          = "#999999",
  "Statin Use"         = "#66C2A5",
  "Antibiotics (3 mo)" = "#8DA0CB",
  "TG"                 = "#E78AC3",
  "LDL"                = "#FC8D62",
  "HDL"                = "#B3B3B3"
)

metric_colors <- c(
  "Bray-Curtis"        = "#377EB8",
  "Unweighted UniFrac" = "#984EA3",
  "Weighted UniFrac"   = "#4DAF4A"
)

# ============================================================================
# UNIFIED PUBLICATION THEME (matches Figure 3 — Gut Microbes specs)
# ============================================================================

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
    axis.text          = element_text(size = 7,   color = "black"),
    axis.title         = element_text(size = 8.5, face = "bold"),
    legend.text        = element_text(size = 8,   face = "bold"),
    legend.key.size    = unit(0.4, "cm"),
    strip.text         = element_text(size = 8.5, face = "bold"),
    strip.background   = element_blank(),
    plot.title         = element_text(size = 8.5, face = "bold", hjust = 0),
    plot.tag           = element_text(size = 10,  face = "bold"),
    plot.margin        = margin(5, 5, 5, 5)
  )

# ============================================================================
# 4. LOAD METADATA
# ============================================================================

cat("Loading metadata...\n")
meta <- read.delim(METADATA_FILE, sep = "\t", stringsAsFactors = FALSE,
                   comment.char = "", check.names = FALSE)
colnames(meta)[1] <- "SampleID"

# Rename NAFLD → MASLD
meta$Status <- gsub("NAFLD", "MASLD", meta$Status)
meta$Status <- factor(meta$Status, levels = c("Healthy", "MASLD", "Cirrhosis", "HCC"))

# Prepare confounder variables
meta$Gender <- factor(meta$Gender)
meta$Age_numeric <- as.numeric(meta$Age)
meta$BMI_numeric <- as.numeric(meta$BMI)
meta$PPI_Use <- factor(meta$PPI_Use)
meta$Metformin <- factor(meta$Metformin)
meta$Antibiotic_Use_Last3mo <- factor(meta$Antibiotic_Use_Last3mo)
meta$Rifaximin_Use <- factor(meta$Rifaximin_Use)
meta$Lactulose_Use <- factor(meta$Lactulose_Use)
meta$Statin_Use <- factor(meta$Statin_Use)
meta$TG_numeric <- as.numeric(meta$TG)
meta$LDL_numeric <- as.numeric(meta$LDL)
meta$HDL_numeric <- as.numeric(meta$HDL)

rownames(meta) <- meta$SampleID

# ---- Impute isolated missing confounder values (retain all 205) ----
# Continuous -> cohort median; binary -> cohort mode. Empty strings treated as missing.
.continuous_cols <- c("Age_numeric","BMI_numeric","TG_numeric","LDL_numeric","HDL_numeric")
.binary_cols     <- c("Gender","PPI_Use","Metformin","Rifaximin_Use","Lactulose_Use","Statin_Use")
.get_mode <- function(x){ x <- x[!is.na(x) & x != ""]; names(sort(table(x), decreasing = TRUE))[1] }
for (col in .continuous_cols) {
  miss <- is.na(meta[[col]])
  if (any(miss)) {
    med <- median(meta[[col]], na.rm = TRUE); meta[[col]][miss] <- med
    cat(sprintf("  Imputed %d missing %s -> median %.2f\n", sum(miss), col, med))
  }
}
for (col in .binary_cols) {
  x <- as.character(meta[[col]]); miss <- is.na(x) | x == ""
  if (any(miss)) {
    mv <- .get_mode(x); x[miss] <- mv; meta[[col]] <- factor(x)
    cat(sprintf("  Imputed %d missing %s -> mode '%s'\n", sum(miss), col, mv))
  }
}

cat(sprintf("  Metadata: %d samples, %d variables\n", nrow(meta), ncol(meta)))
cat(sprintf("  Groups: %s\n", paste(table(meta$Status), collapse = ", ")))

# ============================================================================
# 5. DEFINE CONFOUNDERS
# ============================================================================

confounders <- c(
  "Sex"                   = "Gender",
  "Age"                   = "Age_numeric",
  "BMI"                   = "BMI_numeric",
  "PPI Use"               = "PPI_Use",
  "Metformin"             = "Metformin",
  "Rifaximin"             = "Rifaximin_Use",
  "Lactulose"             = "Lactulose_Use",
  "Statin Use"            = "Statin_Use",
  "TG"                    = "TG_numeric",
  "LDL"                   = "LDL_numeric",
  "HDL"                   = "HDL_numeric"
)

# ============================================================================
# 6. FUNCTION: Convert pairwise long format → square distance matrix
# ============================================================================

pairwise_to_dist <- function(filepath) {
  cat(sprintf("  Reading %s...\n", basename(filepath)))
  df <- read.delim(filepath, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  
  # Columns: (index), SubjectID1, SubjectID2, Group1, Group2, Distance
  samples <- sort(unique(c(df$SubjectID1, df$SubjectID2)))
  n <- length(samples)
  
  mat <- matrix(0, nrow = n, ncol = n, dimnames = list(samples, samples))
  
  for (i in seq_len(nrow(df))) {
    s1 <- df$SubjectID1[i]
    s2 <- df$SubjectID2[i]
    d  <- df$Distance[i]
    mat[s1, s2] <- d
    mat[s2, s1] <- d
  }
  
  return(as.dist(mat))
}

# ============================================================================
# 7. LOAD DISTANCE MATRICES & ALIGN SAMPLES
# ============================================================================

cat("\nLoading distance matrices...\n")
dist_list <- lapply(DIST_FILES, pairwise_to_dist)

# Get common samples across all matrices + metadata
common_samples <- Reduce(intersect, lapply(dist_list, function(d) labels(d)))
common_samples <- intersect(common_samples, meta$SampleID)

# Remove samples with NA in any confounder
meta_sub <- meta[common_samples, ]
na_mask <- complete.cases(meta_sub[, unname(confounders)])
common_samples <- common_samples[na_mask]
meta_sub <- meta[common_samples, ]

# Drop factor levels and remove confounders with <2 levels
for (col in unname(confounders)) {
  if (is.factor(meta_sub[[col]])) {
    meta_sub[[col]] <- droplevels(meta_sub[[col]])
  }
}

dropped <- c()
for (nm in names(confounders)) {
  col <- confounders[nm]
  if (is.factor(meta_sub[[col]]) && nlevels(meta_sub[[col]]) < 2) {
    cat(sprintf("  ⚠ Dropping '%s' — only 1 level ('%s') in filtered data\n",
                nm, levels(meta_sub[[col]])))
    dropped <- c(dropped, nm)
  }
}
if (length(dropped) > 0) {
  confounders <- confounders[!names(confounders) %in% dropped]
  cat(sprintf("  Remaining confounders: %d\n", length(confounders)))
}

cat(sprintf("  Common samples (no NA): %d\n", length(common_samples)))
cat(sprintf("  Groups after filtering: %s\n",
            paste(names(table(meta_sub$Status)), table(meta_sub$Status), 
                  sep = "=", collapse = ", ")))

# Subset distance matrices to common samples
dist_sub <- lapply(dist_list, function(d) {
  mat <- as.matrix(d)
  as.dist(mat[common_samples, common_samples])
})

# ============================================================================
# 8. UNIVARIATE PERMANOVA (each variable independently)
# ============================================================================

cat("\n--- Univariate PERMANOVA ---\n")

univariate_results <- data.frame()

all_vars <- c("Disease Status" = "Status", confounders)

for (metric_name in names(dist_sub)) {
  cat(sprintf("\n  Metric: %s\n", metric_name))
  d <- dist_sub[[metric_name]]
  
  for (var_label in names(all_vars)) {
    var_col <- all_vars[[var_label]]
    
    # Build formula
    frm <- as.formula(paste("d ~", var_col))
    
    # Run adonis2
    res <- adonis2(frm, data = meta_sub, permutations = N_PERM, by = "terms")
    
    r2 <- res[var_col, "R2"]
    f_val <- res[var_col, "F"]
    p_val <- res[var_col, "Pr(>F)"]
    
    univariate_results <- rbind(univariate_results, data.frame(
      Metric = metric_name,
      Variable = var_label,
      Variable_Col = var_col,
      R2 = r2,
      F_statistic = f_val,
      P_value = p_val,
      Type = ifelse(var_label == "Disease Status", "Primary", "Confounder"),
      stringsAsFactors = FALSE
    ))
    
    cat(sprintf("    %s: R²=%.4f, F=%.2f, p=%.4f\n", var_label, r2, f_val, p_val))
  }
}

# FDR correction within each metric
univariate_results <- univariate_results %>%
  group_by(Metric) %>%
  mutate(FDR = p.adjust(P_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    Significance = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE        ~ "ns"
    )
  )

# ============================================================================
# 9. ADJUSTED PERMANOVA (marginal effects, Type III)
# ============================================================================

cat("\n--- Adjusted PERMANOVA (marginal/Type III) ---\n")

adjusted_results <- data.frame()

for (metric_name in names(dist_sub)) {
  cat(sprintf("\n  Metric: %s\n", metric_name))
  d <- dist_sub[[metric_name]]
  
  # Full model: Disease + all confounders
  rhs <- paste(c("Status", unname(confounders)), collapse = " + ")
  frm <- as.formula(paste("d ~", rhs))
  
  res <- adonis2(frm, data = meta_sub, permutations = N_PERM, by = "margin")
  
  # Extract each variable
  var_names <- rownames(res)
  var_names <- var_names[!var_names %in% c("Residual", "Total")]
  
  for (v in var_names) {
    # Map back to label
    label <- names(all_vars)[match(v, all_vars)]
    if (is.na(label)) label <- v
    
    adjusted_results <- rbind(adjusted_results, data.frame(
      Metric = metric_name,
      Variable = label,
      Variable_Col = v,
      R2_adjusted = res[v, "R2"],
      F_adjusted = res[v, "F"],
      P_adjusted = res[v, "Pr(>F)"],
      Type = ifelse(v == "Status", "Primary", "Confounder"),
      stringsAsFactors = FALSE
    ))
    
    cat(sprintf("    %s: R²=%.4f, F=%.2f, p=%.4f\n", 
                label, res[v, "R2"], res[v, "F"], res[v, "Pr(>F)"]))
  }
}

# FDR correction within each metric
adjusted_results <- adjusted_results %>%
  group_by(Metric) %>%
  mutate(FDR_adjusted = p.adjust(P_adjusted, method = "BH")) %>%
  ungroup() %>%
  mutate(
    Significance_adjusted = case_when(
      FDR_adjusted < 0.001 ~ "***",
      FDR_adjusted < 0.01  ~ "**",
      FDR_adjusted < 0.05  ~ "*",
      TRUE                 ~ "ns"
    )
  )

# ============================================================================
# 10. MERGE RESULTS: Unadjusted vs Adjusted
# ============================================================================

merged <- univariate_results %>%
  select(Metric, Variable, R2_unadjusted = R2, F_unadjusted = F_statistic,
         P_unadjusted = P_value, FDR_unadjusted = FDR, 
         Sig_unadjusted = Significance, Type) %>%
  left_join(
    adjusted_results %>%
      select(Metric, Variable, R2_adjusted, F_adjusted, P_adjusted,
             FDR_adjusted, Sig_adjusted = Significance_adjusted),
    by = c("Metric", "Variable")
  ) %>%
  mutate(
    R2_reduction_pct = ifelse(!is.na(R2_adjusted),
                              (R2_unadjusted - R2_adjusted) / R2_unadjusted * 100, NA),
    R2_ratio = NA_real_
  )

# Recalculate ratio properly
for (i in seq_len(nrow(merged))) {
  if (merged$Type[i] == "Confounder" & !is.na(merged$R2_adjusted[i])) {
    disease_r2 <- merged$R2_adjusted[merged$Variable == "Disease Status" & 
                                       merged$Metric == merged$Metric[i]]
    if (length(disease_r2) > 0 & merged$R2_adjusted[i] > 0) {
      merged$R2_ratio[i] <- disease_r2 / merged$R2_adjusted[i]
    }
  }
}

# ============================================================================
# 11. FIGURE 4 — TWO-PANEL MAIN FIGURE
# ============================================================================

cat("\n--- Generating Figure 4 ---\n")

# ---------------------------------------------------------------------------
# Panel A: Disease Status R² — Unadjusted vs Adjusted (grouped bar chart)
# ---------------------------------------------------------------------------

# Metric display order: Weighted first (highest R²)
metric_order <- c("Weighted UniFrac", "Unweighted UniFrac", "Bray-Curtis")
metric_labels <- c("Weighted UniFrac" = "Weighted\nUniFrac", 
                    "Unweighted UniFrac" = "Unweighted\nUniFrac",
                    "Bray-Curtis" = "Bray-Curtis")

disease_data <- merged %>%
  filter(Variable == "Disease Status") %>%
  select(Metric, R2_unadjusted, R2_adjusted, Sig_unadjusted, Sig_adjusted) %>%
  pivot_longer(cols = c(R2_unadjusted, R2_adjusted),
               names_to = "Model", values_to = "R2") %>%
  mutate(
    Model = recode(Model, 
                   "R2_unadjusted" = "Unadjusted",
                   "R2_adjusted"   = "Adjusted"),
    Model = factor(Model, levels = c("Unadjusted", "Adjusted")),
    Sig = ifelse(Model == "Unadjusted", Sig_unadjusted, Sig_adjusted),
    Sig_label = ifelse(Sig == "ns", "", Sig),
    R2_pct = R2 * 100,
    Metric = factor(Metric, levels = metric_order)
  )

panelA <- ggplot(disease_data, aes(x = Metric, y = R2_pct, fill = Model)) +
  geom_col(position = position_dodge(0.7), width = 0.5, color = "black", linewidth = 0.3) +
  geom_text(aes(label = Sig_label, y = R2_pct + 0.15),
            position = position_dodge(0.7), vjust = 0, size = 3.2,
            color = "black", fontface = "bold") +
  scale_fill_manual(values = colors_model, name = NULL) +
  scale_x_discrete(labels = metric_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)), limits = c(0, 8)) +
  labs(x = NULL, y = "Variance Explained (%)", tag = "A") +
  theme_pub +
  theme(legend.position = "top",
        legend.margin = margin(b = -3))

# ---------------------------------------------------------------------------
# Panel B: Individual Variable R² (faceted by metric)
# ---------------------------------------------------------------------------

# Group variables logically: Disease → Demographics/Clinical → Medications
# Within each group, sort by mean R² (descending)

demo_clinical <- c("Sex", "Age", "BMI", "TG", "LDL", "HDL")
medications   <- c("PPI Use", "Metformin", "Statin Use", "Rifaximin", "Lactulose")

mean_r2 <- merged %>%
  filter(Type == "Confounder") %>%
  group_by(Variable) %>%
  summarise(mean_r2 = mean(R2_adjusted, na.rm = TRUE))

demo_order <- mean_r2 %>% filter(Variable %in% demo_clinical) %>% arrange(desc(mean_r2)) %>% pull(Variable)
med_order  <- mean_r2 %>% filter(Variable %in% medications)   %>% arrange(desc(mean_r2)) %>% pull(Variable)

var_order <- c("Disease Status", demo_order, med_order)

plot_data_B <- merged %>%
  filter(!is.na(R2_adjusted)) %>%
  mutate(
    Variable = factor(Variable, levels = var_order),
    R2_pct = R2_adjusted * 100,
    Metric = factor(Metric, levels = metric_order),
    # Only show significant stars (hide "ns")
    Sig_label = ifelse(Sig_adjusted == "ns", "", Sig_adjusted)
  )

panelB <- ggplot(plot_data_B, aes(x = Variable, y = R2_pct, fill = Variable)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.6) +
  geom_text(data = plot_data_B %>% filter(Sig_label != ""),
            aes(label = Sig_label, y = R2_pct + 0.03),
            vjust = 0, size = 2.8, color = "black", fontface = "bold") +
  scale_fill_manual(values = colors_vars) +
  facet_wrap(~ Metric, nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = NULL, y = "Variance Explained (%)", tag = "B") +
  theme_pub +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 6.5, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  )

# ---------------------------------------------------------------------------
# Assemble Figure 4 (2-panel)
# ---------------------------------------------------------------------------

fig4 <- (panelA / panelB) +
  plot_layout(heights = c(0.45, 0.55)) +
  plot_annotation(
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

# ============================================================================
# EXPORT — Gut Microbes specifications
# Width: 174 mm (double-column / full page)
# Height: 130 mm (proportional for 2-row layout)
# Resolution: 600 dpi for raster; vector PDF/SVG for editable
# ============================================================================

ggsave(file.path(OUTPUT_DIR, "Figure4_Confounder_PERMANOVA.pdf"),
       fig4, width = 174, height = 130, units = "mm")

ggsave(file.path(OUTPUT_DIR, "Figure4_Confounder_PERMANOVA.png"),
       fig4, width = 174, height = 130, units = "mm",
       dpi = 600, bg = "white")

ggsave(file.path(OUTPUT_DIR, "Figure4_Confounder_PERMANOVA.svg"),
       fig4, width = 174, height = 130, units = "mm",
       device = svglite)

cat("  Figure 4 saved (Gut Microbes format):\n")
cat("    - Figure4_Confounder_PERMANOVA.pdf  (vector, submission-ready)\n")
cat("    - Figure4_Confounder_PERMANOVA.png  (600 dpi raster)\n")
cat("    - Figure4_Confounder_PERMANOVA.svg  (editable vector)\n")

# ============================================================================
# 12. EXCEL EXPORT
# ============================================================================

cat("\n--- Saving Excel results ---\n")

wb <- createWorkbook()

# Sheet 1: Univariate
addWorksheet(wb, "Univariate_PERMANOVA")
writeData(wb, "Univariate_PERMANOVA", univariate_results)

# Sheet 2: Adjusted
addWorksheet(wb, "Adjusted_PERMANOVA")
writeData(wb, "Adjusted_PERMANOVA", adjusted_results)

# Sheet 3: Merged comparison
addWorksheet(wb, "Unadjusted_vs_Adjusted")
writeData(wb, "Unadjusted_vs_Adjusted", merged)

# Sheet 4: Summary table (wide format for manuscript)
summary_wide <- merged %>%
  filter(!is.na(R2_adjusted)) %>%
  select(Variable, Metric, R2_unadjusted, R2_adjusted, 
         FDR_unadjusted, FDR_adjusted, R2_reduction_pct) %>%
  pivot_wider(
    names_from = Metric,
    values_from = c(R2_unadjusted, R2_adjusted, FDR_unadjusted, FDR_adjusted, R2_reduction_pct),
    names_sep = " | "
  )
addWorksheet(wb, "Summary_Wide")
writeData(wb, "Summary_Wide", summary_wide)

# Sheet 5: Parameters
params <- data.frame(
  Parameter = c("N_samples", "N_permutations", "Seed", "N_confounders",
                "Confounders", "Software", "FDR_method", "Date"),
  Value = c(length(common_samples), N_PERM, 42, length(confounders),
            paste(names(confounders), collapse = ", "),
            paste0("vegan ", packageVersion("vegan")),
            "Benjamini-Hochberg",
            as.character(Sys.Date()))
)
addWorksheet(wb, "Parameters")
writeData(wb, "Parameters", params)

saveWorkbook(wb, file.path(OUTPUT_DIR, "PERMANOVA_Confounder_Results.xlsx"), overwrite = TRUE)
cat("  Excel saved.\n")

# ============================================================================
# 13. SAVE RESULTS FOR SUPPLEMENTARY SCRIPT
# ============================================================================

save(univariate_results, adjusted_results, merged, meta_sub, dist_sub,
     common_samples, confounders, all_vars, group_colors, metric_colors,
     file = file.path(OUTPUT_DIR, "PERMANOVA_results_for_supplementary.RData"))

cat("  RData saved for supplementary script.\n")

# ============================================================================
# 14. MANUSCRIPT TEXT
# ============================================================================

cat("\n--- Generating manuscript text ---\n")

sink(file.path(OUTPUT_DIR, "Manuscript_PERMANOVA_Text.txt"))

cat("=== METHODS ===\n\n")
cat("Confounder Assessment\n\n")
cat(sprintf(
  "To evaluate the independent contribution of disease status to microbiome composition, 
we performed PERMANOVA (adonis2, vegan v%s) using three complementary beta diversity 
metrics: Bray-Curtis dissimilarity, unweighted UniFrac, and weighted UniFrac distances. 
Eleven potential confounders were assessed: sex, age, BMI, proton pump inhibitor (PPI) use, 
metformin use, rifaximin use, lactulose use, 
statin use, triglycerides (TG), low-density lipoprotein (LDL), and high-density lipoprotein 
(HDL). Recent antibiotic use was uniformly absent and therefore excluded. Isolated missing covariate values in three HCC participants (serum lipids; metformin status) were imputed with the cohort median (continuous) or mode (binary), retaining all 205 participants. Lipid parameters were included as confounders given the established role of 
dyslipidemia in MASLD pathogenesis. Univariate models tested each variable independently, 
followed by a full marginal model (Type III, by = \"margin\") including disease status and 
all confounders simultaneously to estimate independent (adjusted) variance contributions. 
All tests used %d permutations with a fixed seed for reproducibility. P-values were 
corrected for multiple testing using the Benjamini-Hochberg false discovery rate (FDR) 
within each distance metric.\n\n",
  packageVersion("vegan"), N_PERM
))

cat("=== RESULTS ===\n\n")
cat("Confounder Assessment\n\n")

# Extract key values
for (m in names(dist_sub)) {
  disease_unadi <- merged$R2_unadjusted[merged$Variable == "Disease Status" & merged$Metric == m]
  disease_adj   <- merged$R2_adjusted[merged$Variable == "Disease Status" & merged$Metric == m]
  disease_fdr   <- merged$FDR_adjusted[merged$Variable == "Disease Status" & merged$Metric == m]
  reduction     <- (disease_unadi - disease_adj) / disease_unadi * 100
  
  cat(sprintf("%s: Disease R² = %.4f (unadjusted) → %.4f (adjusted), %.1f%% reduction, FDR = %.4f\n",
              m, disease_unadi, disease_adj, reduction, disease_fdr))
}

cat("\nSignificant confounders (FDR < 0.05 in adjusted model):\n")
sig_conf <- merged %>% 
  filter(Type == "Confounder", FDR_adjusted < 0.05) %>%
  arrange(Metric, FDR_adjusted)
if (nrow(sig_conf) > 0) {
  for (i in seq_len(nrow(sig_conf))) {
    cat(sprintf("  %s — %s: R² = %.4f, FDR = %.4f\n",
                sig_conf$Metric[i], sig_conf$Variable[i], 
                sig_conf$R2_adjusted[i], sig_conf$FDR_adjusted[i]))
  }
} else {
  cat("  None\n")
}

cat("\nNon-significant confounders after FDR correction:\n")
ns_conf <- merged %>%
  filter(Type == "Confounder", FDR_adjusted >= 0.05) %>%
  pull(Variable) %>% unique()
cat(sprintf("  %s\n", paste(ns_conf, collapse = ", ")))

sink()
cat("  Manuscript text saved.\n")

# ============================================================================
# 15. CONSOLE SUMMARY
# ============================================================================

cat("\n")
cat("=" , rep("=", 70), "\n", sep = "")
cat("FIGURE 4 — PERMANOVA CONFOUNDER ANALYSIS COMPLETE\n")
cat("=" , rep("=", 70), "\n", sep = "")
cat(sprintf("Samples analyzed: %d\n", length(common_samples)))
cat(sprintf("Permutations: %d\n", N_PERM))
cat(sprintf("Confounders tested: %d\n", length(confounders)))
cat("\nKey Results per Metric:\n")

for (m in names(dist_sub)) {
  d_un <- merged$R2_unadjusted[merged$Variable == "Disease Status" & merged$Metric == m]
  d_ad <- merged$R2_adjusted[merged$Variable == "Disease Status" & merged$Metric == m]
  d_fdr <- merged$FDR_adjusted[merged$Variable == "Disease Status" & merged$Metric == m]
  red <- (d_un - d_ad) / d_un * 100
  cat(sprintf("  %s: R² %.4f → %.4f (−%.1f%%), FDR = %s\n",
              m, d_un, d_ad, red, 
              ifelse(d_fdr < 0.001, "<0.001", sprintf("%.4f", d_fdr))))
}

cat(sprintf("\nFiles generated in: %s\n", getwd()))
cat("  • Figure4_Confounder_PERMANOVA.pdf / .png / .svg\n")
cat("  • PERMANOVA_Confounder_Results.xlsx (5 sheets)\n")
cat("  • Manuscript_PERMANOVA_Text.txt\n")
cat("  • PERMANOVA_results_for_supplementary.RData\n")
cat("=" , rep("=", 70), "\n", sep = "")
