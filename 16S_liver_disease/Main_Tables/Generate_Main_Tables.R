################################################################################
# Generate_Main_Tables.R
# ──────────────────────────────────────────────────────────────────────────────
# Gut microbiome alterations across MASLD, cirrhosis, and HCC
# in a Saudi Arabian cohort -- MAIN TABLES
#
#   Table 1: Clinical and demographic characteristics by disease group
#   Table 2: Global PERMANOVA (unadjusted and adjusted for 11 covariates)
#
# Run from 16S_liver_disease/Main_Tables (any folder under Desktop works):
#   - metadata and distance files are found in Desktop/16S_reprocess
#   - Table 2 reuses Figure4_Confounder_PERMANOVA.R output
#     (PERMANOVA_results_for_supplementary.RData, found under 16S_liver_disease)
#     so the Bray-Curtis / UniFrac values are IDENTICAL to Figure 4;
#     Jaccard (not in Figure 4) is computed here with the same model.
# Outputs (this script's folder):
#   Table1_Clinical_Characteristics.docx
#   Table2_PERMANOVA.docx
#
# Required packages: readxl, dplyr, flextable, officer, vegan
################################################################################

# ══════════════════════════════════════════════════════════════════════════════
# 0. LIBRARIES & PATHS
# ══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(flextable); library(officer); library(vegan)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f) == 1 && nzchar(f)) return(dirname(normalizePath(f)))
  for (i in rev(seq_len(sys.nframe()))) { of <- sys.frame(i)$ofile; if (!is.null(of)) return(dirname(normalizePath(of))) }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getSourceEditorContext()$path; if (nzchar(p)) return(dirname(normalizePath(p))) }
  getwd()
}
setwd(get_script_dir())
find_proj_dir <- function(name = "16S_reprocess") {
  d <- normalizePath(getwd())
  repeat { cand <- file.path(d, name); if (dir.exists(cand)) return(cand); up <- dirname(d); if (up == d) break; d <- up }
  stop("Could not find a '", name, "' folder above ", getwd())
}
PROJ_DIR      <- find_proj_dir()                                   # Desktop/16S_reprocess
ANALYSIS_DIR <- dirname(getwd())         # 16S_liver_disease
METADATA_FILE <- file.path(PROJ_DIR, "sample_metadata.tsv")
JACCARD_FILE  <- file.path(PROJ_DIR, "qiime2_analysis_v2", "Beta_diversity", "Jaccard", "Jaccardraw_data.tsv")
rd <- list.files(ANALYSIS_DIR, pattern = "^PERMANOVA_results_for_supplementary\\.RData$", recursive = TRUE, full.names = TRUE)
RDATA_FILE    <- if (length(rd)) rd[which.max(file.mtime(rd))] else NA
OUTPUT_DIR    <- file.path(ANALYSIS_DIR, "Main_Tables")             # 16S_liver_disease/Main_Tables
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
OUT_T1 <- file.path(OUTPUT_DIR, "Table1_Clinical_Characteristics.docx")
OUT_T2 <- file.path(OUTPUT_DIR, "Table2_PERMANOVA.docx")
GROUP_ORDER <- c("Healthy", "MASLD", "Cirrhosis", "HCC")
N_PERM <- 999

for (f in c(METADATA_FILE, JACCARD_FILE)) if (!file.exists(f)) stop("Input not found: ", f)
if (is.na(RDATA_FILE)) stop("PERMANOVA_results_for_supplementary.RData not found under ", ANALYSIS_DIR,
                            "\nRun Figure4_Confounder_PERMANOVA.R first.")
cat("Metadata :", METADATA_FILE, "\nFigure 4 :", RDATA_FILE, "\nJaccard  :", JACCARD_FILE, "\n\n")

# Word post-processing (compact layout, no heading numbering) -- shared by both tables
compact_docx <- function(path) {
  tmp <- tempfile("docx_fix_"); dir.create(tmp); utils::unzip(path, exdir = tmp)
  sp <- file.path(tmp, "word", "styles.xml"); sx <- readChar(sp, file.info(sp)$size, useBytes = TRUE)
  sx <- gsub("<w:numPr>\\s*<w:ilvl w:val=\"\\d+\"/>\\s*<w:numId w:val=\"\\d+\"/>\\s*</w:numPr>", "", sx)
  writeChar(sx, sp, eos = NULL, useBytes = TRUE)
  dp <- file.path(tmp, "word", "document.xml"); dx <- readChar(dp, file.info(dp)$size, useBytes = TRUE)
  for (s in c('w:after="100" w:before="100"', 'w:after="40" w:before="40"', 'w:after="20" w:before="20"'))
    dx <- gsub(s, 'w:after="0" w:before="0"', dx, fixed = TRUE)
  dx <- gsub('<w:trHeight w:val="360" w:hRule="auto"/>', '<w:trHeight w:val="160" w:hRule="auto"/>', dx, fixed = TRUE)
  dx <- gsub('<w:pgMar w:header="708" w:bottom="1417" w:top="1417" w:right="1417" w:left="1417" w:footer="708" w:gutter="0"/>',
             '<w:pgMar w:header="360" w:bottom="720" w:top="720" w:right="1000" w:left="1000" w:footer="360" w:gutter="0"/>', dx, fixed = TRUE)
  writeChar(dx, dp, eos = NULL, useBytes = TRUE)
  old <- setwd(tmp); out <- normalizePath(path, mustWork = FALSE); if (file.exists(out)) file.remove(out)
  system2("zip", args = c("-r", "-X", "-q", shQuote(out), "."), stdout = TRUE, stderr = TRUE)
  if (!file.exists(out)) utils::zip(out, files = list.files(".", recursive = TRUE), flags = "-r9Xq")
  setwd(old); unlink(tmp, recursive = TRUE)
}
new_doc <- function() read_docx() %>% body_set_default_section(prop_section(
  page_size = page_size(width = 11, height = 8.5, orient = "landscape"),
  page_margins = page_mar(top = 0.5, bottom = 0.5, left = 0.7, right = 0.7)))

cat("==========================================================\n  TABLE 1\n==========================================================\n")
# ══════════════════════════════════════════════════════════════════════════════
# 2. HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

fmt_med_iqr <- function(x) {
  x <- as.numeric(x[!is.na(x)])
  if (length(x) == 0) return("--")
  sprintf("%.1f (%.1f-%.1f)", median(x), quantile(x, 0.25), quantile(x, 0.75))
}

fmt_n_pct <- function(x, level = "Yes") {
  # Denominator = full group n (missing counted as non-users), matching
  # Supp Table S1 and the imputation used downstream. Using sum(!is.na(x))
  # here mis-counted HCC metformin as 2/16 = 12.5% instead of 2/17 = 11.8%.
  n_total <- length(x)
  n_yes   <- sum(x == level, na.rm = TRUE)
  if (n_total == 0) return("0 (0.0)")
  sprintf("%d (%.1f)", n_yes, 100 * n_yes / n_total)
}

fmt_p <- function(p) {
  if (is.na(p)) return("--")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. LOAD METADATA (Table 1)
# ══════════════════════════════════════════════════════════════════════════════

cat("-- Loading metadata --\n")
meta <- read.delim(METADATA_FILE, check.names = FALSE, stringsAsFactors = FALSE)
meta <- as.data.frame(meta, stringsAsFactors = FALSE)
colnames(meta)[1] <- "SampleID"
meta$Status <- gsub("^NAFLD$", "MASLD", meta$Status)
meta$Status <- factor(meta$Status, levels = GROUP_ORDER)

num_cols <- c("Age", "BMI", "ALT", "AST", "ALP", "Albumin", "Globulin",
              "Bili_T", "GGT", "AFP", "Glu_R", "Hgb_A1C", "TG", "TC",
              "HDL", "LDL", "BUN", "Creatinine", "Postassium", "Sodium",
              "Hgb", "WBC", "Platelet")
for (col in num_cols) {
  if (col %in% colnames(meta)) meta[[col]] <- as.numeric(meta[[col]])
}
cat(sprintf("  %d samples loaded\n", nrow(meta)))

# ══════════════════════════════════════════════════════════════════════════════
# 4. BUILD TABLE 1
# ══════════════════════════════════════════════════════════════════════════════

cat("\n-- Building Table 1 --\n")

var_labs <- c(
  Age = "Age (years)", BMI = "BMI (kg/m2)",
  ALT = "ALT (U/L)", AST = "AST (U/L)", ALP = "ALP (U/L)",
  Albumin = "Albumin (g/L)", Globulin = "Globulin (g/L)",
  Bili_T = "Total bilirubin (umol/L)", GGT = "GGT (U/L)",
  AFP = "AFP (ng/mL)",
  Hgb_A1C = "HbA1c (%)", TG = "Triglycerides (mmol/L)",
  TC = "Total cholesterol (mmol/L)", HDL = "HDL (mmol/L)",
  LDL = "LDL (mmol/L)", BUN = "BUN (mmol/L)",
  Creatinine = "Creatinine (umol/L)", Postassium = "Potassium (mmol/L)",
  Sodium = "Sodium (mmol/L)", Hgb = "Hemoglobin (g/dL)",
  WBC = "WBC (x10^9/L)", Platelet = "Platelets (x10^9/L)"
)

# Microbiome-relevant medications only
med_vars <- c(
  PPI_Use = "PPI",
  Metformin = "Metformin",
  Antibiotic_Use_Last3mo = "Antibiotics (last 3 mo)",
  Rifaximin_Use = "Rifaximin",
  Lactulose_Use = "Lactulose",
  Probiotic_Use_Last2wk = "Probiotics (last 2 wk)",
  Statin_Use = "Statin"
)

rows <- list()
raw_pvals <- c()     # collect all raw p-values for FDR correction
pval_indices <- c()  # track which rows have testable p-values
ri <- 1

# n row
rows[[ri]] <- data.frame(Variable = "n", Healthy = "70", MASLD = "59",
                         Cirrhosis = "59", HCC = "17", P_raw = NA, Test = "",
                         stringsAsFactors = FALSE)
ri <- ri + 1

# Gender
set.seed(2024)
p_sex <- tryCatch(chisq.test(table(meta$Status, meta$Gender),
                             simulate.p.value = TRUE, B = 10000)$p.value,
                  error = function(e) NA)
rows[[ri]] <- data.frame(
  Variable = "Gender -- Male, n (%)",
  Healthy = fmt_n_pct(meta$Gender[meta$Status == "Healthy"], "M"),
  MASLD = fmt_n_pct(meta$Gender[meta$Status == "MASLD"], "M"),
  Cirrhosis = fmt_n_pct(meta$Gender[meta$Status == "Cirrhosis"], "M"),
  HCC = fmt_n_pct(meta$Gender[meta$Status == "HCC"], "M"),
  P_raw = p_sex,
  Test = "Chi2", stringsAsFactors = FALSE)
ri <- ri + 1

# Disease etiology (Reviewer 2, comment 2) -- descriptive only, not tested
# (etiology partly defines the groups). Source columns:
#   Causes_of_Diseases   : "Viral infection_*" = viral; otherwise non-viral
#   Causes_of_Diseases_2 : HCV / HBV subtype (incl. "HCV+NAFLD", "HBV+NASH")
viral <- grepl("^Viral", meta$Causes_of_Diseases)
hcv   <- viral & grepl("HCV", meta$Causes_of_Diseases_2)
hbv   <- viral & grepl("HBV", meta$Causes_of_Diseases_2)
etio_row <- function(label, flag) {
  cells <- sapply(GROUP_ORDER, function(g) {
    if (g == "Healthy") return("--")
    s <- meta$Status == g
    sprintf("%d (%.1f)", sum(flag[s]), 100 * mean(flag[s]))
  })
  data.frame(Variable = label, Healthy = cells[["Healthy"]], MASLD = cells[["MASLD"]],
             Cirrhosis = cells[["Cirrhosis"]], HCC = cells[["HCC"]],
             P_raw = NA, Test = "--", stringsAsFactors = FALSE)
}
rows[[ri]] <- data.frame(Variable = "Disease etiology, n (%)", Healthy = "", MASLD = "",
                         Cirrhosis = "", HCC = "", P_raw = NA, Test = "", stringsAsFactors = FALSE)
ri <- ri + 1
rows[[ri]] <- etio_row("  Non-viral", !viral); ri <- ri + 1
rows[[ri]] <- etio_row("  Viral hepatitis", viral);        ri <- ri + 1
rows[[ri]] <- etio_row("    HCV", hcv);                    ri <- ri + 1
rows[[ri]] <- etio_row("    HBV", hbv);                    ri <- ri + 1
stopifnot(sum(viral) == 25, sum(hcv) == 17, sum(hbv) == 8)

# Continuous variables
for (var in names(var_labs)) {
  if (!var %in% colnames(meta)) next
  vals <- suppressWarnings(as.numeric(as.character(meta[[var]])))
  p <- tryCatch(kruskal.test(vals ~ meta$Status)$p.value,
                error = function(e) NA)
  rows[[ri]] <- data.frame(
    Variable = var_labs[var],
    Healthy = fmt_med_iqr(vals[meta$Status == "Healthy"]),
    MASLD = fmt_med_iqr(vals[meta$Status == "MASLD"]),
    Cirrhosis = fmt_med_iqr(vals[meta$Status == "Cirrhosis"]),
    HCC = fmt_med_iqr(vals[meta$Status == "HCC"]),
    P_raw = p,
    Test = "KW", stringsAsFactors = FALSE)
  ri <- ri + 1
}

# Medications header
rows[[ri]] <- data.frame(Variable = "Medications, n (%)", Healthy = "",
                         MASLD = "", Cirrhosis = "", HCC = "",
                         P_raw = NA, Test = "", stringsAsFactors = FALSE)
ri <- ri + 1

# Medication variables
for (col in names(med_vars)) {
  if (!col %in% colnames(meta)) next
  if (length(unique(meta[[col]][!is.na(meta[[col]])])) <= 1) {
    p <- NA
  } else {
    set.seed(2024)
    p <- tryCatch(chisq.test(table(meta$Status, meta[[col]]),
                             simulate.p.value = TRUE, B = 10000)$p.value,
                  error = function(e) NA)
  }
  rows[[ri]] <- data.frame(
    Variable = paste0("  ", med_vars[col]),
    Healthy = fmt_n_pct(meta[[col]][meta$Status == "Healthy"]),
    MASLD = fmt_n_pct(meta[[col]][meta$Status == "MASLD"]),
    Cirrhosis = fmt_n_pct(meta[[col]][meta$Status == "Cirrhosis"]),
    HCC = fmt_n_pct(meta[[col]][meta$Status == "HCC"]),
    P_raw = p,
    Test = ifelse(is.na(p), "--", "Chi2"), stringsAsFactors = FALSE)
  ri <- ri + 1
}

st1 <- bind_rows(rows)

# ── Report RAW omnibus p-values (no multiple-testing correction) ──
# Baseline characteristics are descriptive; raw p keeps Table 1 and Supplementary
# Table S1 identical and avoids a panel-dependent BH family.
st1$P <- sapply(st1$P_raw, function(p) {
  if (is.na(p)) return("--")
  paste0(fmt_p(p), stars(p))
})

# Drop raw columns
st1 <- st1 %>% select(Variable, Healthy, MASLD, Cirrhosis, HCC, P, Test)

cat(sprintf("  %d rows (%d tests, raw p-values)\n", nrow(st1), sum(!is.na(st1$P) & st1$P != "--")))

# ══════════════════════════════════════════════════════════════════════════════
# 5. FORMAT TABLE 1
# ══════════════════════════════════════════════════════════════════════════════

ft_st1 <- flextable(st1) %>%
  set_header_labels(Variable = "Variable",
                    Healthy = "Healthy\n(n = 70)",
                    MASLD = "MASLD\n(n = 59)",
                    Cirrhosis = "Cirrhosis\n(n = 59)",
                    HCC = "HCC\n(n = 17)",
                    P = "P-value", Test = "Test") %>%
  fontsize(size = 8, part = "all") %>%
  font(fontname = "Arial", part = "all") %>%
  bg(part = "header", bg = "#2C3E50") %>%
  color(part = "header", color = "white") %>%
  bold(part = "header") %>%
  align(j = 2:7, align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "body") %>%
  width(j = 1, width = 2.4) %>%
  width(j = 2:5, width = 1.6) %>%
  width(j = 6, width = 0.9) %>%
  width(j = 7, width = 0.5) %>%
  border_outer(border = fp_border(color = "#999999")) %>%
  border_inner_h(border = fp_border(color = "#DDDDDD", width = 0.5)) %>%
  padding(padding = 2, part = "body")

# Bold section headers and significant p-values
ft_st1 <- bold(ft_st1, i = which(st1$Variable %in% c("n", "Medications, n (%)", "Disease etiology, n (%)")), j = 1)
for (i in seq_len(nrow(st1))) {
  if (grepl("\\*", st1$P[i])) ft_st1 <- bold(ft_st1, i = i, j = 6)
}


# ══════════════════════════════════════════════════════════════════════════════
# 6. EXPORT TABLE 1
# ══════════════════════════════════════════════════════════════════════════════
doc1 <- new_doc() %>%
  body_add_par("Table 1. Clinical and Demographic Characteristics by Disease Group", style = "heading 2") %>%
  body_add_par(paste0(
    "Data presented as median (IQR) for continuous variables and n (%) for categorical variables. ",
    "P values: Kruskal-Wallis (KW) test for continuous variables and chi-squared test with Monte Carlo simulation ",
    "(10,000 replicates, seed 2024) for categorical variables; P values are unadjusted. ",
    "Disease etiology is reported descriptively and was not tested because it partly defines the clinical groups; ",
    "two participants with viral cirrhosis had coexisting MASLD and were classified as viral. ",
    "Antibiotic (3 months) and probiotic (2 weeks) use were exclusion criteria (0% in all groups), so no test was performed. ",
    "AFP was not measured in healthy controls. --, not applicable or not tested. ",
    "Detailed liver disease severity (Child-Pugh class, ALBI grade, MELD), MASLD cardiometabolic criteria, ",
    "FibroScan steatosis and fibrosis grades, and additional medications are given in Supplementary Table S1. ",
    "*p < 0.05, **p < 0.01, ***p < 0.001."), style = "Normal") %>%
  body_add_flextable(ft_st1)
invisible(body_end_section_landscape(doc1)); print(doc1, target = OUT_T1); compact_docx(OUT_T1)
cat("  Saved:", OUT_T1, "\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# 7. TABLE 2 -- GLOBAL PERMANOVA
# ══════════════════════════════════════════════════════════════════════════════
cat("==========================================================\n  TABLE 2\n==========================================================\n")
load(RDATA_FILE)   # merged, meta_sub, dist_sub, confounders (from Figure 4)

# Jaccard distance aligned to the Figure 4 samples
jl  <- read.delim(JACCARD_FILE)
ids <- rownames(meta_sub)
Mj  <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
k   <- jl$SubjectID1 %in% ids & jl$SubjectID2 %in% ids
Mj[cbind(match(jl$SubjectID1[k], ids), match(jl$SubjectID2[k], ids))] <- jl$Distance[k]
Mj[cbind(match(jl$SubjectID2[k], ids), match(jl$SubjectID1[k], ids))] <- jl$Distance[k]
dist_all <- c(dist_sub, list("Jaccard" = as.dist(Mj)))

rhs_full <- paste(c("Status", unname(confounders)), collapse = " + ")
omega2 <- function(tab, term) {             # df-adjusted effect size from an adonis2 table
  ms_res <- tab["Residual", "SumOfSqs"] / tab["Residual", "Df"]
  max(0, (tab[term, "SumOfSqs"] - tab[term, "Df"] * ms_res) / (tab["Total", "SumOfSqs"] + ms_res))
}

t2_rows <- lapply(c("Weighted UniFrac", "Unweighted UniFrac", "Bray-Curtis", "Jaccard"), function(m) {
  d <- dist_all[[m]]
  # SS, F and R2 do not depend on permutations; 9 permutations only to obtain the table quickly
  full0 <- adonis2(as.formula(paste("d ~", rhs_full)), data = meta_sub, permutations = 9, by = "margin")
  if (m %in% merged$Metric) {               # identical to Figure 4
    r <- merged[merged$Metric == m & merged$Variable == "Disease Status", ]
    out <- data.frame(Metric = m, F_u = r$F_unadjusted, R2_u = r$R2_unadjusted, P_u = r$P_unadjusted,
                      F_a = r$F_adjusted, R2_a = r$R2_adjusted, P_a = r$P_adjusted, Q_a = r$FDR_adjusted)
  } else {                                  # Jaccard: same models, computed here
    set.seed(42)
    uni  <- adonis2(d ~ Status, data = meta_sub, permutations = N_PERM, by = "terms")   # by = "terms": vegan >= 2.7 otherwise labels the row "Model"
    full <- adonis2(as.formula(paste("d ~", rhs_full)), data = meta_sub, permutations = N_PERM, by = "margin")
    terms <- setdiff(rownames(full), c("Residual", "Total"))
    q <- p.adjust(full[terms, "Pr(>F)"], "BH"); names(q) <- terms
    out <- data.frame(Metric = m, F_u = uni["Status", "F"], R2_u = uni["Status", "R2"], P_u = uni["Status", "Pr(>F)"],
                      F_a = full["Status", "F"], R2_a = full["Status", "R2"], P_a = full["Status", "Pr(>F)"], Q_a = q[["Status"]])
  }
  out$Omega2 <- omega2(full0, "Status")
  out
})
t2 <- bind_rows(t2_rows) %>% mutate(Reduction = 100 * (R2_u - R2_a) / R2_u)
if (anyNA(t2[, c("F_u", "R2_u", "P_u", "F_a", "R2_a", "P_a")])) stop("Table 2 has missing values -- check the adonis2 output above")
print(t2, digits = 3)

fmt_pp <- function(p) if (is.na(p)) "NA" else if (p < 0.001) "<0.001" else sprintf("%.3f", p)
st2 <- t2 %>% transmute(
  Metric = sub("Bray-Curtis", "Bray\u2013Curtis", Metric),
  F_u = sprintf("%.2f", F_u), R2_u = sprintf("%.3f", R2_u), P_u = paste0(sapply(P_u, fmt_pp), sapply(P_u, stars)),
  F_a = sprintf("%.2f", F_a), R2_a = sprintf("%.3f", R2_a), P_a = paste0(sapply(P_a, fmt_pp), sapply(P_a, stars)),
  Q_a = sprintf("%.3f", Q_a), Om = sprintf("%.4f", Omega2), Red = sprintf("%.1f", Reduction))

ft_st2 <- flextable(st2) %>%
  set_header_labels(Metric = "Distance metric", F_u = "Pseudo-F", R2_u = "R\u00b2", P_u = "P",
                    F_a = "Pseudo-F", R2_a = "R\u00b2", P_a = "P", Q_a = "FDR q", Om = "\u03c9\u00b2",
                    Red = "R\u00b2 reduction (%)") %>%
  add_header_row(values = c("", "Unadjusted (disease status only)", "Adjusted for 11 covariates (marginal)", ""),
                 colwidths = c(1, 3, 5, 1), top = TRUE) %>%
  fontsize(size = 9, part = "all") %>% font(fontname = "Arial", part = "all") %>%
  bg(part = "header", bg = "#2C3E50") %>% color(part = "header", color = "white") %>% bold(part = "header") %>%
  align(j = 2:10, align = "center", part = "all") %>% align(j = 1, align = "left", part = "body") %>%
  width(j = 1, width = 1.7) %>% width(j = 2:9, width = 0.85) %>% width(j = 10, width = 1.1) %>%
  border_outer(border = fp_border(color = "#999999")) %>%
  border_inner_h(border = fp_border(color = "#DDDDDD", width = 0.5)) %>% padding(padding = 2, part = "body")

doc2 <- new_doc() %>%
  body_add_par("Table 2. Global PERMANOVA of beta diversity across disease groups (healthy, MASLD, cirrhosis, HCC).",
               style = "heading 2") %>%
  body_add_flextable(ft_st2) %>%
  body_add_par(paste0(
    "PERMANOVA was performed with adonis2 (vegan v", as.character(packageVersion("vegan")), ") on distances from ",
    "samples rarefied to 13,000 reads (n = 205). P values were obtained by permutation (", N_PERM, " permutations). ",
    "Unadjusted: distance ~ disease status. Adjusted: disease status and 11 covariates (sex, age, BMI, TG, LDL, HDL, ",
    "and PPI, metformin, statin, rifaximin, and lactulose use) fitted simultaneously with marginal sums of squares ",
    "(by = \"margin\"), so each R\u00b2 is the variance uniquely attributable to that term. FDR q, Benjamini-Hochberg ",
    "correction across the 12 terms of the adjusted model within each metric. \u03c9\u00b2, degrees-of-freedom-adjusted ",
    "effect size of disease status in the adjusted model; disease status had the largest \u03c9\u00b2 of any term for every ",
    "metric (Supplementary Figure S3C). R\u00b2 reduction, percentage decrease in disease R\u00b2 after adjustment. ",
    "Weighted UniFrac, unweighted UniFrac, and Bray-Curtis values are identical to Figure 4; Jaccard was analysed with ",
    "the same models. With n = 205, small R\u00b2 values can reach statistical significance; statistical significance ",
    "does not imply a large biological effect. *p < 0.05, **p < 0.01, ***p < 0.001."), style = "Normal")
invisible(body_end_section_landscape(doc2)); print(doc2, target = OUT_T2); compact_docx(OUT_T2)
cat("  Saved:", OUT_T2, "\n\nDONE\n")
