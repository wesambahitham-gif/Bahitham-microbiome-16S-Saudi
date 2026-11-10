################################################################################
# Generate_Supplementary_Tables.R
# ──────────────────────────────────────────────────────────────────────────────
# Gut microbiome alterations across MASLD, cirrhosis, and HCC
# in a Saudi Arabian cohort (Microbiology Spectrum revision)
#
# Generates publication-ready Word document with supplementary tables:
#   ST1:  Extended Clinical & Demographic Characteristics by Disease Group
#   ST2:  Alpha Diversity Indices and Pairwise Comparisons
#   ST3:  ANCOM-BC2 Full Differential Abundance Results (Figure 6 script)
#   ST4:  ANCOM-BC2 / DESeq2 Concordance (Supp Fig S2 script)
#   ST5:  Spearman Correlations, pooled and group-adjusted (Figure 7 script)
#   ST6:  PERMANOVA Variance Partitioning
#   (ST7 = STORMS checklist, separate file)
#   ST8:  Sensitivity analysis excluding viral etiology         (R2 comment 2)
#   ST9:  Age, BMI and sex: PERMANOVA decomposition               (R2 comments 3-5)
#   ST10: Covariate-adjusted differential abundance              (R2 comment 6)
#   ST11: Dispersion-robust PERMANOVA (Welch W*d, balanced)     (R2 comment 9)
#
# RUN ORDER: Figure6 -> Figure7 -> Figure6 -> Supp Fig S2 ->
#            Supplementary_Analyses_Revision.R -> THIS SCRIPT
#
# Required packages:
#   install.packages(c("readxl", "dplyr", "tidyr", "flextable",
#                      "officer", "dunn.test", "vegan"))
################################################################################

# ══════════════════════════════════════════════════════════════════════════════
# 0. LIBRARIES
# ══════════════════════════════════════════════════════════════════════════════

library(readxl)
library(dplyr)
library(tidyr)
library(flextable)
library(officer)
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

if (!requireNamespace("dunn.test", quietly = TRUE)) install.packages("dunn.test")
library(dunn.test)
if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan")
library(vegan)
# Guard against masking (MASS, ANCOMBC deps, plyr, AnnotationDbi also export select/filter)
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename; mutate <- dplyr::mutate

cat("==========================================================\n")
cat("  Supplementary Tables S1-S6 -- Generation Script\n")
cat("==========================================================\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# 1. CONFIGURATION -- UPDATE FILE PATHS HERE
# ══════════════════════════════════════════════════════════════════════════════

# Paths are RELATIVE: inputs are found from this script's folder; output saved here
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
PROJ_DIR      <- find_proj_dir()                                   # Desktop/16S_reprocess
DATA_DIR      <- file.path(PROJ_DIR, "qiime2_analysis_v2")         # reprocessed data
REVIEWERS_DIR <- file.path(dirname(PROJ_DIR), "Reviewers")         # Desktop/Reviewers
ALPHA_DIR     <- file.path(DATA_DIR, "Alpha_diversity")

# Input files
METADATA_FILE  <- file.path(PROJ_DIR, "sample_metadata.tsv")
EXCL_FILE      <- file.path(PROJ_DIR, "sample_exclusions.txt")     # optional; shared with all scripts
ANCOMBC2_DIR   <- file.path(DATA_DIR, "ANCOMBC2_R")                # Figure 6 / 7 outputs
REV_DIR        <- file.path(REVIEWERS_DIR, "Supp_Analyses")        # Supplementary_Analyses_Revision.R outputs

# Alpha diversity: QIIME 2 exports (rarefied at 13,000), read directly
SHANNON_FILE   <- file.path(ALPHA_DIR, "Shannon", "Shannon_Liver.tsv")
CHAO1_FILE     <- file.path(ALPHA_DIR, "Chao1", "Chao1.tsv")
OBSERVED_FILE  <- file.path(ALPHA_DIR, "Observed feature", "Observed_feature_metadata.tsv")
SIMPSON_FILE   <- file.path(ALPHA_DIR, "Simpson", "Simpson_liver_micorbiome.tsv")
FAITHPD_FILE   <- file.path(ALPHA_DIR, "Faith-pd", "faith-pd.tsv")

# Differential abundance
ANCOMBC_FILE   <- file.path(ANCOMBC2_DIR, "ANCOMBC2_all_results.csv")            # Figure6_ANCOMBC2_Heatmap.R
s4_hits <- list.files(REVIEWERS_DIR, pattern = "^Table_S4_Concordance_Details\\.xlsx$",
                      recursive = TRUE, full.names = TRUE)                       # Supp Fig S2 script output
DESEQ2_FILE    <- if (length(s4_hits)) s4_hits[which.max(file.mtime(s4_hits))] else
                  file.path(REVIEWERS_DIR, "Table_S4_Concordance_Details.xlsx")

# Correlations & PERMANOVA
SPEARMAN_FILE  <- file.path(ANCOMBC2_DIR, "Spearman_pooled_and_partial.csv")     # Figure7_Spearman_Heatmap.R
# Revision tables (Supplementary_Analyses_Revision.R)
REV03_FILE     <- file.path(REV_DIR, "Rev03_PERMANOVA_Robustness.xlsx")           # ST8
REV02_FILE     <- file.path(REV_DIR, "Rev02_DA_Adjusted_Sensitivity.xlsx")        # ST8, ST10
S9_FILE        <- file.path(REV_DIR, "Supp_S9_Age_BMI_Sex_PERMANOVA.xlsx")        # ST9
# Table S6 (PERMANOVA) is computed IN-SCRIPT from the QIIME2 distance exports + metadata
BETA_DIR <- file.path(DATA_DIR, "Beta_diversity")
BRAY_TSV <- file.path(BETA_DIR, "Bray-curtis",        "Bray-Curtis_raw_data.tsv")
UUF_TSV  <- file.path(BETA_DIR, "UnWeighted_UniFrac",  "Unweighted.tsv")
WUF_TSV  <- file.path(BETA_DIR, "Weighted_UniFrac",    "Weighted_rowdata.tsv")

# Output (this script's folder)
OUTPUT_DIR  <- getwd()
OUTPUT_FILE <- file.path(OUTPUT_DIR, "Supplementary_Tables_S1-S6.docx")
XLSX_S4     <- file.path(OUTPUT_DIR, "Supplementary_Table_S4_Differential_Abundance.xlsx")

# Constants
GROUP_ORDER <- c("Healthy", "MASLD", "Cirrhosis", "HCC")
FDR_THRESH  <- 0.05

# ── Verify all input files exist ──
all_files <- c(METADATA_FILE, SHANNON_FILE, CHAO1_FILE, OBSERVED_FILE,
               SIMPSON_FILE, FAITHPD_FILE, ANCOMBC_FILE, DESEQ2_FILE,
               SPEARMAN_FILE, BRAY_TSV, UUF_TSV, WUF_TSV, REV03_FILE, REV02_FILE, S9_FILE)
missing <- all_files[!file.exists(all_files)]
if (length(missing) > 0) {
  cat("\n!!! MISSING FILES !!!\n")
  for (f in missing) cat("  NOT FOUND:", f, "\n")
  cat("\nRun the figure / analysis scripts that create these files first.\n")
  stop("Missing input files. See list above.")
}
cat("  All", length(all_files), "input files found.\n")

# ══════════════════════════════════════════════════════════════════════════════
# 2. HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

fmt_med_iqr <- function(x) {
  x <- as.numeric(x[!is.na(x)])
  if (length(x) == 0) return("--")
  sprintf("%.1f (%.1f-%.1f)", median(x), quantile(x, 0.25), quantile(x, 0.75))
}

fmt_n_pct <- function(x, level = "Yes") {
  n_total <- length(x)          # full group n (missing = non-user), identical to Table 1
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

epsilon_sq <- function(H, n, k) { (H - k + 1) / (n - k) }

# ══════════════════════════════════════════════════════════════════════════════
# 3. LOAD METADATA
# ══════════════════════════════════════════════════════════════════════════════

cat("-- Loading metadata --\n")
meta <- read.delim(METADATA_FILE, sep = "\t", stringsAsFactors = FALSE,
                   check.names = FALSE)
colnames(meta)[1] <- "SampleID"
meta$Status <- gsub("^NAFLD$", "MASLD", meta$Status)
meta$Status <- factor(meta$Status, levels = GROUP_ORDER)
meta <- meta[!grepl("^#", meta$SampleID), ]
if (file.exists(EXCL_FILE)) {
  excl <- trimws(readLines(EXCL_FILE)); meta <- meta[!meta$SampleID %in% excl[nzchar(excl)], ]
  cat("  Excluded samples from sample_exclusions.txt\n")
}
GROUP_N <- setNames(as.integer(table(meta$Status)), GROUP_ORDER)

num_cols <- c("Age", "BMI", "ALT", "AST", "ALP", "Albumin", "Globulin",
              "Bili_T", "GGT", "AFP", "Glu_R", "Hgb_A1C", "TG", "TC",
              "HDL", "LDL", "BUN", "Creatinine", "Postassium", "Sodium",
              "Hgb", "WBC", "Platelet")
for (col in num_cols) {
  if (col %in% colnames(meta)) meta[[col]] <- as.numeric(meta[[col]])
}
cat(sprintf("  %d samples loaded\n", nrow(meta)))


###############################################################################
#   ST1: EXTENDED CLINICAL & DEMOGRAPHIC CHARACTERISTICS                      #
###############################################################################

cat("\n-- ST1: Clinical Characteristics --\n")

var_labs <- c(
  Age = "Age (years)", BMI = "BMI (kg/m2)",
  ALT = "ALT (U/L)", AST = "AST (U/L)", ALP = "ALP (U/L)",
  Albumin = "Albumin (g/L)", Globulin = "Globulin (g/L)",
  Bili_T = "Total bilirubin (umol/L)", GGT = "GGT (U/L)",
  AFP = "AFP (ng/mL)", Glu_R = "Random glucose (mmol/L)",
  Hgb_A1C = "HbA1c (%)", TG = "Triglycerides (mmol/L)",
  TC = "Total cholesterol (mmol/L)", HDL = "HDL (mmol/L)",
  LDL = "LDL (mmol/L)", BUN = "BUN (mmol/L)",
  Creatinine = "Creatinine (umol/L)", Postassium = "Potassium (mmol/L)",
  Sodium = "Sodium (mmol/L)", Hgb = "Hemoglobin (g/dL)",
  WBC = "WBC (x10^9/L)", Platelet = "Platelets (x10^9/L)"
)

med_vars <- c(
  Metformin = "Metformin", Insulin = "Insulin",
  Statin_Use = "Statin", PPI_Use = "PPI",
  Rifaximin_Use = "Rifaximin", Lactulose_Use = "Lactulose",
  BetaBlocker_Use = "Beta-blocker", Diuretic = "Diuretic",
  Antibiotic_Use_Last3mo = "Antibiotics (3 mo)",
  antihypertensive = "Antihypertensive",
  VIT_D_Supplement = "Vitamin D"
)

rows <- list()
ri <- 1

# n row
rows[[ri]] <- data.frame(Variable = "n", Healthy = as.character(GROUP_N[["Healthy"]]),
                          MASLD = as.character(GROUP_N[["MASLD"]]), Cirrhosis = as.character(GROUP_N[["Cirrhosis"]]),
                          HCC = as.character(GROUP_N[["HCC"]]), P = "", Test = "",
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
  P = paste0(fmt_p(p_sex), stars(p_sex)),
  Test = "Chi2", stringsAsFactors = FALSE)
ri <- ri + 1

# ── Disease etiology (M7): non-viral (metabolic / cirrhosis-derived) vs viral (HBV/HCV) ──
# Healthy has no liver disease (shown as "--"); MASLD is metabolic by definition.
etio_viral <- grepl("Viral", meta$Causes_of_Diseases, ignore.case = TRUE)
etio_np <- function(grp, viral) {
  idx <- meta$Status == grp
  n <- sum(idx); v <- sum(idx & etio_viral)
  k <- if (viral) v else n - v
  sprintf("%d (%.1f)", k, 100 * k / n)
}
rows[[ri]] <- data.frame(Variable = "Disease etiology, n (%)",
                         Healthy = "", MASLD = "", Cirrhosis = "", HCC = "",
                         P = "", Test = "", stringsAsFactors = FALSE)
ri <- ri + 1
rows[[ri]] <- data.frame(Variable = "  Non-viral",
                         Healthy = "--", MASLD = etio_np("MASLD", FALSE),
                         Cirrhosis = etio_np("Cirrhosis", FALSE),
                         HCC = etio_np("HCC", FALSE),
                         P = "", Test = "", stringsAsFactors = FALSE)
ri <- ri + 1
rows[[ri]] <- data.frame(Variable = "  Viral hepatitis",
                         Healthy = "--", MASLD = etio_np("MASLD", TRUE),
                         Cirrhosis = etio_np("Cirrhosis", TRUE),
                         HCC = etio_np("HCC", TRUE),
                         P = "", Test = "", stringsAsFactors = FALSE)
ri <- ri + 1
for (vt in c("HCV", "HBV")) {
  sub_np <- function(grp) { idx <- meta$Status == grp; k <- sum(idx & etio_viral & grepl(vt, meta$Causes_of_Diseases_2))
                            sprintf("%d (%.1f)", k, 100 * k / sum(idx)) }
  rows[[ri]] <- data.frame(Variable = paste0("    ", vt), Healthy = "--", MASLD = sub_np("MASLD"),
                           Cirrhosis = sub_np("Cirrhosis"), HCC = sub_np("HCC"),
                           P = "", Test = "", stringsAsFactors = FALSE)
  ri <- ri + 1
}

# ── Liver disease severity in cirrhosis and HCC (Reviewer 2, points 18-19) ──
# Child-Pugh class as documented in the medical reports (overall class only);
# ALBI = log10(bilirubin umol/L) x 0.66 + albumin (g/L) x -0.085 (Johnson et al.,
# J Clin Oncol 2015); grade 1 <= -2.60, grade 2 > -2.60 to <= -1.39, grade 3 > -1.39.
# MELD (UNOS 2016-before-Na formula) computable only where INR was recorded (HCC).
# Healthy and MASLD shown as "--".
sev_header <- function(lbl) data.frame(Variable = lbl, Healthy = "", MASLD = "", Cirrhosis = "", HCC = "",
                                       P = "", Test = "", stringsAsFactors = FALSE)
sev_row <- function(lbl, f) data.frame(Variable = lbl, Healthy = "--", MASLD = "--",
                                       Cirrhosis = f("Cirrhosis"), HCC = f("HCC"),
                                       P = "", Test = "", stringsAsFactors = FALSE)
cp <- trimws(as.character(meta$ChildPugh_Class))
cp_np <- function(level) function(grp) {
  idx <- meta$Status == grp; n <- sum(idx)
  k <- if (is.na(level)) sum(idx & (is.na(cp) | cp == "")) else sum(idx & cp %in% level)
  sprintf("%d (%.1f)", k, 100 * k / n)
}
albi <- log10(meta$Bili_T) * 0.66 + meta$Albumin * -0.085
albi_grade <- cut(albi, c(-Inf, -2.60, -1.39, Inf), labels = c("1", "2", "3"))
albi_np <- function(level) function(grp) {
  idx <- meta$Status == grp; k <- sum(idx & albi_grade %in% level)
  sprintf("%d (%.1f)", k, 100 * k / sum(idx))
}
decomp_tx <- with(meta, (Lactulose_Use %in% "Yes") | (Rifaximin_Use %in% "Yes") | (Diuretic %in% "Yes"))
inr_num  <- suppressWarnings(as.numeric(meta$INR))
meld     <- 3.78 * log(pmax(1, meta$Bili_T / 17.1)) + 11.2 * log(pmax(1, inr_num)) +
            9.57 * log(pmin(4, pmax(1, meta$Creatinine / 88.4))) + 6.43
meld[is.na(inr_num)] <- NA
rows[[ri]] <- sev_header("Liver disease severity"); ri <- ri + 1
rows[[ri]] <- sev_header("  Child-Pugh class, n (%)"); ri <- ri + 1
for (lv in c("A", "B", "C")) {
  rows[[ri]] <- sev_row(paste0("    Class ", lv), cp_np(paste("Class", lv))); ri <- ri + 1
}
rows[[ri]] <- sev_row("    Not recorded", cp_np(NA)); ri <- ri + 1
rows[[ri]] <- sev_row("  ALBI score", function(grp) { x <- albi[meta$Status == grp]; x <- x[!is.na(x)]
  sprintf("%.2f (%.2f to %.2f)", median(x), quantile(x, 0.25), quantile(x, 0.75)) }); ri <- ri + 1
rows[[ri]] <- sev_header("  ALBI grade, n (%)"); ri <- ri + 1
for (lv in c("1", "2", "3")) {
  rows[[ri]] <- sev_row(paste0("    Grade ", lv), albi_np(lv)); ri <- ri + 1
}
rows[[ri]] <- sev_row("  MELD score", function(grp) {
  x <- meld[meta$Status == grp]
  if (all(is.na(x))) "NA (INR not recorded)" else sprintf("%s; n = %d", fmt_med_iqr(x[!is.na(x)]), sum(!is.na(x))) }); ri <- ri + 1
rows[[ri]] <- sev_row("    MELD >= 15", function(grp) {
  x <- meld[meta$Status == grp]
  if (all(is.na(x))) "--" else sprintf("%d (%.1f)", sum(x >= 15, na.rm = TRUE), 100 * mean(x >= 15, na.rm = TRUE)) }); ri <- ri + 1
rows[[ri]] <- sev_row("  Lactulose, rifaximin or diuretic, n (%)", function(grp) {
  idx <- meta$Status == grp; k <- sum(idx & decomp_tx); sprintf("%d (%.1f)", k, 100 * k / sum(idx)) }); ri <- ri + 1
afp_num <- suppressWarnings(as.numeric(meta$AFP))
rows[[ri]] <- sev_header("  AFP category, n (%)"); ri <- ri + 1
for (ab in list(c("    <= 20 ng/mL", -Inf, 20), c("    20-400 ng/mL", 20, 400), c("    > 400 ng/mL", 400, Inf))) {
  lo <- as.numeric(ab[2]); hi <- as.numeric(ab[3])
  rows[[ri]] <- sev_row(ab[1], function(grp) {
    x <- afp_num[meta$Status == grp]; x <- x[!is.na(x)]
    if (!length(x)) return("--")
    k <- sum(x > lo & x <= hi); sprintf("%d/%d (%.1f)", k, length(x), 100 * k / length(x)) }); ri <- ri + 1
}
rows[[ri]] <- sev_row("  Before HCC treatment, n (%)", function(grp)
  if (grp == "HCC") sprintf("%d (100.0)", sum(meta$Status == "HCC")) else "--"); ri <- ri + 1
rows[[ri]] <- sev_row("  BCLC stage / tumor burden", function(grp) if (grp == "HCC") "Not recorded" else "--"); ri <- ri + 1

# ── MASLD characterization (Reviewer 2, point 22): 2023 multisociety cardiometabolic criteria,
#    FibroScan steatosis (CAP) and fibrosis (LSM) grades. Shown for the MASLD column only.
yes <- function(v) tolower(as.character(v)) %in% c("yes", "1", "true")
num <- function(v) suppressWarnings(as.numeric(v))
cm <- data.frame(
  bmi  = num(meta$BMI) >= 25,
  glu  = num(meta$Hgb_A1C) >= 5.7 | num(meta$Glu_R) >= 7.8 | yes(meta$Metformin) | yes(meta$Insulin) |
         yes(meta$Sulfonylurea_Use) | yes(meta$DPP4_Inhibitor_Use) | yes(meta$SGLT2_Inhibitor_Use) | yes(meta$GLP1_RA_Use),
  bp   = yes(meta$antihypertensive),
  tg   = num(meta$TG) >= 1.70 | yes(meta$Statin_Use),
  hdl  = (meta$Gender == "M" & num(meta$HDL) <= 1.0) | (meta$Gender == "F" & num(meta$HDL) <= 1.3))
cm[is.na(cm)] <- FALSE
cm_n <- rowSums(cm)
masld_row <- function(lbl, f) data.frame(Variable = lbl, Healthy = "--", MASLD = f(meta$Status == "MASLD"),
                                         Cirrhosis = "--", HCC = "--", P = "", Test = "", stringsAsFactors = FALSE)
np_str <- function(k, n) sprintf("%d (%.1f)", k, 100 * k / n)
rows[[ri]] <- sev_header("MASLD characterization"); ri <- ri + 1
rows[[ri]] <- sev_header("  Cardiometabolic criteria met, n (%)"); ri <- ri + 1
for (cr in list(c("bmi", "    BMI >= 25 kg/m2"), c("glu", "    Prediabetes/diabetes or glucose-lowering drug"),
                c("bp", "    Antihypertensive treatment"), c("tg", "    TG >= 1.70 mmol/L or lipid-lowering drug"),
                c("hdl", "    Low HDL (<= 1.0 M; <= 1.3 F mmol/L)"))) {
  rows[[ri]] <- masld_row(cr[2], function(i) np_str(sum(cm[[cr[1]]][i]), sum(i))); ri <- ri + 1
}
rows[[ri]] <- masld_row("    >= 2 criteria", function(i) np_str(sum(cm_n[i] >= 2), sum(i))); ri <- ri + 1
st <- trimws(as.character(meta$Steatosis)); fb <- trimws(as.character(meta$Fibrosis))
rows[[ri]] <- sev_header("  Steatosis grade (FibroScan CAP), n (%)"); ri <- ri + 1
for (lv in c("S1", "S2", "S3")) { rows[[ri]] <- masld_row(paste0("    ", lv), function(i) np_str(sum(st[i] %in% lv), sum(i))); ri <- ri + 1 }
rows[[ri]] <- masld_row("    Not recorded", function(i) np_str(sum(i & (is.na(st) | st == "")), sum(i))); ri <- ri + 1
fb_grp <- ifelse(fb %in% c("F0", "F1", "F0-F1"), "F0-F1", ifelse(fb %in% c("F2", "F2-3"), "F2", fb))
rows[[ri]] <- sev_header("  Fibrosis stage (FibroScan LSM), n (%)"); ri <- ri + 1
for (lv in c("F0-F1", "F2", "F3", "F4")) { rows[[ri]] <- masld_row(paste0("    ", lv), function(i) np_str(sum(fb_grp[i] %in% lv), sum(i))); ri <- ri + 1 }
rows[[ri]] <- masld_row("    Not recorded", function(i) np_str(sum(i & (is.na(fb) | fb == "")), sum(i))); ri <- ri + 1
rows[[ri]] <- masld_row("    Significant fibrosis (>= F2), of graded", function(i) {
  g <- i & !(is.na(fb) | fb == ""); np_str(sum(fb_grp[g] %in% c("F2", "F3", "F4")), sum(g)) }); ri <- ri + 1

# Continuous variables
for (var in names(var_labs)) {
  if (!var %in% colnames(meta)) next
  vals <- meta[[var]]
  p <- tryCatch(kruskal.test(vals ~ meta$Status)$p.value,
                error = function(e) NA)
  rows[[ri]] <- data.frame(
    Variable = var_labs[var],
    Healthy = fmt_med_iqr(vals[meta$Status == "Healthy"]),
    MASLD = fmt_med_iqr(vals[meta$Status == "MASLD"]),
    Cirrhosis = fmt_med_iqr(vals[meta$Status == "Cirrhosis"]),
    HCC = fmt_med_iqr(vals[meta$Status == "HCC"]),
    P = paste0(fmt_p(p), stars(p)),
    Test = "KW", stringsAsFactors = FALSE)
  ri <- ri + 1
}

# Medications header
rows[[ri]] <- data.frame(Variable = "Medications, n (%)", Healthy = "",
                          MASLD = "", Cirrhosis = "", HCC = "",
                          P = "", Test = "", stringsAsFactors = FALSE)
ri <- ri + 1

# Medication variables
for (col in names(med_vars)) {
  if (!col %in% colnames(meta)) next
  # Reviewer 2 comment 30: no test when the exposure has a single level (e.g. antibiotics 0% in all groups)
  if (length(unique(meta[[col]][!is.na(meta[[col]]) & meta[[col]] != ""])) <= 1) {
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
    P = paste0(fmt_p(p), stars(p)),
    Test = ifelse(is.na(p), "--", "Chi2"), stringsAsFactors = FALSE)
  ri <- ri + 1
}

st1 <- bind_rows(rows)
cat(sprintf("  %d rows\n", nrow(st1)))

# Build flextable
ft_st1 <- flextable::flextable(st1) %>%
  flextable::set_header_labels(Variable = "Variable",
                    Healthy = "Healthy\n(n = 70)",
                    MASLD = "MASLD\n(n = 59)",
                    Cirrhosis = "Cirrhosis\n(n = 59)",
                    HCC = "HCC\n(n = 17)",
                    P = "P-value", Test = "Test") %>%
  flextable::fontsize(size = 8, part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::bg(part = "header", bg = "#2C3E50") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::bold(part = "header") %>%
  flextable::align(j = 2:7, align = "center", part = "all") %>%
  flextable::align(j = 1, align = "left", part = "body") %>%
  flextable::width(j = 1, width = 2.4) %>%
  flextable::width(j = 2:5, width = 1.6) %>%
  flextable::width(j = 6, width = 0.9) %>%
  flextable::width(j = 7, width = 0.5) %>%
  flextable::border_outer(border = officer::fp_border(color = "#999999")) %>%
  flextable::border_inner_h(border = officer::fp_border(color = "#DDDDDD", width = 0.5)) %>%
  flextable::padding(padding = 2, part = "body")

# Bold section headers and significant p-values
ft_st1 <- flextable::bold(ft_st1, i = which(st1$Variable %in% c("n", "Disease etiology, n (%)", "Liver disease severity", "Medications, n (%)")), j = 1)
for (i in seq_len(nrow(st1))) {
  if (grepl("\\*", st1$P[i])) ft_st1 <- flextable::bold(ft_st1, i = i, j = 6)
}


###############################################################################
#   ST2: ALPHA DIVERSITY PAIRWISE COMPARISONS                                 #
###############################################################################

cat("\n-- ST2: Alpha Diversity --\n")

# Load all 5 metrics from Excel files (Raw_Data sheet)
# Each Excel has: id, Status, <metric_value> in 3rd column
alpha_list <- list()

read_alpha_excel <- function(filepath, metric_name) {
  # Reads the QIIME 2 alpha export (id, metadata columns, metric value in the LAST column)
  cat(sprintf("  Loading %s...\n", metric_name))
  df <- read.delim(filepath, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
  df <- df[df[[1]] != "#q2:types", ]
  df$Status <- gsub("^NAFLD$", "MASLD", df$Status)
  data.frame(SampleID = df[[1]], Status = df$Status,
             value = as.numeric(df[[ncol(df)]]), stringsAsFactors = FALSE)
}

alpha_list[["Shannon"]]           <- read_alpha_excel(SHANNON_FILE, "Shannon")
alpha_list[["Chao1"]]             <- read_alpha_excel(CHAO1_FILE, "Chao1")
alpha_list[["Observed Features"]] <- read_alpha_excel(OBSERVED_FILE, "Observed Features")
alpha_list[["Simpson"]]           <- read_alpha_excel(SIMPSON_FILE, "Simpson")
alpha_list[["Faith PD"]]          <- read_alpha_excel(FAITHPD_FILE, "Faith PD")

# Panel A: Descriptive + Kruskal-Wallis
desc_rows <- list()
di <- 1

for (metric in names(alpha_list)) {
  df <- alpha_list[[metric]]
  df$Status <- factor(df$Status, levels = GROUP_ORDER)
  kw <- kruskal.test(value ~ Status, data = df)
  n_total <- nrow(df)
  k <- length(GROUP_ORDER)
  eps2 <- epsilon_sq(kw$statistic, n_total, k)
  
  for (g in GROUP_ORDER) {
    gv <- df$value[df$Status == g]
    desc_rows[[di]] <- data.frame(
      Metric = metric, Group = g, n = length(gv),
      Median = sprintf("%.3f", median(gv)),
      IQR = sprintf("%.3f", IQR(gv)),
      Range = sprintf("%.3f - %.3f", min(gv), max(gv)),
      stringsAsFactors = FALSE)
    di <- di + 1
  }
  desc_rows[[di]] <- data.frame(
    Metric = metric, Group = "Kruskal-Wallis", n = NA,
    Median = sprintf("H = %.3f", kw$statistic),
    IQR = paste0("p = ", fmt_p(kw$p.value), stars(kw$p.value)),
    Range = sprintf("eps2 = %.4f", eps2),
    stringsAsFactors = FALSE)
  di <- di + 1
  cat(sprintf("  %s: H=%.3f, p=%s\n", metric, kw$statistic, fmt_p(kw$p.value)))
}

st2a <- bind_rows(desc_rows)

# Panel B: Pairwise Dunn's test
dunn_rows <- list()
pi <- 1

for (metric in names(alpha_list)) {
  df <- alpha_list[[metric]]
  df$Status <- factor(df$Status, levels = GROUP_ORDER)
  
  dt <- dunn.test::dunn.test(df$value, df$Status, method = "bh",
                              kw = FALSE, table = FALSE)
  
  for (j in seq_along(dt$comparisons)) {
    groups <- trimws(unlist(strsplit(dt$comparisons[j], " - ")))
    dunn_rows[[pi]] <- data.frame(
      Metric = metric, Group1 = groups[1], Group2 = groups[2],
      Z = sprintf("%.3f", dt$Z[j]),
      FDR = sprintf("%.4f", dt$P.adjusted[j]),
      Sig = stars(dt$P.adjusted[j]),
      stringsAsFactors = FALSE)
    if (dunn_rows[[pi]]$Sig == "") dunn_rows[[pi]]$Sig <- "ns"
    pi <- pi + 1
  }
}

st2b <- bind_rows(dunn_rows)

# Reshape pairwise to wide format
st2b_wide <- st2b %>%
  dplyr::mutate(Comparison = paste0(Group1, " vs ", Group2),
         Result = paste0("Z=", Z, ", FDR=", FDR, " ", Sig)) %>%
  dplyr::select(Metric, Comparison, Result) %>%
  pivot_wider(names_from = Metric, values_from = Result)

# Flextables
ft_st2a <- flextable::flextable(st2a) %>%
  flextable::fontsize(size = 8, part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::bg(part = "header", bg = "#2C3E50") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::bold(part = "header") %>%
  flextable::align(j = 2:6, align = "center", part = "all") %>%
  flextable::width(j = 1, width = 1.4) %>%
  flextable::width(j = 2, width = 1.2) %>%
  flextable::width(j = 3, width = 0.5) %>%
  flextable::width(j = 4:6, width = 1.5) %>%
  flextable::border_outer(border = officer::fp_border(color = "#999999")) %>%
  flextable::border_inner_h(border = officer::fp_border(color = "#DDDDDD", width = 0.5)) %>%
  flextable::padding(padding = 2, part = "body") %>%
  flextable::merge_v(j = 1)

kw_idx <- which(st2a$Group == "Kruskal-Wallis")
ft_st2a <- flextable::bold(ft_st2a, i = kw_idx)
ft_st2a <- flextable::bg(ft_st2a, i = kw_idx, bg = "#FFF3CD")

ft_st2b <- flextable::flextable(st2b_wide) %>%
  flextable::fontsize(size = 7, part = "all") %>%
  flextable::fontsize(size = 8, part = "header") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::bg(part = "header", bg = "#2C3E50") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::bold(part = "header") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1, align = "left", part = "body") %>%
  flextable::width(j = 1, width = 1.8) %>%
  flextable::autofit() %>%
  flextable::border_outer(border = officer::fp_border(color = "#999999")) %>%
  flextable::border_inner_h(border = officer::fp_border(color = "#DDDDDD", width = 0.5)) %>%
  flextable::padding(padding = 2, part = "body")

# Bold significant pairwise comparisons
for (i in seq_len(nrow(st2b_wide))) {
  for (j in 2:ncol(st2b_wide)) {
    val <- st2b_wide[[j]][i]
    if (!is.na(val) && grepl("\\*", val)) ft_st2b <- flextable::bold(ft_st2b, i = i, j = j)
  }
}

# ── ST2 Panel C: effect sizes vs Healthy + minimum detectable effect (Alpha_Effect_Sizes.R) ──
find_rev <- function(fn) {
  # search this script's folder, Desktop/Reviewers, then the whole Desktop (newest copy wins)
  for (root in unique(c(getwd(), REVIEWERS_DIR, dirname(PROJ_DIR)))) {
    h <- list.files(root, pattern = paste0("^", fn, "$"), recursive = TRUE, full.names = TRUE)
    if (length(h)) { f <- h[which.max(file.mtime(h))]; cat("  Using", f, "\n"); return(f) }
  }
  stop("Not found: ", gsub("\\\\", "", fn), "\nSearched: ", getwd(), ", ", REVIEWERS_DIR, ", ", dirname(PROJ_DIR),
       "\nRun Alpha_Effect_Sizes.R first, or copy its CSV files into this script's folder.")
}
es  <- read.csv(find_rev("Alpha_effect_sizes_vs_Healthy\\.csv"), check.names = FALSE)
mde <- read.csv(find_rev("Alpha_minimum_detectable_effect\\.csv"), check.names = FALSE)
st2c <- es %>%
  dplyr::mutate(Group = factor(Group, levels = c("MASLD", "Cirrhosis", "HCC")),
                Cell = sprintf("%+.1f%% (%.1f to %.1f); %+.2f", Pct_diff, CI_low, CI_high, Cliffs_delta)) %>%
  dplyr::select(Metric, Group, Cell) %>%
  tidyr::pivot_wider(names_from = Group, values_from = Cell)
mde_row <- as.data.frame(as.list(c(Metric = "MDE at 80% power: Cohen's d (Cliff's delta)",
  setNames(sprintf("d = %.2f (%.2f)", mde$MDE_Cohens_d, mde$MDE_Cliffs_delta), mde$Group))), check.names = FALSE)
st2c <- dplyr::bind_rows(st2c, mde_row[, names(st2c)])


###############################################################################
#   ST3: ANCOM-BC2 FULL RESULTS (from Figure6_ANCOMBC2_Heatmap.R)             #
###############################################################################

cat("\n-- ST3: ANCOM-BC2 Results --\n")

# Shared flextable style for the long result tables
std_ft <- function(df, widths = NULL, size = 7) {
  ft <- flextable::flextable(df) %>%
    flextable::fontsize(size = size, part = "all") %>% flextable::fontsize(size = 8, part = "header") %>%
    flextable::font(fontname = "Arial", part = "all") %>%
    flextable::bg(part = "header", bg = "#2C3E50") %>% flextable::color(part = "header", color = "white") %>%
    flextable::bold(part = "header") %>% flextable::align(j = 2:ncol(df), align = "center", part = "all") %>%
    flextable::border_outer(border = officer::fp_border(color = "#999999")) %>%
    flextable::border_inner_h(border = officer::fp_border(color = "#DDDDDD", width = 0.5)) %>%
    flextable::padding(padding = 1, part = "body")
  if (!is.null(widths)) for (j in seq_along(widths)) ft <- flextable::width(ft, j = j, width = widths[j])
  ft
}

# ST2 Panel C table (built here because it uses std_ft)
ft_st2c <- std_ft(st2c, widths = c(2.4, 2.3, 2.3, 2.3)) %>%
  flextable::set_header_labels(MASLD = "MASLD vs Healthy", Cirrhosis = "Cirrhosis vs Healthy", HCC = "HCC vs Healthy") %>%
  flextable::italic(i = nrow(st2c), part = "body") %>%
  flextable::bg(i = nrow(st2c), bg = "#F2F2F2", part = "body")
fmt_e <- function(x) ifelse(is.na(x), "--", formatC(x, format = "e", digits = 2))
clean_display <- function(x) gsub("_", " ", gsub("\\[|\\]", "", x))

anc_all <- read.csv(ANCOMBC_FILE, stringsAsFactors = FALSE)
anc_un  <- anc_all %>% dplyr::filter(Model == "Unadjusted")
n_genera_tested <- length(unique(anc_un$Taxon))
st3_tables <- list()
sheet_names_ancom <- c("MASLD vs Healthy", "Cirrhosis vs Healthy", "HCC vs Healthy")

for (sn in sheet_names_ancom) {
  df <- anc_un %>% dplyr::filter(Contrast == sn) %>%
    dplyr::mutate(Direction = case_when(Significant & LFC > 0 ~ "Enriched",
                                 Significant & LFC < 0 ~ "Depleted", TRUE ~ "ns")) %>%
    dplyr::arrange(Q) %>%
    dplyr::filter(Q < 0.05)   # Word shows q < 0.05 only; all genera are in the Excel workbook (sheet S4A)
  cat(sprintf("  %s: %d genera, %d significant\n", sn, nrow(df), sum(df$Direction != "ns")))
  df_fmt <- df %>% transmute(Genus = clean_display(Genus), Log2FC = sprintf("%.3f", LFC), SE = sprintf("%.3f", SE),
                             W = sprintf("%.3f", W), P_value = fmt_e(P), FDR_Q = fmt_e(Q),
                             Sensitivity = ifelse(is.na(Passed_SS), "--", ifelse(Passed_SS, "Pass", "Fail")),
                             Direction)
  ft <- std_ft(df_fmt, widths = c(2.0, 0.8, 0.7, 0.8, 1.0, 1.0, 0.9, 0.9)) %>%
    flextable::set_header_labels(Log2FC = "LFC (ln)", W = "W-statistic", P_value = "P-value", FDR_Q = "FDR q-value",
                      Sensitivity = "Pseudo-count\nsensitivity") %>%
    flextable::italic(j = 1, part = "body")
  for (i in seq_len(nrow(df))) {
    ft <- flextable::color(ft, i = i, j = 2, color = ifelse(df$LFC[i] > 0, "#C0392B", "#2980B9"))
    if (df$Direction[i] != "ns") {
      ft <- flextable::color(ft, i = i, j = 8, color = ifelse(df$Direction[i] == "Enriched", "#C0392B", "#2980B9"))
      ft <- flextable::bold(ft, i = i, j = c(6, 8))
    } else ft <- flextable::color(ft, i = i, j = 8, color = "#888888")
  }
  st3_tables[[sn]] <- ft
}


###############################################################################
#   ST4: ANCOM-BC2 / DESeq2 CONCORDANCE (from Supp Fig S2 script)             #
###############################################################################

cat("\n-- ST4: Cross-method concordance --\n")

st4_raw <- read_excel(DESEQ2_FILE, sheet = "Concordance_Details")
lfc_sig <- function(lfc, q, sig) ifelse(is.na(lfc), "--",
  paste0(sprintf("%.2f", lfc), ifelse(sig %in% TRUE, sapply(q, stars), "")))
st4_display <- st4_raw %>% transmute(
  Genus = clean_display(Genus),
  A_MASLD = lfc_sig(ANC_LFC_MASLD, ANC_Q_MASLD, ANC_Sig_MASLD),
  A_Cirr  = lfc_sig(ANC_LFC_Cirrhosis, ANC_Q_Cirrhosis, ANC_Sig_Cirrhosis),
  A_HCC   = lfc_sig(ANC_LFC_HCC, ANC_Q_HCC, ANC_Sig_HCC),
  D_MASLD = lfc_sig(DES_LFC_MASLD, DES_Q_MASLD, DES_Sig_MASLD),
  D_Cirr  = lfc_sig(DES_LFC_Cirrhosis, DES_Q_Cirrhosis, DES_Sig_Cirrhosis),
  D_HCC   = lfc_sig(DES_LFC_HCC, DES_Q_HCC, DES_Sig_HCC),
  C_MASLD = ifelse(Concordant_MASLD %in% TRUE, "Yes", "No"),
  C_Cirr  = ifelse(Concordant_Cirrhosis %in% TRUE, "Yes", "No"),
  C_HCC   = ifelse(Concordant_HCC %in% TRUE, "Yes", "No"))
n_conc <- sum(st4_raw$Contrasts_concordant >= 1)
cat(sprintf("  %d ANCOM-BC2 genera; %d concordant in >= 1 contrast\n", nrow(st4_display), n_conc))

ft_st4 <- std_ft(st4_display, widths = c(1.9, rep(0.75, 6), rep(0.72, 3))) %>%
  flextable::set_header_labels(A_MASLD = "MASLD", A_Cirr = "Cirrhosis", A_HCC = "HCC",
                    D_MASLD = "MASLD", D_Cirr = "Cirrhosis", D_HCC = "HCC",
                    C_MASLD = "MASLD", C_Cirr = "Cirrhosis", C_HCC = "HCC") %>%
  flextable::add_header_row(values = c("", "ANCOM-BC2 LFC (ln)", "DESeq2 log2FC", "Concordant"),
                 colwidths = c(1, 3, 3, 3), top = TRUE) %>%
  flextable::bg(i = 1, part = "header", bg = "#2C3E50") %>% flextable::color(part = "header", color = "white") %>%
  flextable::italic(j = 1, part = "body")
for (i in seq_len(nrow(st4_display))) for (j in 8:10)
  if (st4_display[[j]][i] == "Yes") ft_st4 <- flextable::bold(flextable::color(ft_st4, i = i, j = j, color = "#27AE60"), i = i, j = j)


###############################################################################
#   ST5: SPEARMAN CORRELATIONS, POOLED AND GROUP-ADJUSTED (Figure 7 script)   #
###############################################################################

cat("\n-- ST5: Spearman Correlations --\n")

sp <- read.csv(SPEARMAN_FILE, stringsAsFactors = FALSE)
param_order <- c("ALT", "AST", "GGT", "ALP", "Albumin", "Bilirubin", "Platelet", "BMI", "HbA1c", "TG", "HDL")
mk_wide <- function(rho, q) {
  sp %>% dplyr::mutate(display = paste0(sprintf("%.2f", .data[[rho]]),
                                 ifelse(.data[[q]] < 0.001, "***", ifelse(.data[[q]] < 0.01, "**",
                                 ifelse(.data[[q]] < 0.05, "*", ""))))) %>%
    dplyr::select(Genus, Clinical_Parameter, display) %>%
    pivot_wider(names_from = Clinical_Parameter, values_from = display) %>%
    dplyr::select(Genus, any_of(param_order))
}
st5a_df <- mk_wide("Spearman_Rho", "FDR_Q")
st5b_df <- mk_wide("Partial_Spearman_Rho", "Partial_FDR_Q")
cat(sprintf("  %d genera x %d parameters; pooled significant %d, partial significant %d\n",
            nrow(st5a_df), ncol(st5a_df) - 1, sum(sp$FDR_Q < 0.05, na.rm = TRUE),
            sum(sp$Partial_FDR_Q < 0.05, na.rm = TRUE)))
color_rho <- function(ft, df) {
  for (i in seq_len(nrow(df))) for (j in 2:ncol(df)) {
    v <- df[[j]][i]
    if (!is.na(v) && grepl("\\*", v)) {
      r <- as.numeric(sub("\\*+$", "", v))
      ft <- flextable::bold(flextable::color(ft, i = i, j = j, color = ifelse(r > 0, "#C0392B", "#2980B9")), i = i, j = j)
    }
  }
  ft
}
ft_st5a <- color_rho(std_ft(st5a_df, widths = c(1.8, rep(0.6, ncol(st5a_df) - 1))) %>% flextable::italic(j = 1, part = "body"), st5a_df)
ft_st5b <- color_rho(std_ft(st5b_df, widths = c(1.8, rep(0.6, ncol(st5b_df) - 1))) %>% flextable::italic(j = 1, part = "body"), st5b_df)

###############################################################################
#   ST6: PERMANOVA VARIANCE PARTITIONING                                      #
###############################################################################

cat("\n-- ST6: PERMANOVA (computed in-script from distance matrices) --\n")

# Rebuild a square distance matrix from a QIIME2 long-format pairwise export
read_long_dm <- function(path) {
  d <- read.delim(path, check.names = FALSE)
  ids <- sort(unique(c(d$SubjectID1, d$SubjectID2)))
  M <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  M[cbind(match(d$SubjectID1, ids), match(d$SubjectID2, ids))] <- d$Distance
  M[cbind(match(d$SubjectID2, ids), match(d$SubjectID1, ids))] <- d$Distance
  as.dist(M)
}
dm_list <- list(`Bray-Curtis`        = read_long_dm(BRAY_TSV),
                `Unweighted UniFrac` = read_long_dm(UUF_TSV),
                `Weighted UniFrac`   = read_long_dm(WUF_TSV))

# Metadata + the SAME imputation used by the adjusted model:
# 3 HCC participants missing LDL/HDL/TG (cohort median) and metformin (mode).
meta_p <- read.delim(METADATA_FILE, check.names = FALSE)
rownames(meta_p) <- meta_p$`#SampleID`
for (cc in c("Age", "BMI", "LDL", "TG", "HDL"))
  meta_p[[cc]] <- as.numeric(as.character(meta_p[[cc]]))
for (cc in c("LDL", "TG", "HDL"))
  meta_p[[cc]][is.na(meta_p[[cc]])] <- median(meta_p[[cc]], na.rm = TRUE)
.mode <- function(x) { x <- x[!is.na(x) & x != ""]; names(sort(table(x), decreasing = TRUE))[1] }
meta_p$Metformin[is.na(meta_p$Metformin) | meta_p$Metformin == ""] <- .mode(meta_p$Metformin)
meta_p$Status <- factor(meta_p$Status)
meta_p <- meta_p[rownames(meta_p) %in% meta$SampleID, ]          # same exclusions as everywhere else
dm_list <- lapply(dm_list, function(D) { m <- as.matrix(D); k <- intersect(rownames(m), rownames(meta_p)); as.dist(m[k, k]) })

# Variable label -> metadata column (display order is preserved in the table)
.perm_vars <- c("Disease Status" = "Status", "Sex" = "Gender", "Age" = "Age",
                "BMI" = "BMI", "PPI Use" = "PPI_Use", "Metformin" = "Metformin",
                "Rifaximin" = "Rifaximin_Use", "Lactulose" = "Lactulose_Use",
                "Statin Use" = "Statin_Use", "TG" = "TG", "LDL" = "LDL", "HDL" = "HDL")

# Univariate PERMANOVA per metric; BH across the 12 variables within each metric
set.seed(2024)
perm_uni <- do.call(rbind, lapply(names(dm_list), function(mn) {
  D  <- dm_list[[mn]]
  md <- meta_p[labels(D), , drop = FALSE]            # align metadata to matrix order
  res <- do.call(rbind, lapply(names(.perm_vars), function(lab) {
    a <- adonis2(as.formula(paste0("D ~ `", .perm_vars[[lab]], "`")),
                 data = md, permutations = 999, by = "terms")
    data.frame(Metric = mn, Variable = lab,
               R2 = a$R2[1], F_statistic = a$F[1], FDR = a$`Pr(>F)`[1])
  }))
  res$FDR <- p.adjust(res$FDR, "BH")
  res
}))
perm_uni$Variable <- factor(perm_uni$Variable, levels = names(.perm_vars))
perm_uni <- perm_uni[order(match(perm_uni$Metric, names(dm_list)), perm_uni$Variable), ]

# Provenance guard: weighted-UniFrac disease MUST equal Table 2 / Figure 3.
.chk <- subset(perm_uni, Metric == "Weighted UniFrac" & Variable == "Disease Status")
cat(sprintf("  %d rows computed; weighted UniFrac disease R2=%.4f, F=%.3f\n",
            nrow(perm_uni), .chk$R2, .chk$F_statistic))

st6_df <- perm_uni %>%
  dplyr::mutate(R2_fmt = sprintf("%.4f", R2),
         F_fmt = sprintf("%.3f", F_statistic),
         FDR_fmt = paste0(sapply(FDR, fmt_p), sapply(FDR, stars))) %>%
  dplyr::select(Metric, Variable, R2_fmt, F_fmt, FDR_fmt) %>%
  pivot_wider(names_from = Metric,
              values_from = c(R2_fmt, F_fmt, FDR_fmt),
              names_glue = "{Metric}_{.value}")

# Reorder columns logically: Variable, then per metric R2/F/FDR
metrics <- unique(perm_uni$Metric)
col_order <- "Variable"
for (m in metrics) {
  col_order <- c(col_order,
                 paste0(m, "_R2_fmt"),
                 paste0(m, "_F_fmt"),
                 paste0(m, "_FDR_fmt"))
}
col_order <- col_order[col_order %in% colnames(st6_df)]
st6_df <- st6_df[, col_order]

# ── Clean column names: replace raw R names with formatted labels ──
clean_labels <- setNames(colnames(st6_df), colnames(st6_df))
clean_labels["Variable"] <- "Variable"
for (m in metrics) {
  clean_labels[paste0(m, "_R2_fmt")]  <- "R\u00B2"
  clean_labels[paste0(m, "_F_fmt")]   <- "F"
  clean_labels[paste0(m, "_FDR_fmt")] <- "FDR"
}

ft_st6 <- flextable::flextable(st6_df) %>%
  flextable::set_header_labels(values = clean_labels) %>%
  flextable::fontsize(size = 7, part = "all") %>%
  flextable::fontsize(size = 8, part = "header") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::bg(part = "header", bg = "#2C3E50") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::bold(part = "header") %>%
  flextable::align(j = 2:ncol(st6_df), align = "center", part = "all") %>%
  flextable::width(j = 1, width = 1.1) %>%
  flextable::width(j = 2:ncol(st6_df), width = 0.95) %>%
  flextable::border_outer(border = officer::fp_border(color = "#999999")) %>%
  flextable::border_inner_h(border = officer::fp_border(color = "#DDDDDD", width = 0.5)) %>%
  flextable::padding(padding = 1, part = "body")

# Add spanning header for metrics
if (length(metrics) >= 1) {
  header_vals <- c("")
  header_widths <- c(1)
  for (m in metrics) {
    header_vals <- c(header_vals, m)
    header_widths <- c(header_widths, 3)
  }
  ft_st6 <- flextable::add_header_row(ft_st6, values = header_vals,
                            colwidths = header_widths, top = TRUE)
  ft_st6 <- flextable::bg(ft_st6, i = 1, part = "header", bg = "#2C3E50")
  ft_st6 <- flextable::color(ft_st6, i = 1, part = "header", color = "white")
  ft_st6 <- flextable::bold(ft_st6, i = 1, part = "header")
}

# Highlight Disease Status row
ds_row <- which(grepl("Disease|Status", st6_df$Variable, ignore.case = TRUE))
if (length(ds_row) > 0) {
  ft_st6 <- flextable::bold(ft_st6, i = ds_row)
  ft_st6 <- flextable::bg(ft_st6, i = ds_row, bg = "#FEF9E7")
}


###############################################################################
#   ST8: SENSITIVITY ANALYSIS EXCLUDING VIRAL ETIOLOGY (R2 comment 2)         #
###############################################################################

cat("\n-- ST8: Non-viral sensitivity --\n")

eff <- read_excel(REV03_FILE, sheet = "A_Effect_sizes_all_terms") %>% dplyr::filter(Term == "Status")
sens <- read_excel(REV03_FILE, sheet = "D_Sensitivity_subsets")
st8a <- eff %>% dplyr::select(Metric, Model, R2, P) %>%
  dplyr::mutate(Model = ifelse(Model == "Unadjusted", "Full_unadj", "Full_adj")) %>%
  pivot_wider(names_from = Model, values_from = c(R2, P)) %>%
  left_join(sens %>% dplyr::filter(grepl("Non-viral", Subset)) %>%
              dplyr::select(Metric, NV_R2_unadj = R2_unadj, NV_P_unadj = P_unadj, NV_R2_adj = R2_adj,
                     NV_P_adj = P_adj, NV_Welch = Welch_Wd_P), by = "Metric") %>%
  transmute(Metric,
            `Full: R2 (p)`          = sprintf("%.3f (%s)", R2_Full_unadj, sapply(P_Full_unadj, fmt_p)),
            `Non-viral: R2 (p)`     = sprintf("%.3f (%s)", NV_R2_unadj, sapply(NV_P_unadj, fmt_p)),
            `Full adj: R2 (p)`      = sprintf("%.3f (%s)", R2_Full_adj, sapply(P_Full_adj, fmt_p)),
            `Non-viral adj: R2 (p)` = sprintf("%.3f (%s)", NV_R2_adj, sapply(NV_P_adj, fmt_p)),
            `Non-viral Welch W*d p` = sapply(NV_Welch, fmt_p))
nv_groups <- unique(sens$Groups[grepl("Non-viral", sens$Subset)])

da_all <- read_excel(REV02_FILE, sheet = "All_results")
m1_sig <- da_all %>% dplyr::filter(Model == "M1_Full_unadjusted", ANC_sig) %>% dplyr::select(Contrast, Taxon)
nv_cell <- function(lfc, q, sig) ifelse(is.na(lfc), "Not tested",
  paste0(sprintf("%.2f", lfc), " (", fmt_e(q), ")", ifelse(sig %in% TRUE, "*", ifelse(!is.na(q) & q < 0.05, "\u2020", ""))))
nv_mod <- c(M1_Full_unadjusted = "Full", M4u_NonViral_unadjusted = "NV_unadj", M4_NonViral_adj_core = "NV_adj")
st8b <- m1_sig %>%
  dplyr::left_join(da_all %>% dplyr::filter(Model == "M1_Full_unadjusted") %>% dplyr::select(Contrast, Taxon, Genus),
                   by = c("Contrast", "Taxon"))
for (m in names(nv_mod)) {
  d <- da_all %>% dplyr::filter(Model == m) %>%
    dplyr::transmute(Contrast, Taxon, !!nv_mod[[m]] := nv_cell(ANC_lfc, ANC_q, ANC_sig),
                     !!paste0(nv_mod[[m]], "_DES") := ifelse(is.na(DES_q), "--", fmt_e(DES_q)))
  st8b <- dplyr::left_join(st8b, d, by = c("Contrast", "Taxon"))
}
st8b <- st8b %>% dplyr::mutate(dplyr::across(-c(Contrast, Taxon, Genus), ~ ifelse(is.na(.x), "Not tested", .x))) %>%
  dplyr::transmute(Genus = clean_display(Genus), Contrast,
            `Full cohort: LFC (q)` = Full, `Non-viral unadjusted: LFC (q)` = NV_unadj,
            `Non-viral adjusted: LFC (q)` = NV_adj, `DESeq2 non-viral unadjusted q` = NV_unadj_DES) %>%
  dplyr::arrange(factor(Contrast, levels = c("MASLD vs Healthy", "Cirrhosis vs Healthy", "HCC vs Healthy")), Genus)
ft_st8a <- std_ft(st8a, widths = c(1.9, rep(1.25, 4), 1.0))
ft_st8b <- std_ft(st8b, widths = c(1.9, 1.5, 1.5, 1.6, 1.6, 1.3)) %>% flextable::italic(j = 1, part = "body")
cat(sprintf("  %d metrics; %d genus contrasts\n", nrow(st8a), nrow(st8b)))


###############################################################################
#   ST9: AGE, BMI AND SEX (R2 comments 3-5)                                   #
###############################################################################

cat("\n-- ST9: Age, BMI, sex --\n")

s9 <- read_excel(S9_FILE, sheet = "S9_Age_BMI_Sex")
st9 <- s9 %>% transmute(Metric, Variable,
  `Univariate R2 (p)`    = ifelse(is.na(Univariate_R2), "--", sprintf("%.3f (%s)", Univariate_R2, sapply(Univariate_P, fmt_p))),
  `Adjusted R2 (p)`      = sprintf("%.3f (%s)", Adjusted_R2, sapply(Adjusted_P, fmt_p)),
  `Adjusted omega2`      = sprintf("%.4f", Adjusted_Omega2),
  `Patients only R2 (p)` = sprintf("%.3f (%s)", PatientsOnly_R2, sapply(PatientsOnly_P, fmt_p)),
  `Controls only p`      = ifelse(is.na(ControlsOnly_P), "--", sapply(ControlsOnly_P, fmt_p)))
ft_st9 <- std_ft(st9, widths = c(1.8, 1.1, 1.2, 1.2, 0.9, 1.3, 0.9)) %>% flextable::merge_v(j = 1)
ft_st9 <- flextable::bg(ft_st9, i = which(st9$Variable == "Disease status"), bg = "#FEF9E7")


###############################################################################
#   ST10: COVARIATE-ADJUSTED DIFFERENTIAL ABUNDANCE (R2 comment 6)            #
###############################################################################

cat("\n-- ST10: Adjusted differential abundance --\n")

mod_lab <- c(M1_Full_unadjusted = "Unadjusted", M2b_Full_adj_noAge = "Sex, BMI, medications",
             M2_Full_adj_core = "+ Age", M3_Full_adj_full = "All 11 covariates")
st10 <- da_all %>% dplyr::filter(Model %in% names(mod_lab)) %>%
  semi_join(m1_sig, by = c("Contrast", "Taxon")) %>%
  dplyr::mutate(cell = ifelse(is.na(ANC_lfc), "--",
                       paste0(sprintf("%.2f", ANC_lfc), ifelse(ANC_sig %in% TRUE, sapply(ANC_q, stars), ""))),
         Model = factor(mod_lab[Model], levels = mod_lab)) %>%
  dplyr::select(Genus, Contrast, Model, cell) %>% dplyr::arrange(Model) %>%
  pivot_wider(names_from = Model, values_from = cell) %>%
  left_join(da_all %>% dplyr::filter(Model == "M3_Full_adj_full") %>%
              transmute(Contrast, Genus, `DESeq2, all 11: q` = fmt_e(DES_q)), by = c("Contrast", "Genus")) %>%
  dplyr::mutate(Genus = clean_display(Genus)) %>% dplyr::arrange(Contrast, Genus)
ft_st10 <- std_ft(st10, widths = c(1.9, 1.5, rep(1.0, length(mod_lab)), 1.1)) %>%
  flextable::add_header_row(values = c("", "", "ANCOM-BC2 LFC (ln) (significance)", ""),
                 colwidths = c(1, 1, length(mod_lab), 1), top = TRUE) %>%
  flextable::bg(i = 1, part = "header", bg = "#2C3E50") %>% flextable::color(part = "header", color = "white") %>%
  flextable::italic(j = 1, part = "body")
cat(sprintf("  %d genus contrasts x %d models\n", nrow(st10), length(mod_lab)))


###############################################################################
#   ST11: DISPERSION-ROBUST TESTS OF GROUP SEPARATION (R2 comment 9)          #
###############################################################################

cat("\n-- ST11: Dispersion-robust PERMANOVA --\n")

wd  <- read_excel(REV03_FILE, sheet = "B1_Welch_Wd_star")
bal <- read_excel(REV03_FILE, sheet = "B2_Balanced_subsampling")
wcell <- function(w, p) sprintf("%.2f (%s)", w, sapply(p, fmt_p))
st11 <- wd %>%
  dplyr::mutate(P_show = ifelse(Comparison == "Global (4 groups)", P, P_BH_pairwise),
                Comparison = dplyr::recode(Comparison, "Global (4 groups)" = "Global",
                                           "MASLD vs Healthy" = "MASLD", "Cirrhosis vs Healthy" = "Cirrhosis",
                                           "HCC vs Healthy" = "HCC"),
                cell = wcell(Wd, P_show)) %>%
  dplyr::select(Metric, Comparison, cell) %>%
  tidyr::pivot_wider(names_from = Comparison, values_from = cell) %>%
  dplyr::left_join(bal %>% dplyr::transmute(Metric,
                     `Balanced R2, median (IQR)` = sprintf("%.3f (%s)", R2_median, R2_IQR),
                     `Balanced PERMANOVA p<0.05 (%)` = sprintf("%.1f", Pct_P_lt_0.05),
                     `Balanced PERMDISP p<0.05 (%)` = sprintf("%.1f", Dispersion_Pct_P_lt_0.05)), by = "Metric") %>%
  dplyr::rename(`W*d global (p)` = Global, `W*d MASLD vs H (q)` = MASLD,
                `W*d Cirrhosis vs H (q)` = Cirrhosis, `W*d HCC vs H (q)` = HCC)
ft_st11 <- std_ft(st11, widths = c(1.7, rep(1.05, 4), 1.35, 1.0, 1.0))
n_iter  <- unique(bal$Iterations); n_pg <- unique(bal$n_per_group)
cat(sprintf("  %d metrics\n", nrow(st11)))

###############################################################################
#                                                                             #
#   EXPORT TO WORD                                                            #
#                                                                             #
###############################################################################

cat("\n==========================================================\n")

cat("  Exporting to Word\n")
cat("==========================================================\n\n")

doc <- officer::read_docx()

# ── Disable heading auto-numbering in default template ──
# The officer default template links Heading styles to a multilevel
# numbering definition, producing "1.1. Table S1" etc. Override by
# removing the numPr linkage from the heading styles.
doc <- doc %>%
  officer::body_set_default_section(
    officer::prop_section(
      page_size = officer::page_size(width = 11, height = 8.5, orient = "landscape"),
      page_margins = officer::page_mar(top = 0.5, bottom = 0.5,
                              left = 0.7, right = 0.7)
    )
  )

# ── ST1 ──
doc <- doc %>%
  officer::body_add_par("Table S1. Clinical and Demographic Characteristics by Disease Group",
               style = "heading 2") %>%
  officer::body_add_par(paste0(
    "Data presented as median (IQR) for continuous variables and n (%) ",
    "for categorical variables. P-values: Kruskal-Wallis (KW) for continuous, ",
    "Chi-squared with Monte Carlo simulation for categorical. ",
    "*P < 0.05, **P < 0.01, ***P < 0.001. Liver disease severity (cirrhosis and HCC only): Child-Pugh class as ",
    "documented in the medical reports (overall class only); ALBI score = log10(bilirubin, umol/L) x 0.66 + ",
    "albumin (g/L) x -0.085 (grade 1 <= -2.60; grade 2 > -2.60 to <= -1.39; grade 3 > -1.39; Johnson et al. 2015). ",
    "MELD = 3.78 ln(bilirubin, mg/dL) + 11.2 ln(INR) + 9.57 ln(creatinine, mg/dL) + 6.43, calculable only where INR ",
    "was recorded (not recorded in cirrhosis). All HCC samples were collected at diagnosis, before HCC-directed treatment. ",
    "MASLD: 2023 multisociety cardiometabolic criteria (blood pressure criterion assessed as antihypertensive treatment); ",
    "steatosis and fibrosis graded by FibroScan (CAP; liver stiffness, >= F2 = >= 8.2 kPa); no liver biopsies were performed."
  ), style = "Normal") %>%
  flextable::body_add_flextable(ft_st1)

# ── ST2 ──
doc <- doc %>%
  officer::body_add_par("Table S2. Alpha Diversity Indices Across Disease Groups",
               style = "heading 2") %>%
  officer::body_add_par("Panel A: Descriptive Statistics and Global Tests",
               style = "heading 3") %>%
  officer::body_add_par(paste0(
    "Five alpha diversity metrics with group-level descriptive statistics ",
    "and Kruskal-Wallis test. eps2 = epsilon-squared effect size. Samples rarefied to 13,000 reads. ",
    "Chao1 approximates observed features because DADA2 removes singleton sequence variants (no sample contained a ",
    "singleton ASV before rarefaction; after rarefaction, Chao1 exceeded observed features by a median of 1.0%)."
  ), style = "Normal") %>%
  flextable::body_add_flextable(ft_st2a) %>%
  officer::body_add_par("Panel B: Pairwise Comparisons (Dunn's Test, BH-corrected)",
               style = "heading 3") %>%
  officer::body_add_par(paste0("Pairwise tests are shown for completeness only; no global Kruskal-Wallis test ",
                               "was significant after Benjamini-Hochberg correction across the five metrics ",
                               "(all q >= 0.089), so pairwise differences are not interpreted. ",
                               "*FDR < 0.05, **FDR < 0.01, ***FDR < 0.001."),
               style = "Normal") %>%
  flextable::body_add_flextable(ft_st2b) %>%
  officer::body_add_par("Panel C: Effect Sizes Relative to Healthy Controls", style = "heading 3") %>%
  officer::body_add_par(paste0("Percentage difference in medians (95% CI from 2,000 bootstrap resamples); Cliff's delta ",
                               "(positive = higher in the disease group). Samples rarefied to 13,000 reads. Last row: minimum ",
                               "detectable effect at 80% power (two-sided alpha = 0.05) for each comparison with healthy controls (n = 70)."),
               style = "Normal") %>%
  flextable::body_add_flextable(ft_st2c)

# ── S3: COMMUNITY-LEVEL ROBUSTNESS (old S11, S6, S9, S8A) ──
doc <- doc %>%
  officer::body_add_par("Table S3. Robustness of Community-Level Differences (PERMANOVA)", style = "heading 2") %>%
  officer::body_add_par("Panel A: Dispersion-robust tests of group separation", style = "heading 3") %>%
  officer::body_add_par(paste0("W*d = Welch-type heteroscedastic PERMANOVA statistic (Hamidi et al. 2019), which does not ",
                      "assume equal multivariate dispersion; 999 permutations; pairwise comparisons vs Healthy BH-adjusted. ",
                      "Balanced = PERMANOVA and PERMDISP repeated in ", n_iter, " random subsamples of ", n_pg,
                      " participants per group (PERMANOVA is robust to unequal dispersion when group sizes are equal; ",
                      "Anderson & Walsh 2013); 199 permutations each."), style = "Normal") %>%
  flextable::body_add_flextable(ft_st11) %>%
  officer::body_add_par("Panel B: Univariate PERMANOVA variance partitioning", style = "heading 3") %>%
  officer::body_add_par(paste0(
    "Univariate PERMANOVA (adonis2, 999 permutations) of each beta-diversity distance matrix on disease status and on ",
    "each candidate confounder, fitted one variable at a time on the full n = 205 cohort (isolated missing LDL/HDL/TG ",
    "and metformin values in three HCC participants imputed by cohort median/mode, as in the adjusted model). FDR by ",
    "Benjamini-Hochberg across the 12 variables within each metric. Disease Status row highlighted. ",
    "*FDR < 0.05, **FDR < 0.01, ***FDR < 0.001."), style = "Normal") %>%
  flextable::body_add_flextable(ft_st6) %>%
  officer::body_add_par("Panel C: Age, BMI and sex", style = "heading 3") %>%
  officer::body_add_par(paste0("Univariate = variable alone. Adjusted = marginal effect in the model with disease status and ",
                      "11 covariates (all participants). Omega-squared is comparable across terms with different degrees ",
                      "of freedom. Patients only = same model in MASLD, cirrhosis and HCC (ages overlap). Controls only = ",
                      "age, BMI and sex within healthy controls. 999 permutations."), style = "Normal") %>%
  flextable::body_add_flextable(ft_st9) %>%
  officer::body_add_par("Panel D: Sensitivity analysis excluding participants with viral hepatitis", style = "heading 3") %>%
  officer::body_add_par(paste0("Non-viral subset: ", nv_groups, ". Adjusted = marginal PERMANOVA with 11 covariates ",
                      "(999 permutations). Welch W*d = heteroscedasticity-robust PERMANOVA (Hamidi et al. 2019)."),
               style = "Normal") %>%
  flextable::body_add_flextable(ft_st8a)

# ── S4: DIFFERENTIAL ABUNDANCE (old S3, S10, S4, S8B) ──
doc <- doc %>%
  officer::body_add_par("Table S4. Differential Abundance Analyses", style = "heading 2") %>%
  officer::body_add_par(paste0("Full results for every genus are provided in the accompanying Excel workbook ",
                      "(Supplementary_Table_S4_Differential_Abundance.xlsx; sheets S4A-S4D)."), style = "Normal") %>%
  officer::body_add_par(sprintf("Panel A: ANCOM-BC2 (%d genera tested; q < 0.05 shown)", n_genera_tested), style = "heading 3") %>%
  officer::body_add_par(paste0(
    "ANCOM-BC2 (ANCOMBC R package; Status ~ Healthy reference; prevalence filter 10%; ",
    "structural-zero detection; pseudo-count sensitivity analysis). LFC = natural-log fold change relative to Healthy; ",
    "positive = enriched in disease. FDR = Benjamini-Hochberg across genera within each contrast. ",
    "Genera with q < 0.05 are shown; Direction is given only when the genus also passed the pseudo-count sensitivity ",
    "analysis (q < 0.05 with Direction = ns indicates a non-robust result). All genera are listed in sheet S4A of the Excel workbook. ",
    "Red = enriched, blue = depleted."), style = "Normal")
for (sn in sheet_names_ancom) {
  doc <- doc %>% officer::body_add_par(sn, style = "heading 3") %>% flextable::body_add_flextable(st3_tables[[sn]])
}
doc <- doc %>%
  officer::body_add_par("Panel B: Covariate-adjusted differential abundance", style = "heading 3") %>%
  officer::body_add_par(paste0("ANCOM-BC2 LFC (ln) vs Healthy for genera significant in the unadjusted model, across ",
                      "increasingly adjusted models (continuous covariates z-scored). Stars = q < 0.05 and passing the ",
                      "sensitivity analysis. Age is near-collinear with control status (controls 22-39 years), which ",
                      "inflates standard errors once it is included. Last column: DESeq2 q-value in the 11-covariate model."),
               style = "Normal") %>%
  flextable::body_add_flextable(ft_st10) %>%
  officer::body_add_par("Panel C: Cross-method concordance of ANCOM-BC2 and DESeq2", style = "heading 3") %>%
  officer::body_add_par(sprintf(paste0(
    "Fold changes for the %d genera significant in ANCOM-BC2 in >= 1 contrast (ANCOM-BC2: natural log; DESeq2: log2). DESeq2 was fitted ",
    "as one model including all groups (Healthy reference; poscounts size factors). Concordant = q < 0.05 ",
    "in both methods with the same direction (%d genera in >= 1 contrast). DESeq2 is a complementary method ",
    "applied to the same data, not independent validation. *q<0.05, **q<0.01, ***q<0.001 (Benjamini-Hochberg)."),
    nrow(st4_display), n_conc), style = "Normal") %>%
  flextable::body_add_flextable(ft_st4) %>%
  officer::body_add_par("Panel D: Sensitivity analysis excluding participants with viral hepatitis", style = "heading 3") %>%
  officer::body_add_par(paste0("ANCOM-BC2 LFC (ln) (BH q) for genera significant in the full cohort, refitted after excluding ",
                      "viral cases (", nv_groups, "). Adjusted = age, sex, BMI, PPI, metformin and statin use. ",
                      "* = q < 0.05 and passing the pseudo-count sensitivity analysis; \u2020 = q < 0.05 but not robust to ",
                      "the pseudo-count choice. Not tested = below the 10% prevalence filter in the reduced dataset."),
               style = "Normal") %>%
  flextable::body_add_flextable(ft_st8b)

# ── ST5 ──
doc <- doc %>%
  officer::body_add_par("Table S5. Spearman Correlations: Differentially Abundant Genera vs Clinical Parameters",
               style = "heading 2") %>%
  officer::body_add_par("Panel A: Pooled across all participants", style = "heading 3") %>%
  officer::body_add_par(paste0("Spearman rho; BH-FDR across all genus x parameter tests. Pooled correlations include ",
                      "between-group differences. Red = positive, blue = negative. *FDR<0.05, **FDR<0.01, ***FDR<0.001."),
               style = "Normal") %>%
  flextable::body_add_flextable(ft_st5a) %>%
  officer::body_add_par("Panel B: Partial correlations adjusted for disease group", style = "heading 3") %>%
  officer::body_add_par(sprintf("Partial Spearman rho adjusted for disease group (ppcor); BH-FDR. %d correlations significant.",
                       sum(sp$Partial_FDR_Q < 0.05, na.rm = TRUE)), style = "Normal") %>%
  flextable::body_add_flextable(ft_st5b)

# ── S6: STORMS CHECKLIST ──
doc <- doc %>%
  officer::body_add_par("Table S6. STORMS Checklist", style = "heading 2") %>%
  officer::body_add_par(paste0("Strengthening The Organization and Reporting of Microbiome Studies (STORMS; Mirzayi et al. 2021). ",
                      "The completed checklist is provided as a separate file."), style = "Normal")


###############################################################################
#   EXCEL WORKBOOK FOR TABLE S4 (full, unformatted values)                    #
###############################################################################
cat("  Writing", basename(XLSX_S4), "\n")
s4a_x <- anc_un %>% dplyr::mutate(Genus = clean_display(Genus),
  Direction = dplyr::case_when(Significant & LFC > 0 ~ "Enriched", Significant & LFC < 0 ~ "Depleted", TRUE ~ "ns")) %>%
  dplyr::select(Contrast, Genus, dplyr::any_of(c("Taxon")), LFC_ln = LFC, SE, W, P, Q_BH = Q,
                Passed_pseudocount_sensitivity = Passed_SS, Direction) %>%
  dplyr::arrange(factor(Contrast, levels = sheet_names_ancom), Q_BH)
s4c_x <- st4_raw %>% dplyr::mutate(Genus = clean_display(Genus))
s4b_x <- da_all %>% dplyr::filter(Model %in% c(names(mod_lab))) %>%
  dplyr::mutate(Model = unname(mod_lab[Model]), Genus = clean_display(Genus))
s4d_x <- da_all %>% dplyr::filter(Model %in% names(nv_mod)) %>%
  dplyr::mutate(Genus = clean_display(Genus))
wb <- openxlsx::createWorkbook()
hs <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#2C3E50", fontColour = "white")
add_sheet <- function(nm, title, df) {
  openxlsx::addWorksheet(wb, nm)
  openxlsx::writeData(wb, nm, title, startRow = 1)
  openxlsx::writeData(wb, nm, df, startRow = 3, headerStyle = hs)
  openxlsx::freezePane(wb, nm, firstActiveRow = 4)
  openxlsx::setColWidths(wb, nm, cols = seq_len(ncol(df)), widths = "auto")
}
add_sheet("S4A_ANCOMBC2_all", sprintf("Table S4A. ANCOM-BC2, all %d genera, unadjusted (LFC natural log vs Healthy; BH q).", n_genera_tested), s4a_x)
add_sheet("S4B_Adjusted", "Table S4B. ANCOM-BC2 and DESeq2 results across covariate-adjusted models (all genera).", s4b_x)
add_sheet("S4C_DESeq2_concordance", "Table S4C. ANCOM-BC2 (ln) vs DESeq2 (log2) for ANCOM-BC2-significant genera.", s4c_x)
add_sheet("S4D_NonViral", "Table S4D. ANCOM-BC2 and DESeq2 in the full cohort and after excluding viral cases (all genera).", s4d_x)
openxlsx::saveWorkbook(wb, XLSX_S4, overwrite = TRUE)

# Apply landscape to entire document
officer::body_end_section_landscape(doc)

# Save initial version
print(doc, target = OUTPUT_FILE)

# ── Post-process DOCX to fix officer template defaults ──
# officer's default template links Heading styles to multilevel numbering
# (producing "1.1. Table S1") and uses generous cell spacing. Fix via XML.
cat("  Post-processing DOCX for compact layout...\n")

tmp_dir <- tempfile("docx_fix_")
dir.create(tmp_dir)
utils::unzip(OUTPUT_FILE, exdir = tmp_dir)

# Fix 1: Remove heading auto-numbering from styles.xml
styles_path <- file.path(tmp_dir, "word", "styles.xml")
styles_xml  <- readChar(styles_path, file.info(styles_path)$size, useBytes = TRUE)
styles_xml  <- gsub(
  "<w:numPr>\\s*<w:ilvl w:val=\"\\d+\"/>\\s*<w:numId w:val=\"\\d+\"/>\\s*</w:numPr>",
  "", styles_xml
)
writeChar(styles_xml, styles_path, eos = NULL, useBytes = TRUE)

# Fix 2: Compact cell spacing and row heights in document.xml
doc_path <- file.path(tmp_dir, "word", "document.xml")
doc_xml  <- readChar(doc_path, file.info(doc_path)$size, useBytes = TRUE)

# Reduce paragraph spacing inside table cells to zero
doc_xml <- gsub('w:after="100" w:before="100"', 'w:after="0" w:before="0"', doc_xml)
doc_xml <- gsub('w:after="40" w:before="40"',   'w:after="0" w:before="0"', doc_xml)
doc_xml <- gsub('w:after="20" w:before="20"',   'w:after="0" w:before="0"', doc_xml)

# Reduce row heights for compact tables
doc_xml <- gsub(
  '<w:trHeight w:val="360" w:hRule="auto"/>',
  '<w:trHeight w:val="160" w:hRule="auto"/>', doc_xml
)

# Unify first-section margins (officer default = 1417 twips = 1")
# to match our compact settings (720 = 0.5", 1000 = ~0.7")
doc_xml <- gsub(
  '<w:pgMar w:header="708" w:bottom="1417" w:top="1417" w:right="1417" w:left="1417" w:footer="708" w:gutter="0"/>',
  '<w:pgMar w:header="360" w:bottom="720" w:top="720" w:right="1000" w:left="1000" w:footer="360" w:gutter="0"/>',
  doc_xml
)
writeChar(doc_xml, doc_path, eos = NULL, useBytes = TRUE)

# Repack the DOCX using system zip for Word compatibility
# (utils::zip can produce archives Word rejects)
old_wd <- setwd(tmp_dir)
all_files <- list.files(".", recursive = TRUE, all.files = FALSE)
out_path <- normalizePath(OUTPUT_FILE, mustWork = FALSE)
if (file.exists(out_path)) file.remove(out_path)

# DOCX requires [Content_Types].xml at root with no path prefix
# Use system zip with -X (no extra attributes) for maximum compatibility
zip_result <- system2("zip", args = c("-r", "-X", "-q", shQuote(out_path), "."),
                      stdout = TRUE, stderr = TRUE)
if (!file.exists(out_path)) {
  # Fallback: try utils::zip if system zip unavailable
  utils::zip(out_path, files = all_files, flags = "-r9Xq")
}
setwd(old_wd)
unlink(tmp_dir, recursive = TRUE)

cat("  Post-processing complete.\n")

cat("\n==========================================================\n")
cat("  DONE -- saved to:", OUTPUT_FILE, "\n")
cat("  S1: Clinical characteristics\n  S2: Alpha diversity\n")
cat("  S3: PERMANOVA robustness (A dispersion-robust, B univariate, C age/BMI/sex, D non-viral)\n")
cat("  S4: Differential abundance (A ANCOM-BC2, B adjusted, C DESeq2 concordance, D non-viral)\n")
cat("      + Excel:", XLSX_S4, "\n")
cat("  S5: Spearman (pooled + partial)\n  S6: STORMS checklist (separate file)\n")
cat("==========================================================\n")
