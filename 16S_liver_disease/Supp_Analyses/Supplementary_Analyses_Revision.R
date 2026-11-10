################################################################################
# Supplementary_Analyses_Revision.R
# Paper 1 — 16S rRNA Gut Microbiome | Saudi Arabian Liver Disease Cohort
# Microbiology Spectrum revision: ALL supplementary / sensitivity analyses
# ------------------------------------------------------------------------------
# One script, eight sections. Switch sections on/off in RUN below.
# Each section runs in its own environment (no variable leakage between them)
# and writes one Excel workbook to REV_DIR.
#
#   Section               Reviewer points            Output workbook
#   1 BATCH               R2-20                      Rev01_Batch_Effects.xlsx
#   2 PERMANOVA_ROBUSTNESS R2-2,3,5,9,10,11; R3-1    Rev03_PERMANOVA_Robustness.xlsx  (-> Table S8)
#   3 CONFOUNDERS         R2-3,4,5                   Supp_S9_Age_BMI_Sex_PERMANOVA.xlsx (-> Table S9)
#   4 DA_ADJUSTED         R2-2,3,6,7,8,24,25         Rev02_DA_Adjusted_Sensitivity.xlsx (-> Table S10)
#   5 DYSBIOSIS_INDEX     R2-12,13                   Rev04_Dysbiosis_Index_NonCircular.xlsx
#   6 RAREFACTION         R2-26,27                   Rev05_Rarefaction_Sensitivity.xlsx
#   7 CLINICAL            R2-5,18,19,22,30           Rev06_Clinical_Characterisation.xlsx
#   8 HCC_UNCERTAINTY     R2-17                      Rev07_HCC_Uncertainty.xlsx
#
# Main-figure scripts (Figure 6, Figure 7, Supp Fig S2) are separate and write
# to qiime2_analysis/ANCOMBC2_R; Generate_Supplementary_Tables.R reads both.
#
# Duplicate / excluded samples: list SampleIDs (one per line) in
#   <BASE_DIR>/sample_exclusions.txt  -> applied here AND in Figures 6, 7, S2.
################################################################################

RUN <- list(
  BATCH                = FALSE,   # needs batch_map.tsv (run Rev_01_batch_map.sh on the FASTQ folder first)
  PERMANOVA_ROBUSTNESS = TRUE,
  CONFOUNDERS          = TRUE,
  DA_ADJUSTED          = TRUE,    # ~10-20 min
  DYSBIOSIS_INDEX      = FALSE,   # superseded: Figure8_Dysbiosis_Index.R computes the LOO / split-half / partial DI checks
  RAREFACTION          = TRUE,    # needs QIIME2 exports (Step 6d: asv_table.tsv, tree.nwk, stats.tsv); ~30-60 min
  CLINICAL             = TRUE,
  HCC_UNCERTAINTY      = TRUE
)

# ==============================================================================
# 0. CONFIGURATION + HELPERS
# ==============================================================================
suppressMessages({
  library(vegan)
  library(dplyr)
  library(openxlsx)
})
select <- dplyr::select   # ppcor/MASS (section 5) would otherwise mask dplyr::select
filter <- dplyr::filter

# ============================================================================
# 1. PATHS
# ============================================================================

# Paths are RELATIVE: the script finds Desktop/16S_reprocess from any sub-folder
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
BASE_DIR  <- find_proj_dir()                                  # Desktop/16S_reprocess
QIIME_DIR <- file.path(BASE_DIR, "qiime2_analysis_v2")        # reprocessed data (V4, primers removed)

META_FILE         <- file.path(BASE_DIR, "sample_metadata.tsv")
GENUS_COUNTS_FILE <- file.path(QIIME_DIR, "Desq/DeSeq_l6/feature-table_L6.tsv")   # unrarefied genus counts (DESeq2 input)

DIST_FILES <- list(
  "Bray-Curtis"        = file.path(QIIME_DIR, "Beta_diversity/Bray-curtis/Bray-Curtis_raw_data.tsv"),
  "Jaccard"            = file.path(QIIME_DIR, "Beta_diversity/Jaccard/Jaccardraw_data.tsv"),
  "Unweighted UniFrac" = file.path(QIIME_DIR, "Beta_diversity/UnWeighted_UniFrac/Unweighted.tsv"),
  "Weighted UniFrac"   = file.path(QIIME_DIR, "Beta_diversity/Weighted_UniFrac/Weighted_rowdata.tsv")
)

# Exported by the QIIME2 commands in Rev_05 (create these first)
ASV_TABLE_FILE  <- file.path(QIIME_DIR, "revision_exports/asv_table/asv_table.tsv")
TREE_FILE       <- file.path(QIIME_DIR, "revision_exports/tree/tree.nwk")
DADA2_STATS     <- file.path(QIIME_DIR, "revision_exports/dada2_stats/stats.tsv")

# Output folder for ALL revision results
REV_DIR <- file.path(dirname(getwd()), "Supp_Analyses")   # 16S_liver_disease/Supp_Analyses (read by Supp Fig S3 + Tables)
dir.create(REV_DIR, showWarnings = FALSE, recursive = TRUE)

BATCH_MAP_FILE <- file.path(REV_DIR, "batch_map.tsv")   # produced by Rev_01_batch_map.sh
EXCL_FILE      <- file.path(BASE_DIR, "sample_exclusions.txt")   # shared with Figures 6, 7, S2

# ============================================================================
# 2. METADATA COLUMN NAMES
# ============================================================================

ETIOLOGY_COL   <- "Causes_of_Diseases"        # values: NAFLD_Cirrhosis, Cirrhosis_HCC, Viral infection_Cirrhosis/_HCC
VIRAL_PATTERN  <- "Viral"                     # 20 cirrhosis + 5 HCC
CHILD_PUGH_COL <- "ChildPugh_Class"           # "Class A/B/C"; blank = missing
STEATOSIS_COL  <- "Steatosis"
FIBROSIS_COL   <- "Fibrosis"

GROUPS <- c("Healthy", "MASLD", "Cirrhosis", "HCC")
GROUP_COLORS <- c(Healthy = "#4DAF4A", MASLD = "#377EB8", Cirrhosis = "#984EA3", HCC = "#E41A1C")

# Covariate sets
#   CORE  — used for taxon-level (ANCOM-BC2 / DESeq2) adjustment. Rifaximin (n = 5 users)
#           and lactulose (0 users in Healthy and MASLD) are excluded: too sparse /
#           structurally confounded with group to estimate stably at the taxon level.
#   FULL  — the 11 covariates of the published adjusted PERMANOVA (Table 2).
COVARS_CORE <- c("Age_z", "Gender", "BMI_z", "PPI_Use", "Metformin", "Statin_Use")
COVARS_FULL <- c(COVARS_CORE, "Rifaximin_Use", "Lactulose_Use", "TG_z", "LDL_z", "HDL_z")

# Age-overlap sensitivity: controls at or above this age vs MASLD at or below max control age
AGE_OVERLAP_MIN <- 40
MIN_N_SUBSET    <- 10     # minimum per-group n for a subset analysis (non-viral HCC n = 12)

N_PERM <- 999
SEED   <- 2024

# ============================================================================
# 3. HELPERS
# ============================================================================

msg <- function(...) cat(sprintf(...), "\n")

check_cols <- function(df, cols, context = "") {
  miss <- setdiff(cols, names(df))
  if (length(miss)) {
    warning(sprintf("[%s] Column(s) not found: %s\n  Available columns:\n  %s",
                    context, paste(miss, collapse = ", "),
                    paste(names(df), collapse = ", ")), call. = FALSE)
  }
  length(miss) == 0
}

# Load metadata exactly as Figure4_Confounder_PERMANOVA.R does (same imputation),
# plus z-scored continuous covariates and a Viral flag.
load_meta <- function() {
  meta <- read.delim(META_FILE, sep = "\t", stringsAsFactors = FALSE,
                     comment.char = "", check.names = FALSE)
  colnames(meta)[1] <- "SampleID"
  meta <- meta[!grepl("^#", meta$SampleID), ]
  meta$Status <- factor(gsub("NAFLD", "MASLD", meta$Status), levels = GROUPS)
  if (file.exists(EXCL_FILE)) {                      # optional: one SampleID per line
    excl <- trimws(readLines(EXCL_FILE)); excl <- excl[nzchar(excl)]
    meta <- meta[!meta$SampleID %in% excl, ]
    msg("  Excluded %d samples listed in %s", length(excl), basename(EXCL_FILE))
  }

  meta$Gender        <- factor(meta$Gender)
  meta$Age_numeric   <- suppressWarnings(as.numeric(meta$Age))
  meta$BMI_numeric   <- suppressWarnings(as.numeric(meta$BMI))
  meta$TG_numeric    <- suppressWarnings(as.numeric(meta$TG))
  meta$LDL_numeric   <- suppressWarnings(as.numeric(meta$LDL))
  meta$HDL_numeric   <- suppressWarnings(as.numeric(meta$HDL))
  for (v in c("PPI_Use", "Metformin", "Statin_Use", "Rifaximin_Use", "Lactulose_Use"))
    meta[[v]] <- factor(meta[[v]])

  # --- Imputation (identical rule to Figure 4: median / mode) ---
  get_mode <- function(x) { x <- x[!is.na(x) & x != ""]; names(sort(table(x), decreasing = TRUE))[1] }
  for (col in c("Age_numeric", "BMI_numeric", "TG_numeric", "LDL_numeric", "HDL_numeric")) {
    miss <- is.na(meta[[col]])
    if (any(miss)) meta[[col]][miss] <- median(meta[[col]], na.rm = TRUE)
  }
  for (col in c("Gender", "PPI_Use", "Metformin", "Rifaximin_Use", "Lactulose_Use", "Statin_Use")) {
    x <- as.character(meta[[col]]); miss <- is.na(x) | x == ""
    if (any(miss)) { x[miss] <- get_mode(x); meta[[col]] <- factor(x) }
  }

  # z-scored continuous covariates (log2FC per 1 SD; improves model stability)
  meta$Age_z <- as.numeric(scale(meta$Age_numeric))
  meta$BMI_z <- as.numeric(scale(meta$BMI_numeric))
  meta$TG_z  <- as.numeric(scale(meta$TG_numeric))
  meta$LDL_z <- as.numeric(scale(meta$LDL_numeric))
  meta$HDL_z <- as.numeric(scale(meta$HDL_numeric))

  # Viral etiology flag (NA if column missing)
  if (check_cols(meta, ETIOLOGY_COL, "load_meta: etiology")) {
    meta$Viral <- grepl(VIRAL_PATTERN, meta[[ETIOLOGY_COL]], ignore.case = TRUE)
  } else {
    meta$Viral <- NA
  }
  rownames(meta) <- meta$SampleID
  meta
}

# Robust version of the pairwise long-format -> dist converter used in Figures 3/4
pairwise_to_dist <- function(filepath) {
  df  <- read.delim(filepath, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  ids <- sort(unique(c(df$SubjectID1, df$SubjectID2)))
  m   <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  i   <- match(df$SubjectID1, ids); j <- match(df$SubjectID2, ids)
  m[cbind(i, j)] <- df$Distance
  m[cbind(j, i)] <- df$Distance
  as.dist(m)
}

load_dists <- function() lapply(DIST_FILES, pairwise_to_dist)

subset_dist <- function(d, ids) as.dist(as.matrix(d)[ids, ids])

# Genus-level count matrix (taxa x samples), full taxonomy strings as rownames
load_genus_counts <- function(file = GENUS_COUNTS_FILE) {
  raw <- read.table(file, header = TRUE, row.names = 1, sep = "\t", skip = 1,
                    check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
  m <- as.matrix(raw); storage.mode(m) <- "numeric"
  m[rowSums(m) > 0, , drop = FALSE]
}

# Same genus-label cleaner as Figure 6 / Supp S2 scripts
clean_genus <- function(taxon) {
  parts <- strsplit(taxon, ";")[[1]]
  for (offset in 1:7) {
    if (length(parts) >= offset) {
      name <- parts[length(parts) - offset + 1]
      name <- trimws(gsub("^[a-z]__", "", name))
      name <- gsub("_+", "_", name)
      if (!is.na(name) && nchar(name) >= 2 &&
          !name %in% c("unclassified", "uncultured", "metagenome", "NA", "Unknown") &&
          grepl("[a-zA-Z]", name)) {
        return(if (offset > 1) paste0("[", name, "]") else name)
      }
    }
  }
  paste0("Taxon_", substr(gsub("[^a-zA-Z0-9]", "", taxon), 1, 8))
}

# Drop covariates with no variation inside a subset
varying <- function(df, vars) vars[vapply(vars, function(v) length(unique(df[[v]])) > 1, logical(1))]

# Omega-squared and R2-per-df for every term of an adonis2 table
# omega2 = (SS_term - df_term * MS_res) / (SS_total + MS_res)
perm_effect_sizes <- function(tab) {
  tab  <- as.data.frame(tab)
  ms_r <- tab["Residual", "SumOfSqs"] / tab["Residual", "Df"]
  ss_t <- tab["Total", "SumOfSqs"]
  terms <- setdiff(rownames(tab), c("Residual", "Total"))
  data.frame(
    Term      = terms,
    Df        = tab[terms, "Df"],
    R2        = tab[terms, "R2"],
    R2_per_df = tab[terms, "R2"] / tab[terms, "Df"],
    Omega2    = pmax(0, (tab[terms, "SumOfSqs"] - tab[terms, "Df"] * ms_r) / (ss_t + ms_r)),
    F         = tab[terms, "F"],
    P         = tab[terms, "Pr(>F)"],
    row.names = NULL
  )
}

# Wilson 95% CI for a proportion
wilson_ci <- function(x, n, z = 1.96) {
  if (n == 0) return(c(NA, NA))
  p <- x / n; den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den
  hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(ctr - hw, ctr + hw)
}

write_sheets <- function(sheets, file) {
  wb <- createWorkbook()
  for (nm in names(sheets)) {
    sh <- substr(nm, 1, 31)
    addWorksheet(wb, sh)
    writeData(wb, sh, sheets[[nm]])
    setColWidths(wb, sh, cols = seq_len(max(1, ncol(sheets[[nm]]))), widths = "auto")
  }
  saveWorkbook(wb, file, overwrite = TRUE)
  msg("  Saved: %s", file)
}

msg("Config loaded. Output folder: %s", REV_DIR)

# ==============================================================================
# 1. BATCH
#    Sequencing batch effects (R2-20). Needs batch_map.tsv from Rev_01_batch_map.sh
# ==============================================================================
if (isTRUE(RUN$BATCH)) local({
  msg("\n################ SECTION 1: BATCH ################")

  suppressMessages(library(permute))

  # ---------------------------------------------------------------------------
  # 1. Parse batch map
  # ---------------------------------------------------------------------------
  bm <- read.delim(BATCH_MAP_FILE, header = FALSE, stringsAsFactors = FALSE,
                   col.names = c("SampleID", "Header", "UDP", "File"))
  hdr <- strsplit(sub("^@", "", bm$Header), ":")
  bm$Instrument <- vapply(hdr, `[`, "", 1)
  bm$RunNumber  <- vapply(hdr, `[`, "", 2)
  bm$Flowcell   <- vapply(hdr, `[`, "", 3)
  bm$Lane       <- vapply(hdr, `[`, "", 4)
  bm$Run        <- paste(bm$Instrument, bm$RunNumber, bm$Flowcell, sep = "_")
  bm$UDP_num    <- suppressWarnings(as.integer(sub("UDP", "", bm$UDP)))
  bm$IndexBlock <- ifelse(is.na(bm$UDP_num), NA, paste0("Block_", LETTERS[ceiling(bm$UDP_num / 96)]))

  meta <- load_meta()
  unmatched_fastq <- setdiff(bm$SampleID, meta$SampleID)
  unmatched_meta  <- setdiff(meta$SampleID, bm$SampleID)
  if (length(unmatched_fastq)) msg("FASTQ IDs not in metadata (%d): %s", length(unmatched_fastq), paste(head(unmatched_fastq, 20), collapse = ", "))
  if (length(unmatched_meta))  msg("Metadata IDs without FASTQ (%d): %s", length(unmatched_meta),  paste(head(unmatched_meta, 20),  collapse = ", "))

  meta <- merge(meta, bm[, c("SampleID", "Run", "Flowcell", "IndexBlock")], by = "SampleID")
  rownames(meta) <- meta$SampleID
  msg("Samples with batch info: %d", nrow(meta))

  # ---------------------------------------------------------------------------
  # 2. Batch x group balance
  # ---------------------------------------------------------------------------
  balance <- function(var) {
    tab <- table(meta$Status, meta[[var]])
    print(tab)
    if (ncol(tab) < 2) {
      msg("  -> %s has a single level: no batch confounding possible for this variable.", var)
      return(list(tab = as.data.frame.matrix(tab), test = data.frame(Variable = var, Levels = 1, P = NA, CramersV = NA)))
    }
    set.seed(SEED)
    ft <- fisher.test(tab, simulate.p.value = TRUE, B = 1e5)
    chi <- suppressWarnings(chisq.test(tab))
    V <- sqrt(as.numeric(chi$statistic) / (sum(tab) * (min(dim(tab)) - 1)))
    list(tab = as.data.frame.matrix(tab),
         test = data.frame(Variable = var, Levels = ncol(tab), Fisher_MC_P = ft$p.value, CramersV = V))
  }

  msg("\n=== Group x sequencing run ===");   bal_run   <- balance("Run")
  msg("\n=== Group x UDP index block ==="); bal_block <- balance("IndexBlock")

  # ---------------------------------------------------------------------------
  # 3. PERMANOVA: batch alone, disease + batch, disease with permutations within batch
  # ---------------------------------------------------------------------------
  dists <- lapply(load_dists(), subset_dist, ids = meta$SampleID)

  perm_batch <- function(var) {
    if (length(unique(meta[[var]])) < 2) return(NULL)
    out <- list()
    for (m in names(dists)) {
      d <- dists[[m]]
      set.seed(SEED)
      a_batch <- adonis2(as.formula(paste("d ~", var)), data = meta, permutations = N_PERM)
      a_both  <- adonis2(as.formula(paste("d ~ Status +", var)), data = meta, permutations = N_PERM, by = "margin")
      ctrl    <- how(nperm = N_PERM, blocks = factor(meta[[var]]))
      a_strat <- adonis2(d ~ Status, data = meta, permutations = ctrl)
      bd      <- permutest(betadisper(d, factor(meta[[var]])), permutations = N_PERM)
      out[[m]] <- data.frame(
        Metric = m, Batch_variable = var,
        Batch_R2 = a_batch$R2[1], Batch_P = a_batch$`Pr(>F)`[1],
        Disease_R2_adj_for_batch = a_both["Status", "R2"], Disease_P_adj_for_batch = a_both["Status", "Pr(>F)"],
        Batch_R2_adj_for_disease = a_both[var, "R2"],      Batch_P_adj_for_disease = a_both[var, "Pr(>F)"],
        Disease_P_perm_within_batch = a_strat$`Pr(>F)`[1],
        Batch_dispersion_P = bd$tab$`Pr(>F)`[1]
      )
      msg("  %-20s %-10s batch R2=%.3f (p=%.3f) | disease|batch R2=%.3f (p=%.3f) | within-batch perm p=%.3f",
          m, var, a_batch$R2[1], a_batch$`Pr(>F)`[1], a_both["Status", "R2"], a_both["Status", "Pr(>F)"], a_strat$`Pr(>F)`[1])
    }
    bind_rows(out)
  }

  msg("\n=== PERMANOVA with batch ===")
  perm_res <- bind_rows(perm_batch("Run"), perm_batch("IndexBlock"))

  # ---------------------------------------------------------------------------
  # 4. Save
  # ---------------------------------------------------------------------------
  write_sheets(list(
    Batch_map          = bm,
    Group_x_Run        = cbind(Group = rownames(bal_run$tab), bal_run$tab),
    Group_x_IndexBlock = cbind(Group = rownames(bal_block$tab), bal_block$tab),
    Balance_tests      = bind_rows(bal_run$test, bal_block$test),
    PERMANOVA_batch    = if (is.null(perm_res) || !nrow(perm_res)) data.frame(Note = "Single batch level - not testable") else perm_res
  ), file.path(REV_DIR, "Rev01_Batch_Effects.xlsx"))

  msg("\nINTERPRETATION GUIDE")
  msg("  * Single run           -> state 'all libraries sequenced on one MiSeq run' (strongest answer).")
  msg("  * Multiple runs, balanced (Fisher p > 0.05) and disease p unchanged with batch -> report as robust.")
  msg("  * Groups confounded with run (Cramer's V > 0.5) -> STOP and discuss before responding.")
})

# ==============================================================================
# 2. PERMANOVA_ROBUSTNESS
#    Effect sizes, Welch Wd*, balanced subsampling, Aitchison, non-viral & patients-only subsets (R2-2,3,5,9,10,11; R3-1) -> Supplementary Table S8
# ==============================================================================
if (isTRUE(RUN$PERMANOVA_ROBUSTNESS)) local({
  msg("\n################ SECTION 2: PERMANOVA_ROBUSTNESS ################")


  N_SUB      <- 500    # balanced-subsampling iterations (500 x 5 metrics ~ 10-20 min)
  N_PERM_SUB <- 199    # permutations inside each subsample

  meta  <- load_meta()
  dists <- load_dists()

  # C. Add robust Aitchison distance on unrarefied genus counts
  cts <- load_genus_counts()
  cts <- cts[rowSums(cts > 0) >= ceiling(0.10 * ncol(cts)), ]            # 10% prevalence
  dists[["Aitchison (robust, unrarefied)"]] <- vegdist(t(cts), method = "robust.aitchison")

  ids   <- Reduce(intersect, c(lapply(dists, labels), list(meta$SampleID)))
  meta  <- meta[ids, ]
  dists <- lapply(dists, subset_dist, ids = ids)
  msg("Samples: %d | metrics: %s", length(ids), paste(names(dists), collapse = ", "))

  PERMA_COVARS <- c("Gender", "Age_numeric", "BMI_numeric", "PPI_Use", "Metformin", "Rifaximin_Use",
                    "Lactulose_Use", "Statin_Use", "TG_numeric", "LDL_numeric", "HDL_numeric")

  # ---------------------------------------------------------------------------
  # A. Effect sizes
  # ---------------------------------------------------------------------------
  msg("\n=== A. Effect sizes (omega2, R2/df) ===")
  eff <- list()
  for (m in names(dists)) {
    d <- dists[[m]]
    set.seed(SEED)
    un  <- adonis2(d ~ Status, data = meta, permutations = N_PERM, by = "terms")   # by = "terms": row named Status (vegan >= 2.7 default differs)
    adj <- adonis2(as.formula(paste("d ~ Status +", paste(PERMA_COVARS, collapse = " + "))),
                   data = meta, permutations = N_PERM, by = "margin")
    eff[[m]] <- bind_rows(cbind(Metric = m, Model = "Unadjusted", perm_effect_sizes(un)),
                          cbind(Metric = m, Model = "Adjusted (marginal)", perm_effect_sizes(adj)))
    s <- eff[[m]] %>% dplyr::filter(Term == "Status")
    msg("  %-32s Status: R2 %.3f -> %.3f | omega2 %.3f -> %.3f", m, s$R2[1], s$R2[2], s$Omega2[1], s$Omega2[2])
  }
  eff_tab <- bind_rows(eff)
  # Rank of disease among all terms in the adjusted model, by omega2 and by R2/df
  rank_tab <- eff_tab %>% dplyr::filter(Model == "Adjusted (marginal)") %>% group_by(Metric) %>%
    dplyr::mutate(Rank_omega2 = rank(-Omega2), Rank_R2_per_df = rank(-R2_per_df)) %>%
    dplyr::filter(Term == "Status") %>% dplyr::select(Metric, R2, Omega2, R2_per_df, Rank_omega2, Rank_R2_per_df)
  print(as.data.frame(rank_tab))

  # ---------------------------------------------------------------------------
  # B1. Welch-type Wd* test
  # ---------------------------------------------------------------------------
  wd_star_stat <- function(D2, g) {
    g   <- droplevels(g); idx <- split(seq_along(g), g); k <- length(idx)
    n   <- lengths(idx)
    within <- vapply(idx, function(ii) sum(D2[ii, ii]) / 2, 0)          # sum of d^2 over pairs a<b
    s2  <- (within / n) / (n - 1)                                        # within-group variance
    w   <- n / s2; W <- sum(w)
    num <- 0
    for (i in 1:(k - 1)) for (j in (i + 1):k) {
      cd2 <- mean(D2[idx[[i]], idx[[j]]]) - within[i] / n[i]^2 - within[j] / n[j]^2   # ||c_i - c_j||^2
      num <- num + w[i] * w[j] * cd2
    }
    num <- num / W / (k - 1)
    den <- 1 + (2 * (k - 2) / (k^2 - 1)) * sum((1 - w / W)^2 / (n - 1))
    unname(num / den)
  }
  wd_star_test <- function(d, g, nperm = N_PERM) {
    D2 <- as.matrix(d)^2; g <- droplevels(factor(g))
    obs <- wd_star_stat(D2, g)
    set.seed(SEED)
    perm <- replicate(nperm, wd_star_stat(D2, sample(g)))
    c(Wd = obs, P = (1 + sum(perm >= obs)) / (nperm + 1))
  }

  msg("\n=== B1. Welch Wd* test (global + Healthy vs each) ===")
  welch <- list()
  for (m in names(dists)) {
    d <- dists[[m]]
    g <- wd_star_test(d, meta$Status)
    rows <- data.frame(Metric = m, Comparison = "Global (4 groups)", Wd = g[["Wd"]], P = g[["P"]])
    for (lv in GROUPS[-1]) {
      keep <- meta$Status %in% c("Healthy", lv)
      t2 <- wd_star_test(subset_dist(d, meta$SampleID[keep]), meta$Status[keep])
      rows <- rbind(rows, data.frame(Metric = m, Comparison = paste(lv, "vs Healthy"), Wd = t2[["Wd"]], P = t2[["P"]]))
    }
    rows$P_BH_pairwise <- NA
    pw <- rows$Comparison != "Global (4 groups)"
    rows$P_BH_pairwise[pw] <- p.adjust(rows$P[pw], "BH")
    welch[[m]] <- rows
    msg("  %-32s global Wd*=%.2f p=%.3f", m, g[["Wd"]], g[["P"]])
  }
  welch_tab <- bind_rows(welch); rownames(welch_tab) <- NULL

  # ---------------------------------------------------------------------------
  # B2. Balanced subsampling (17 per group)
  # ---------------------------------------------------------------------------
  msg("\n=== B2. Balanced subsampling: %d x (17 per group) ===", N_SUB)
  n_min <- min(table(meta$Status))
  set.seed(SEED)
  sub_ids <- replicate(N_SUB, unlist(lapply(split(meta$SampleID, meta$Status), sample, n_min)), simplify = FALSE)
  bal <- list()
  for (m in names(dists)) {
    dm <- as.matrix(dists[[m]])
    r <- t(vapply(sub_ids, function(s) {
      md <- meta[s, ]; d <- as.dist(dm[s, s])
      a  <- adonis2(d ~ Status, data = md, permutations = N_PERM_SUB)
      bd <- permutest(betadisper(d, md$Status), permutations = N_PERM_SUB)
      c(R2 = a$R2[1], P = a$`Pr(>F)`[1], Disp_P = bd$tab$`Pr(>F)`[1])
    }, numeric(3)))
    bal[[m]] <- data.frame(Metric = m, n_per_group = n_min, Iterations = N_SUB,
                           R2_median = median(r[, "R2"]), R2_IQR = paste(round(quantile(r[, "R2"], c(.25, .75)), 3), collapse = "-"),
                           P_median = median(r[, "P"]), Pct_P_lt_0.05 = 100 * mean(r[, "P"] < 0.05),
                           Dispersion_Pct_P_lt_0.05 = 100 * mean(r[, "Disp_P"] < 0.05))
    msg("  %-32s median R2=%.3f, median p=%.3f, %%p<0.05 = %.0f%%", m, bal[[m]]$R2_median, bal[[m]]$P_median, bal[[m]]$Pct_P_lt_0.05)
  }
  bal_tab <- bind_rows(bal)

  # ---------------------------------------------------------------------------
  # D. Sensitivity subsets
  # ---------------------------------------------------------------------------
  run_subset <- function(keep, label, covars) {
    md <- meta[keep, ]; md$Status <- droplevels(md$Status)
    if (any(table(md$Status) < MIN_N_SUBSET)) { msg("  [SKIP] %s (group n < %d)", label, MIN_N_SUBSET); return(NULL) }
    covars <- varying(md, covars)
    bind_rows(lapply(names(dists), function(m) {
      d <- subset_dist(dists[[m]], md$SampleID)
      set.seed(SEED)
      un  <- adonis2(d ~ Status, data = md, permutations = N_PERM)
      adj <- adonis2(as.formula(paste("d ~ Status +", paste(covars, collapse = " + "))), data = md,
                     permutations = N_PERM, by = "margin")
      bd  <- permutest(betadisper(d, md$Status), permutations = N_PERM)
      wd  <- wd_star_test(d, md$Status)
      data.frame(Subset = label, Groups = paste(names(table(md$Status)), table(md$Status), sep = "=", collapse = ", "),
                 Metric = m, R2_unadj = un$R2[1], P_unadj = un$`Pr(>F)`[1],
                 R2_adj = adj["Status", "R2"], P_adj = adj["Status", "Pr(>F)"],
                 PERMDISP_P = bd$tab$`Pr(>F)`[1], Welch_Wd_P = wd[["P"]], row.names = NULL)
    }))
  }

  msg("\n=== D. Sensitivity subsets ===")
  sens <- list()
  if (!all(is.na(meta$Viral))) sens$nonviral <- run_subset(!(meta$Viral %in% TRUE), "Non-viral (viral cirrhosis/HCC excluded)", PERMA_COVARS)
  sens$disease  <- run_subset(meta$Status != "Healthy", "Disease-only (MASLD, cirrhosis, HCC)", PERMA_COVARS)
  ctrl_old <- meta$Status == "Healthy" & meta$Age_numeric >= AGE_OVERLAP_MIN
  if (sum(ctrl_old) >= MIN_N_SUBSET) {
    maxage <- max(meta$Age_numeric[meta$Status == "Healthy"])
    sens$ageov <- run_subset(ctrl_old | (meta$Status == "MASLD" & meta$Age_numeric <= maxage),
                             sprintf("Age-overlap (controls >= %d y vs MASLD <= %d y)", AGE_OVERLAP_MIN, maxage), PERMA_COVARS)
  } else msg("  [SKIP] age-overlap: only %d controls aged >= %d", sum(ctrl_old), AGE_OVERLAP_MIN)
  sens_tab <- bind_rows(sens)
  if (nrow(sens_tab)) print(sens_tab[, c("Subset", "Metric", "R2_unadj", "P_unadj", "R2_adj", "P_adj", "PERMDISP_P", "Welch_Wd_P")])

  # ---------------------------------------------------------------------------
  # Save
  # ---------------------------------------------------------------------------
  write_sheets(list(A_Effect_sizes_all_terms = eff_tab, A_Disease_rank = rank_tab,
                    B1_Welch_Wd_star = welch_tab, B2_Balanced_subsampling = bal_tab,
                    D_Sensitivity_subsets = sens_tab,
                    Notes = data.frame(Note = c(
                      "Omega2 = (SS_term - df_term*MS_res)/(SS_total + MS_res); unbiased and comparable across terms with different df.",
                      "Wd*: Welch-type heteroscedastic PERMANOVA statistic (Hamidi et al. 2019, Microbiome 7:51); permutation p (999).",
                      sprintf("Balanced subsampling: %d draws of %d samples per group; %d permutations each.", N_SUB, n_min, N_PERM_SUB),
                      "Aitchison: vegan robust.aitchison on unrarefied genus counts, taxa present in >= 10% of samples."))),
               file.path(REV_DIR, "Rev03_PERMANOVA_Robustness.xlsx"))
})

# ==============================================================================
# 3. CONFOUNDERS
#    Age, BMI and sex: univariate, adjusted, patients-only, controls-only (R2-3,4,5) -> Supplementary Table S9
# ==============================================================================
if (isTRUE(RUN$CONFOUNDERS)) local({
  msg("\n################ SECTION 3: CONFOUNDERS ################")
  # Reviewer 2 comments 3-5 (age, BMI, sex) -> Supplementary Table S9
  meta  <- load_meta()
  dists <- load_dists()
  cts <- load_genus_counts(); cts <- cts[rowSums(cts > 0) >= ceiling(0.10 * ncol(cts)), ]
  dists[["Aitchison (robust, unrarefied)"]] <- vegdist(t(cts), "robust.aitchison")
  cv <- c("Gender", "Age_numeric", "BMI_numeric", "PPI_Use", "Metformin", "Rifaximin_Use",
          "Lactulose_Use", "Statin_Use", "TG_numeric", "LDL_numeric", "HDL_numeric")
  msg("Age range by group:"); print(sapply(split(meta$Age_numeric, meta$Status), range))
  msg("Patients aged <= %d (max control age): %d", max(meta$Age_numeric[meta$Status == "Healthy"]),
      sum(meta$Status != "Healthy" & meta$Age_numeric <= max(meta$Age_numeric[meta$Status == "Healthy"])))
  rows <- list()
  for (m in names(dists)) {
    ids <- intersect(labels(dists[[m]]), meta$SampleID); md <- meta[ids, ]; dd <- subset_dist(dists[[m]], ids)
    set.seed(SEED)
    uni  <- lapply(c(Age = "Age_numeric", BMI = "BMI_numeric", Sex = "Gender"), function(v)
            adonis2(as.formula(paste("dd ~", v)), data = md, permutations = N_PERM))
    full <- adonis2(as.formula(paste("dd ~ Status +", paste(cv, collapse = " + "))), data = md,
                    permutations = N_PERM, by = "margin")
    es   <- perm_effect_sizes(full)
    pt   <- md$Status != "Healthy"; mp <- md[pt, ]; mp$Status <- droplevels(mp$Status); dp <- subset_dist(dd, mp$SampleID)
    pat  <- adonis2(as.formula(paste("dp ~ Status +", paste(varying(mp, cv), collapse = " + "))), data = mp,
                    permutations = N_PERM, by = "margin")
    h    <- md$Status == "Healthy"; mh <- md[h, ]; dh <- subset_dist(dd, mh$SampleID)
    ctl  <- adonis2(dh ~ Age_numeric + BMI_numeric + Gender, data = mh, permutations = N_PERM, by = "margin")
    for (v in c("Status", "Age_numeric", "BMI_numeric", "Gender")) {
      lab <- c(Status = "Disease status", Age_numeric = "Age", BMI_numeric = "BMI", Gender = "Sex")[[v]]
      u <- if (v == "Status") NULL else uni[[c(Age_numeric = "Age", BMI_numeric = "BMI", Gender = "Sex")[[v]]]]
      rows[[length(rows) + 1]] <- data.frame(
        Metric = m, Variable = lab,
        Univariate_R2 = if (is.null(u)) NA else u$R2[1], Univariate_P = if (is.null(u)) NA else u$`Pr(>F)`[1],
        Adjusted_R2 = full[v, "R2"], Adjusted_Omega2 = es$Omega2[es$Term == v], Adjusted_P = full[v, "Pr(>F)"],
        PatientsOnly_R2 = pat[v, "R2"], PatientsOnly_P = pat[v, "Pr(>F)"],
        ControlsOnly_P = if (v == "Status") NA else ctl[v, "Pr(>F)"], row.names = NULL)
    }
    msg("  %-32s done", m)
  }
  conf <- bind_rows(rows)
  print(conf %>% dplyr::mutate(across(where(is.numeric), ~ signif(.x, 3))))
  write_sheets(list(S9_Age_BMI_Sex = conf,
                    Notes = data.frame(Note = c(
                      "Adjusted = marginal PERMANOVA with disease status + 11 covariates (n = all samples).",
                      "PatientsOnly = same model in MASLD, cirrhosis and HCC only (ages overlap).",
                      "ControlsOnly = age + BMI + sex within healthy controls.",
                      "Omega2 is comparable across terms with different degrees of freedom."))),
               file.path(REV_DIR, "Supp_S9_Age_BMI_Sex_PERMANOVA.xlsx"))
})

# ==============================================================================
# 4. DA_ADJUSTED
#    ANCOM-BC2 + DESeq2: unadjusted, adjusted, non-viral, patients-only models (R2-2,3,6,7,8,24,25) -> Supplementary Table S10
# ==============================================================================
if (isTRUE(RUN$DA_ADJUSTED)) local({
  msg("\n################ SECTION 4: DA_ADJUSTED ################")

  suppressMessages({
    library(ANCOMBC)                    # >= 2.12
    library(DESeq2)
  })

  # ANCOM-BC2 settings (reported verbatim in Methods)
  ANC <- list(prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05, struc_zero = TRUE,
              neg_lb = FALSE, alpha = 0.05, p_adj_method = "BH", pseudo_sens = TRUE)   # neg_lb FALSE: HCC n = 17 (package advises TRUE only for large groups)
  # DESeq2 settings (unchanged from Supp Fig S2 script, but ONE model with all groups)
  DESEQ_MIN_COUNT <- 3; DESEQ_MIN_PCT <- 0.05

  KEY_TAXA <- c("Bacteroides", "Clostridium_sensu_stricto_1", "Lactobacillus", "Streptococcus",
                "Negativibacillus", "[Rhodospirillales]", "[Bacteroidales]")

  # ---------------------------------------------------------------------------
  # 1. Data
  # ---------------------------------------------------------------------------
  meta <- load_meta()
  cts  <- load_genus_counts()
  common <- intersect(colnames(cts), meta$SampleID)
  cts  <- cts[, common]; meta <- meta[common, ]
  msg("Samples: %d | genus-level taxa: %d", ncol(cts), nrow(cts))
  print(table(meta$Status))

  # ---------------------------------------------------------------------------
  # 2. Model runners
  # ---------------------------------------------------------------------------
  run_ancombc2 <- function(m, cts_sub, covars, ref, label) {
    m$Status <- relevel(droplevels(m$Status), ref = ref)
    covars   <- varying(m, covars)
    fml      <- paste(c("Status", covars), collapse = " + ")
    msg("\n[ANCOM-BC2] %s | n = %d | formula: %s", label, nrow(m), fml)
    set.seed(SEED)
    # Matrix input (taxa x samples) + meta_data: avoids the mia/TSE conversion layer
    out <- ancombc2(data = cts_sub[, m$SampleID], taxa_are_rows = TRUE, meta_data = m,
                    fix_formula = fml, rand_formula = NULL,
                    p_adj_method = ANC$p_adj_method, pseudo_sens = ANC$pseudo_sens,
                    prv_cut = ANC$prv_cut, lib_cut = ANC$lib_cut, s0_perc = ANC$s0_perc,
                    group = "Status", struc_zero = ANC$struc_zero, neg_lb = ANC$neg_lb,
                    alpha = ANC$alpha, n_cl = 1, verbose = FALSE,
                    global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE)
    res <- out$res
    lv  <- setdiff(levels(m$Status), ref)
    long <- bind_rows(lapply(lv, function(l) {
      g <- function(prefix) { col <- paste0(prefix, "Status", l); if (col %in% names(res)) res[[col]] else NA }
      data.frame(Model = label, Contrast = paste(l, "vs", ref), Taxon = res$taxon,
                 Genus = vapply(res$taxon, clean_genus, ""),
                 ANC_lfc = g("lfc_"), ANC_se = g("se_"), ANC_p = g("p_"), ANC_q = g("q_"),
                 ANC_diff = g("diff_"), ANC_passed_ss = g("passed_ss_"),
                 stringsAsFactors = FALSE)
    }))
    long$ANC_sig <- long$ANC_q < ANC$alpha & (is.na(long$ANC_passed_ss) | long$ANC_passed_ss)
    list(long = long, formula = fml, n_taxa_tested = nrow(res), zero_ind = out$zero_ind)
  }

  run_deseq2 <- function(m, cts_sub, covars, ref, label) {
    m$Status <- relevel(droplevels(m$Status), ref = ref)
    covars   <- varying(m, covars)
    x <- round(cts_sub[, m$SampleID])
    keep <- rowSums(x >= DESEQ_MIN_COUNT) >= ceiling(DESEQ_MIN_PCT * ncol(x))
    x <- x[keep, ]
    design <- as.formula(paste("~", paste(c(covars, "Status"), collapse = " + ")))
    msg("[DESeq2]    %s | taxa after filter: %d | design: %s", label, nrow(x), deparse(design))
    dds <- DESeqDataSetFromMatrix(x, as.data.frame(m), design)
    dds <- estimateSizeFactors(dds, type = "poscounts")   # = geometric-mean-of-positives (as before)
    dds <- DESeq(dds, fitType = "local", quiet = TRUE)
    bind_rows(lapply(setdiff(levels(m$Status), ref), function(l) {
      r <- results(dds, contrast = c("Status", l, ref), cooksCutoff = FALSE, independentFiltering = FALSE)
      data.frame(Model = label, Contrast = paste(l, "vs", ref), Taxon = rownames(r),
                 DES_lfc = r$log2FoldChange, DES_se = r$lfcSE, DES_p = r$pvalue, DES_q = r$padj,
                 stringsAsFactors = FALSE)
    }))
  }

  run_model <- function(m, covars, ref, label) {
    if (any(table(droplevels(m$Status)) < MIN_N_SUBSET)) {
      msg("\n[SKIP] %s: a group has n < %d", label, MIN_N_SUBSET); print(table(droplevels(m$Status)))
      return(NULL)
    }
    a <- run_ancombc2(m, cts, covars, ref, label)
    d <- run_deseq2(m, cts, covars, ref, label)
    tab <- left_join(a$long, d, by = c("Model", "Contrast", "Taxon"))
    tab$DES_sig    <- !is.na(tab$DES_q) & tab$DES_q < 0.05
    tab$Concordant <- tab$ANC_sig & tab$DES_sig & sign(tab$ANC_lfc) == sign(tab$DES_lfc)
    list(table = tab, formula = a$formula, n_taxa = a$n_taxa_tested, n = nrow(m),
         groups = paste(names(table(droplevels(m$Status))), table(droplevels(m$Status)), sep = "=", collapse = ", "))
  }

  # ---------------------------------------------------------------------------
  # 3. Define subsets
  # ---------------------------------------------------------------------------
  models <- list()
  models$M1_Full_unadjusted <- run_model(meta, character(0), "Healthy", "M1_Full_unadjusted")
  models$M2_Full_adj_core   <- run_model(meta, COVARS_CORE,  "Healthy", "M2_Full_adj_core")
  models$M3_Full_adj_full   <- run_model(meta, COVARS_FULL,  "Healthy", "M3_Full_adj_full")
  # Age is near-collinear with Healthy (controls 22-39 y; 6 patients <= 39 y): SEs double when it enters.
  models$M2b_Full_adj_noAge <- run_model(meta, setdiff(COVARS_CORE, "Age_z"), "Healthy", "M2b_Full_adj_noAge")

  if (all(is.na(meta$Viral))) {
    msg("\n[SKIP] M4: etiology column '%s' not found - set ETIOLOGY_COL in Rev_00.", ETIOLOGY_COL)
  } else {
    nv <- meta[!(meta$Viral %in% TRUE), ]
    msg("\nNon-viral subset removes %d viral cases", sum(meta$Viral, na.rm = TRUE))
    models$M4u_NonViral_unadjusted <- run_model(nv, character(0), "Healthy", "M4u_NonViral_unadjusted")   # direct comparison with M1
    models$M4_NonViral_adj_core    <- run_model(nv, COVARS_CORE,  "Healthy", "M4_NonViral_adj_core")
  }

  dis <- meta[meta$Status != "Healthy", ]
  models$M5_DiseaseOnly_adj_core <- run_model(dis, COVARS_CORE, "MASLD", "M5_DiseaseOnly_adj_core")

  ctrl_old <- meta$Status == "Healthy" & meta$Age_numeric >= AGE_OVERLAP_MIN
  msg("\nControls aged >= %d: %d (control age range %s)", AGE_OVERLAP_MIN, sum(ctrl_old),
      paste(range(meta$Age_numeric[meta$Status == "Healthy"]), collapse = "-"))
  if (sum(ctrl_old) >= MIN_N_SUBSET) {
    max_ctrl_age <- max(meta$Age_numeric[meta$Status == "Healthy"])
    ao <- meta[ctrl_old | (meta$Status == "MASLD" & meta$Age_numeric <= max_ctrl_age), ]
    models$M6_AgeOverlap_adj_core <- run_model(ao, COVARS_CORE, "Healthy", "M6_AgeOverlap_adj_core")
  } else {
    msg("[SKIP] M6 age-overlap: too few older controls. Report this explicitly as a limitation.")
  }
  models <- models[!vapply(models, is.null, logical(1))]

  # ---------------------------------------------------------------------------
  # 4. Summaries
  # ---------------------------------------------------------------------------
  all_tab <- bind_rows(lapply(models, `[[`, "table"))

  summary_tab <- all_tab %>%
    group_by(Model, Contrast) %>%
    dplyr::summarise(ANCOMBC2_sig = sum(ANC_sig, na.rm = TRUE),
              DESeq2_sig   = sum(DES_sig, na.rm = TRUE),
              Concordant   = sum(Concordant, na.rm = TRUE),
              Concordance_pct_of_ANCOMBC2 = round(100 * Concordant / pmax(1, ANCOMBC2_sig), 1),
              .groups = "drop")
  print(as.data.frame(summary_tab))

  stars <- function(q) ifelse(is.na(q), "", ifelse(q < 0.001, "***", ifelse(q < 0.01, "**", ifelse(q < 0.05, "*", ""))))
  key_tab <- all_tab %>%
    dplyr::filter(Genus %in% KEY_TAXA) %>%
    dplyr::mutate(ANCOMBC2 = sprintf("%.2f%s", ANC_lfc, stars(ifelse(ANC_sig, ANC_q, NA))),
           DESeq2   = sprintf("%.2f%s", DES_lfc, stars(DES_q))) %>%
    dplyr::select(Genus, Model, Contrast, ANCOMBC2, DESeq2, Concordant) %>%
    dplyr::arrange(Genus, Contrast, Model)

  # Robustness: of taxa significant in M1, how many keep q<0.05 and same sign in M2/M3?
  m1 <- all_tab %>% dplyr::filter(Model == "M1_Full_unadjusted", ANC_sig) %>% dplyr::select(Contrast, Taxon, Genus, lfc_M1 = ANC_lfc)
  robust_tab <- m1
  for (mm in intersect(c("M2_Full_adj_core", "M3_Full_adj_full", "M4_NonViral_adj_core"), names(models))) {
    x <- all_tab %>% dplyr::filter(Model == mm) %>% transmute(Contrast, Taxon, sig = ANC_sig, lfc = ANC_lfc)
    names(x)[3:4] <- paste0(c("sig_", "lfc_"), mm)
    robust_tab <- left_join(robust_tab, x, by = c("Contrast", "Taxon"))
  }

  params <- data.frame(
    Parameter = c("Package", "fix_formula (per model)", "group", "p_adj_method", "prv_cut", "lib_cut",
                  "s0_perc", "struc_zero", "neg_lb", "pseudo_sens", "alpha",
                  "FDR scope", "DESeq2 size factors", "DESeq2 filter", "DESeq2 fitType", "DESeq2 cooksCutoff / independentFiltering",
                  "Continuous covariates", "Seed"),
    Value = c(paste0("ANCOMBC ", packageVersion("ANCOMBC"), "; DESeq2 ", packageVersion("DESeq2")),
              paste(names(models), vapply(models, `[[`, "", "formula"), sep = ": ", collapse = " | "),
              "Status", ANC$p_adj_method, ANC$prv_cut, ANC$lib_cut, ANC$s0_perc, ANC$struc_zero, ANC$neg_lb,
              ANC$pseudo_sens, ANC$alpha,
              "BH across taxa, separately within each contrast (one model containing all groups)",
              "poscounts", sprintf(">= %d reads in >= %.0f%% of samples", DESEQ_MIN_COUNT, 100 * DESEQ_MIN_PCT),
              "local", "FALSE / FALSE", "z-scored (log2FC per 1 SD)", SEED))

  model_info <- data.frame(Model = names(models),
                           N = sapply(models, `[[`, "n"),
                           Groups = vapply(models, `[[`, "", "groups"),
                           Taxa_tested_ANCOMBC2 = sapply(models, `[[`, "n_taxa"),
                           Formula = vapply(models, `[[`, "", "formula"))

  write_sheets(list(Summary = summary_tab, Model_info = model_info, Key_7_taxa = key_tab,
                    Robustness_of_M1_hits = robust_tab, All_results = all_tab, Parameters = params),
               file.path(REV_DIR, "Rev02_DA_Adjusted_Sensitivity.xlsx"))
  save(models, all_tab, file = file.path(REV_DIR, "Rev02_DA_models.RData"))

  msg("\nNOTE for M2/M3: PPI/metformin/statin are 0%% in controls and age barely overlaps,")
  msg("so Healthy-vs-disease adjusted LFCs partly rely on extrapolation. M5 (disease-only)")
  msg("and M6 (age-overlap) are the cleaner tests of whether signals reflect disease.")
})

# ==============================================================================
# 5. DYSBIOSIS_INDEX
#    Non-circular Dysbiosis Index: LOO centroid, split-half, partial correlations (R2-12,13)
# ==============================================================================
if (isTRUE(RUN$DYSBIOSIS_INDEX)) local({
  msg("\n################ SECTION 5: DYSBIOSIS_INDEX ################")

  suppressMessages(library(ppcor))

  DI_PSEUDO  <- 1      # <<< CHECK: pseudocount used in original DI script
  DI_PRV     <- 0.10   # <<< CHECK: prevalence filter used in original DI script (0 = none)
  DI_VAR     <- 0.80
  N_SPLIT    <- 200
  ORIG_SCORES <- file.path(QIIME_DIR, "dysbiosis_analysis/Dysbiosis_Index_Scores.csv")
  CLIN_VARS  <- c("BMI", "Age", "GGT", "ALT", "AST", "Hgb_A1C", "TG", "BUN", "HDL", "Platelet")

  meta <- load_meta()
  cts  <- load_genus_counts()
  common <- intersect(colnames(cts), meta$SampleID)
  cts <- cts[, common]; meta <- meta[common, ]
  if (DI_PRV > 0) cts <- cts[rowSums(cts > 0) >= ceiling(DI_PRV * ncol(cts)), ]

  clr <- function(m) { l <- log(m); sweep(l, 1, rowMeans(l)) }
  X   <- clr(t(cts) + DI_PSEUDO)                       # samples x taxa
  pca <- prcomp(X, center = TRUE, scale. = FALSE)
  k   <- which(cumsum(pca$sdev^2) / sum(pca$sdev^2) >= DI_VAR)[1]
  S   <- pca$x[, 1:k, drop = FALSE]
  msg("PCs retained (>= %.0f%% variance): %d", 100 * DI_VAR, k)

  is_h <- meta$Status == "Healthy"
  dist_to <- function(S, ctr) sqrt(rowSums(sweep(S, 2, ctr)^2))
  rescale <- function(x) 100 * (x - min(x)) / (max(x) - min(x))

  # ---------------------------------------------------------------------------
  # 0. Original (in-sample) DI and reproduction check
  # ---------------------------------------------------------------------------
  di_orig_raw <- dist_to(S, colMeans(S[is_h, , drop = FALSE]))
  if (file.exists(ORIG_SCORES)) {
    o <- read.csv(ORIG_SCORES, stringsAsFactors = FALSE, check.names = FALSE)
    idc <- intersect(c("SampleID", "id", "sample_id", "Sample"), names(o))[1]
    if (is.na(idc)) idc <- names(o)[1]
    mm <- match(names(di_orig_raw), o[[idc]])
    if (sum(!is.na(mm)) > 10) {
      r <- cor(di_orig_raw[!is.na(mm)], o$Dysbiosis_Index[mm[!is.na(mm)]], method = "spearman")
      msg("Reproduction check vs published DI: Spearman rho = %.3f %s", r,
          ifelse(r > 0.95, "(OK)", "(LOW - align DI_PSEUDO / DI_PRV with original script)"))
    }
  }

  # ---------------------------------------------------------------------------
  # 1. LOO centroid for controls
  # ---------------------------------------------------------------------------
  di_loo_raw <- di_orig_raw
  hS <- S[is_h, , drop = FALSE]; nh <- nrow(hS); sum_h <- colSums(hS)
  for (i in seq_len(nh)) {
    ctr_i <- (sum_h - hS[i, ]) / (nh - 1)
    di_loo_raw[which(is_h)[i]] <- sqrt(sum((hS[i, ] - ctr_i)^2))
  }
  meta$DI_original <- rescale(di_orig_raw)
  meta$DI_LOO      <- rescale(di_loo_raw)

  grp_tests <- function(v) {
    kw <- kruskal.test(meta[[v]] ~ meta$Status)
    pw <- pairwise.wilcox.test(meta[[v]], meta$Status, p.adjust.method = "BH", exact = FALSE)$p.value
    med <- tapply(meta[[v]], meta$Status, median)
    data.frame(DI = v, Healthy_median = med["Healthy"], MASLD_median = med["MASLD"],
               Cirrhosis_median = med["Cirrhosis"], HCC_median = med["HCC"],
               KW_P = kw$p.value,
               q_MASLD_vs_H = pw["MASLD", "Healthy"], q_Cirr_vs_H = pw["Cirrhosis", "Healthy"],
               q_HCC_vs_H = pw["HCC", "Healthy"], row.names = NULL)
  }
  loo_tab <- bind_rows(grp_tests("DI_original"), grp_tests("DI_LOO"))
  msg("\n=== Original vs LOO DI ==="); print(loo_tab)

  # ---------------------------------------------------------------------------
  # 2. Split-half validation
  # ---------------------------------------------------------------------------
  set.seed(SEED)
  h_ids <- which(is_h)
  split_res <- t(replicate(N_SPLIT, {
    ref  <- sample(h_ids, floor(length(h_ids) / 2)); held <- setdiff(h_ids, ref)
    di   <- dist_to(S, colMeans(S[ref, , drop = FALSE]))
    sapply(c("MASLD", "Cirrhosis", "HCC"), function(g) {
      wilcox.test(di[meta$Status == g], di[held], exact = FALSE)$p.value
    })
  }))
  split_tab <- data.frame(Comparison = paste(colnames(split_res), "vs held-out Healthy"),
                          Median_P = apply(split_res, 2, median),
                          Pct_P_lt_0.05 = 100 * colMeans(split_res < 0.05),
                          Iterations = N_SPLIT, row.names = NULL)
  msg("\n=== Split-half validation ==="); print(split_tab)

  # ---------------------------------------------------------------------------
  # 3. Disease-only comparisons (no healthy centroid involvement in the contrast)
  # ---------------------------------------------------------------------------
  dm <- meta[!is_h, ]; dm$Status <- droplevels(dm$Status)
  kw_d <- kruskal.test(DI_LOO ~ Status, data = dm)
  pw_d <- pairwise.wilcox.test(dm$DI_LOO, dm$Status, p.adjust.method = "BH", exact = FALSE)$p.value
  dis_tab <- data.frame(Test = c("Kruskal-Wallis (MASLD, cirrhosis, HCC)", "Cirrhosis vs MASLD", "HCC vs MASLD", "HCC vs Cirrhosis"),
                        P = c(kw_d$p.value, pw_d["Cirrhosis", "MASLD"], pw_d["HCC", "MASLD"], pw_d["HCC", "Cirrhosis"]))
  msg("\n=== Disease-only DI comparisons ==="); print(dis_tab)

  # ---------------------------------------------------------------------------
  # 4. Clinical correlations: pooled, partial (adjusted for group), within-group
  # ---------------------------------------------------------------------------
  Z <- model.matrix(~ Status, data = meta)[, -1, drop = FALSE]
  cor_rows <- list()
  for (v in CLIN_VARS) {
    if (!v %in% names(meta)) { msg("  [skip] %s not in metadata", v); next }
    y <- suppressWarnings(as.numeric(meta[[v]])); ok <- !is.na(y)
    pooled <- cor.test(meta$DI_LOO[ok], y[ok], method = "spearman", exact = FALSE)
    part   <- pcor.test(meta$DI_LOO[ok], y[ok], Z[ok, , drop = FALSE], method = "spearman")
    within <- sapply(GROUPS, function(g) {
      s <- ok & meta$Status == g
      if (sum(s) < 10 || length(unique(y[s])) < 3) return(NA)
      suppressWarnings(cor(meta$DI_LOO[s], y[s], method = "spearman"))
    })
    cor_rows[[v]] <- data.frame(Variable = v, N = sum(ok),
                                Pooled_rho = unname(pooled$estimate), Pooled_P = pooled$p.value,
                                Partial_rho = part$estimate, Partial_P = part$p.value,
                                t(setNames(within, paste0("Within_", GROUPS, "_rho"))), check.names = FALSE)
  }
  cor_tab <- bind_rows(cor_rows)
  cor_tab$Pooled_q  <- p.adjust(cor_tab$Pooled_P, "BH")
  cor_tab$Partial_q <- p.adjust(cor_tab$Partial_P, "BH")
  msg("\n=== DI-clinical correlations (pooled vs adjusted for group) ===")
  print(cor_tab[, c("Variable", "Pooled_rho", "Pooled_q", "Partial_rho", "Partial_q")])

  write_sheets(list(LOO_vs_original = loo_tab, Split_half = split_tab, Disease_only = dis_tab,
                    Clinical_correlations = cor_tab,
                    Scores = meta[, c("SampleID", "Status", "DI_original", "DI_LOO")],
                    Settings = data.frame(Setting = c("Pseudocount", "Prevalence filter", "Variance retained", "PCs", "Split-half iterations"),
                                          Value = c(DI_PSEUDO, DI_PRV, DI_VAR, k, N_SPLIT))),
               file.path(REV_DIR, "Rev04_Dysbiosis_Index_NonCircular.xlsx"))
})

# ==============================================================================
# 6. RAREFACTION
#    DADA2 retention, rarefaction-depth sensitivity, Chao1 = Observed check (R2-26,27). Needs QIIME2 exports (see section header)
# ==============================================================================
if (isTRUE(RUN$RAREFACTION)) local({
  msg("\n################ SECTION 6: RAREFACTION ################")

  suppressMessages({ library(phyloseq); library(ape) })
  has_picante <- requireNamespace("picante", quietly = TRUE)

  DEPTHS <- c(13000, 5000, 8000, 10000, 20000)   # first = depth used in the paper; 20,000 drops some samples
  N_ITER <- 10

  meta <- load_meta()

  # ---------------------------------------------------------------------------
  # a. DADA2 retention
  # ---------------------------------------------------------------------------
  if (file.exists(DADA2_STATS)) {
    st <- read.delim(DADA2_STATS, check.names = FALSE, stringsAsFactors = FALSE)
    st <- st[!grepl("^#q2", st[[1]]), ]
    names(st)[1] <- "SampleID"
    num <- setdiff(names(st), "SampleID"); st[num] <- lapply(st[num], as.numeric)
    st$Retained_pct <- 100 * st$`non-chimeric` / st$input
    st <- merge(st, meta[, c("SampleID", "Status")], by = "SampleID")
    dada_sum <- st %>% group_by(Status) %>%
      dplyr::summarise(n = n(), Input_median = median(input), Filtered_median = median(filtered),
                Merged_median = median(merged), NonChimeric_median = median(`non-chimeric`),
                NonChimeric_min = min(`non-chimeric`), Retained_pct_median = median(Retained_pct),
                Below_5000 = sum(`non-chimeric` < 5000), Below_10000 = sum(`non-chimeric` < 10000), .groups = "drop")
    msg("\n=== DADA2 retention ==="); print(as.data.frame(dada_sum))
    msg("Lowest 10 samples after DADA2:"); print(head(st[order(st$`non-chimeric`), c("SampleID", "Status", "input", "filtered", "merged", "non-chimeric", "Retained_pct")], 10))
    msg("Kruskal-Wallis: retained reads ~ group, p = %.3g", kruskal.test(`non-chimeric` ~ Status, data = st)$p.value)
  } else { st <- data.frame(); dada_sum <- data.frame(Note = "DADA2 stats not exported yet") }

  # ---------------------------------------------------------------------------
  # b. Rarefaction sensitivity
  # ---------------------------------------------------------------------------
  asv <- read.table(ASV_TABLE_FILE, header = TRUE, row.names = 1, sep = "\t", skip = 1,
                    check.names = FALSE, comment.char = "")
  asv <- as.matrix(asv); storage.mode(asv) <- "numeric"
  tree <- read.tree(TREE_FILE)
  ids  <- intersect(colnames(asv), meta$SampleID)
  asv  <- asv[, ids]; md <- meta[ids, ]
  depth <- colSums(asv)
  msg("\nASV table: %d ASVs x %d samples | depth min %d, median %d", nrow(asv), ncol(asv), min(depth), round(median(depth)))

  # c. Chao1 vs observed: count singletons per sample
  singletons <- colSums(asv == 1)
  msg("Samples containing any singleton ASV: %d of %d (Chao1 = Observed when this is 0)", sum(singletons > 0), ncol(asv))

  retained <- t(sapply(DEPTHS, function(dp) table(factor(md$Status[depth >= dp], levels = GROUPS))))
  rownames(retained) <- DEPTHS
  msg("\nSamples retained per depth:"); print(retained)

  alpha_rows <- list(); beta_rows <- list()
  for (dp in DEPTHS) {
    keep <- names(depth)[depth >= dp]; mdd <- md[keep, ]
    if (any(table(mdd$Status) < 10)) { msg("[skip depth %d] a group falls below 10 samples", dp); next }
    for (it in seq_len(N_ITER)) {
      set.seed(SEED + it)
      rar <- rrarefy(t(asv[, keep]), sample = dp)               # samples x ASVs
      rar <- rar[, colSums(rar) > 0, drop = FALSE]
      # alpha
      al <- data.frame(Observed = rowSums(rar > 0),
                       Shannon  = diversity(rar, "shannon"),
                       Simpson  = diversity(rar, "simpson"))
      if (has_picante) {
        tr <- keep.tip(tree, intersect(tree$tip.label, colnames(rar)))
        al$FaithPD <- picante::pd(rar[, tr$tip.label], tr, include.root = TRUE)$PD
      }
      kw <- sapply(al, function(v) kruskal.test(v ~ mdd$Status)$p.value)
      alpha_rows[[length(alpha_rows) + 1]] <- data.frame(Depth = dp, Iter = it, Metric = names(kw), KW_P = kw,
                                                         BH_q = p.adjust(kw, "BH"))
      # beta
      ps <- phyloseq(otu_table(t(rar), taxa_are_rows = TRUE),
                     phy_tree(keep.tip(tree, intersect(tree$tip.label, colnames(rar)))))
      dl <- list("Bray-Curtis" = vegdist(rar, "bray"),
                 "Jaccard"     = vegdist(rar, "jaccard", binary = TRUE),
                 "Unweighted UniFrac" = UniFrac(ps, weighted = FALSE),
                 "Weighted UniFrac"   = UniFrac(ps, weighted = TRUE, normalized = FALSE))   # unnormalised, as QIIME 2
      for (m in names(dl)) {
        d  <- as.dist(as.matrix(dl[[m]])[keep, keep])
        a  <- adonis2(d ~ Status, data = mdd, permutations = 499)
        bd <- permutest(betadisper(d, mdd$Status), permutations = 499)
        beta_rows[[length(beta_rows) + 1]] <- data.frame(Depth = dp, Iter = it, Metric = m,
                                                         R2 = a$R2[1], P = a$`Pr(>F)`[1], PERMDISP_P = bd$tab$`Pr(>F)`[1])
      }
    }
    msg("  depth %d done (n = %d)", dp, length(keep))
  }
  alpha_all <- bind_rows(alpha_rows); beta_all <- bind_rows(beta_rows)

  alpha_sum <- alpha_all %>% group_by(Depth, Metric) %>%
    dplyr::summarise(KW_P_median = median(KW_P), BH_q_median = median(BH_q), Pct_q_lt_0.05 = 100 * mean(BH_q < 0.05), .groups = "drop")
  beta_sum <- beta_all %>% group_by(Depth, Metric) %>%
    dplyr::summarise(R2_median = median(R2), P_median = median(P), Pct_P_lt_0.05 = 100 * mean(P < 0.05),
              PERMDISP_P_median = median(PERMDISP_P), .groups = "drop")
  msg("\n=== Alpha (median over %d rarefactions) ===", N_ITER); print(as.data.frame(alpha_sum))
  msg("\n=== Beta (median over %d rarefactions) ===",  N_ITER); print(as.data.frame(beta_sum))

  write_sheets(list(DADA2_by_group = dada_sum, DADA2_per_sample = st,
                    Retained_per_depth = cbind(Depth = DEPTHS, as.data.frame(retained)),
                    Alpha_summary = alpha_sum, Beta_summary = beta_sum,
                    Alpha_all_iterations = alpha_all, Beta_all_iterations = beta_all,
                    Singletons = data.frame(SampleID = names(singletons), Singleton_ASVs = singletons)),
               file.path(REV_DIR, "Rev05_Rarefaction_Sensitivity.xlsx"))

})

# ==============================================================================
# 7. CLINICAL
#    Control profile, MELD/Child-Pugh, MASLD criteria, HCC chart-review template (R2-5,18,19,22,30)
# ==============================================================================
if (isTRUE(RUN$CLINICAL)) local({
  msg("\n################ SECTION 7: CLINICAL ################")

  meta <- load_meta()

  yes <- function(x) grepl("^(1|y|yes|true|positive|pos)$", tolower(trimws(as.character(x))))
  num <- function(x) suppressWarnings(as.numeric(x))
  med_iqr <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return("--")
    sprintf("%.1f (%.1f-%.1f)", median(x), quantile(x, .25), quantile(x, .75)) }
  srow  <- function(label, x) data.frame(Variable = label, as.list(setNames(as.character(x), names(x))), check.names = FALSE)
  n_pct <- function(l) { l <- l[!is.na(l)]; sprintf("%d (%.1f)", sum(l), 100 * mean(l)) }

  # ---------------------------------------------------------------------------
  # 1. Control age profile
  # ---------------------------------------------------------------------------
  age_tab <- meta %>% group_by(Status) %>%
    dplyr::summarise(n = n(), Min = min(Age_numeric), Q1 = quantile(Age_numeric, .25), Median = median(Age_numeric),
              Q3 = quantile(Age_numeric, .75), Max = max(Age_numeric),
              n_ge_40 = sum(Age_numeric >= 40), n_ge_50 = sum(Age_numeric >= 50), .groups = "drop")
  msg("=== Age by group ==="); print(as.data.frame(age_tab))

  # ---------------------------------------------------------------------------
  # 2. Severity in cirrhosis / HCC
  # ---------------------------------------------------------------------------
  bili <- num(meta$Bili_T) / 17.1; cr <- num(meta$Creatinine) / 88.4; inr <- num(meta$INR); na <- num(meta$Sodium)
  fl <- function(x, lo = 1) pmax(x, lo)
  meld <- 10 * (0.957 * log(pmin(fl(cr), 4)) + 0.378 * log(fl(bili)) + 1.120 * log(fl(inr)) + 0.643)
  meld <- round(meld)
  nab  <- pmin(pmax(na, 125), 137)
  meldna <- ifelse(meld > 11, round(meld + 1.32 * (137 - nab) - 0.033 * meld * (137 - nab)), meld)
  meta$MELD <- meld; meta$MELD_Na <- meldna

  liv <- meta[meta$Status %in% c("Cirrhosis", "HCC"), ]
  liv$Status <- droplevels(liv$Status)
  sev_rows <- list(
    srow("MELD, median (IQR)",    tapply(liv$MELD,    liv$Status, med_iqr)),
    srow("MELD-Na, median (IQR)", tapply(liv$MELD_Na, liv$Status, med_iqr)),
    srow("MELD >= 15, n (%)",     tapply(liv$MELD >= 15, liv$Status, n_pct)),
    srow("INR available, n",      tapply(!is.na(num(liv$INR)), liv$Status, sum))
  )
  if (check_cols(liv, CHILD_PUGH_COL, "Child-Pugh")) {
    raw    <- toupper(trimws(as.character(liv[[CHILD_PUGH_COL]])))
    pat    <- "(^|[^A-Z])([ABC])([^A-Z]|$)"                     # stand-alone A/B/C ("A", "Class B", "CP-C")
    letter <- ifelse(grepl(pat, raw, perl = TRUE), sub(paste0(".*", pat, ".*"), "\\2", raw, perl = TRUE), NA)
    score  <- suppressWarnings(as.numeric(gsub("[^0-9]", "", raw)))   # numeric Child-Pugh score 5-15
    cp     <- ifelse(!is.na(letter), letter, as.character(cut(score, c(4, 6, 9, 15), labels = c("A", "B", "C"))))
    for (cl in c("A", "B", "C")) sev_rows[[length(sev_rows) + 1]] <-
      srow(paste0("Child-Pugh class ", cl, ", n (%)"), tapply(cp == cl, liv$Status, n_pct))
    sev_rows[[length(sev_rows) + 1]] <- srow("Child-Pugh missing, n", tapply(is.na(cp), liv$Status, sum))
  }
  for (v in c("Lactulose_Use", "Rifaximin_Use", "Diuretic", "BetaBlocker_Use")) {
    if (v %in% names(liv)) sev_rows[[length(sev_rows) + 1]] <-
      srow(paste0(v, " (decompensation proxy), n (%)"), tapply(yes(liv[[v]]), liv$Status, n_pct))
  }
  sev_tab <- bind_rows(sev_rows)
  msg("\n=== Liver disease severity (cirrhosis, HCC) ==="); print(sev_tab)

  # ---------------------------------------------------------------------------
  # 3. MASLD definition audit (2023 Delphi: steatosis + >= 1 of 5 cardiometabolic criteria)
  #    BP is not in the dataset -> antihypertensive use is the proxy.
  # ---------------------------------------------------------------------------
  male <- grepl("^m", tolower(as.character(meta$Gender)))
  dm_tx <- Reduce(`|`, lapply(intersect(c("Metformin", "Insulin", "Sulfonylurea_Use", "DPP4_Inhibitor_Use",
                                           "SGLT2_Inhibitor_Use", "GLP1_RA_Use"), names(meta)), function(v) yes(meta[[v]])))
  crit <- data.frame(
    C1_BMI_ge25 = num(meta$BMI) >= 25,                                     # waist circumference not recorded
    C2_HbA1c_ge5.7_or_T2D_tx = num(meta$Hgb_A1C) >= 5.7 | dm_tx,          # fasting glucose not recorded (Glu_R is random)
    C3_BP_proxy_antihypertensive = if ("antihypertensive" %in% names(meta)) yes(meta$antihypertensive) else NA,
    C4_TG_ge1.70_or_tx = num(meta$TG) >= 1.70 | yes(meta$Statin_Use),
    C5_HDL_low_or_tx = ifelse(male, num(meta$HDL) <= 1.0, num(meta$HDL) <= 1.3) | yes(meta$Statin_Use)
  )
  crit$N_criteria <- rowSums(crit, na.rm = TRUE)
  ms <- meta$Status == "MASLD"
  masld_tab <- data.frame(
    Criterion = c(names(crit)[1:5], ">= 1 criterion", ">= 2 criteria", "0 criteria"),
    MASLD_n_pct = c(vapply(1:5, function(j) n_pct(crit[ms, j]), ""),
                    n_pct(crit$N_criteria[ms] >= 1), n_pct(crit$N_criteria[ms] >= 2), n_pct(crit$N_criteria[ms] == 0)))
  msg("\n=== MASLD cardiometabolic criteria (MASLD group) ==="); print(masld_tab)

  fib_tab <- data.frame()
  for (v in c(STEATOSIS_COL, FIBROSIS_COL)) {
    if (check_cols(meta, v, "MASLD staging")) {
      t1 <- as.data.frame(table(Stage = meta[[v]][ms], useNA = "ifany")); t1$Variable <- v
      fib_tab <- bind_rows(fib_tab, t1)
    }
  }
  if (nrow(fib_tab)) { msg("\n=== Steatosis / fibrosis grades in MASLD ==="); print(fib_tab) }

  # Same criteria in controls (did any control meet cardiometabolic criteria?)
  ctrl_crit <- data.frame(Criterion = c(names(crit)[1:5], ">= 1 criterion"),
                          Healthy_n_pct = c(vapply(1:5, function(j) n_pct(crit[meta$Status == "Healthy", j]), ""),
                                            n_pct(crit$N_criteria[meta$Status == "Healthy"] >= 1)))

  # ---------------------------------------------------------------------------
  # 4. HCC chart-review template (fill manually, then re-run to summarise)
  # ---------------------------------------------------------------------------
  hcc_template <- data.frame(SampleID = meta$SampleID[meta$Status == "HCC"],
                             BCLC_stage = "", Tumour_number = "", Largest_nodule_cm = "",
                             Macrovascular_invasion = "", Extrahepatic_spread = "",
                             Treatment_before_sampling = "",   # naive / resection / ablation / TACE / systemic
                             Days_from_treatment_to_sample = "", Ascites = "", Hepatic_encephalopathy = "")
  cirr_template <- data.frame(SampleID = meta$SampleID[meta$Status == "Cirrhosis"],
                              Ever_decompensated = "", Ascites = "", Hepatic_encephalopathy = "",
                              Variceal_bleeding = "", Child_Pugh_score = "")

  # ---------------------------------------------------------------------------
  # 5. Table S1 fix: zero-exposure rows must not carry a p-value
  # ---------------------------------------------------------------------------
  safe_chisq_p <- function(x, g) {
    tab <- table(g, x)
    if (ncol(tab) < 2) return(NA_real_)      # no variation (e.g. antibiotics 0% everywhere) -> report "--"
    set.seed(SEED); chisq.test(tab, simulate.p.value = TRUE, B = 10000)$p.value
  }
  abx_col <- intersect(c("Antibiotic_Use_Last3mo"), names(meta))
  if (length(abx_col)) msg("\nAntibiotics row p-value (should be NA -> '--'): %s",
                           format(safe_chisq_p(meta[[abx_col]], meta$Status)))

  write_sheets(list(Age_by_group = age_tab, Liver_severity = sev_tab, MASLD_criteria = masld_tab,
                    Controls_criteria = ctrl_crit, MASLD_staging = fib_tab,
                    MELD_per_patient = meta[meta$Status %in% c("Cirrhosis", "HCC"), c("SampleID", "Status", "MELD", "MELD_Na")],
                    TEMPLATE_HCC_chart_review = hcc_template, TEMPLATE_Cirrhosis_chart = cirr_template),
               file.path(REV_DIR, "Rev06_Clinical_Characterisation.xlsx"))
})

# ==============================================================================
# 8. HCC_UNCERTAINTY
#    F/B and Streptococcus bootstrap CIs, leave-one-out (R2-17)
# ==============================================================================
if (isTRUE(RUN$HCC_UNCERTAINTY)) local({
  msg("\n################ SECTION 8: HCC_UNCERTAINTY ################")

  N_BOOT <- 2000

  meta <- load_meta()
  cts  <- load_genus_counts()
  common <- intersect(colnames(cts), meta$SampleID)
  cts <- cts[, common]; meta <- meta[common, ]
  rel <- sweep(cts, 2, colSums(cts), "/") * 100

  phylum_sum <- function(p) colSums(rel[grepl(paste0("p__", p), rownames(rel)), , drop = FALSE])
  firm <- phylum_sum("Firmicutes"); bact <- phylum_sum("Bacteroidota")
  if (all(firm == 0)) firm <- phylum_sum("Bacillota")          # SILVA naming fallback
  fb      <- firm / bact
  log2fb  <- log2((firm + 0.01) / (bact + 0.01))
  strep_r <- grepl("g__Streptococcus($|;)", rownames(cts))
  strep_c <- colSums(cts[strep_r, , drop = FALSE]); strep_rel <- colSums(rel[strep_r, , drop = FALSE])

  boot_ci <- function(x, f = median) {
    set.seed(SEED); b <- replicate(N_BOOT, f(sample(x, replace = TRUE)))
    quantile(b[is.finite(b)], c(.025, .975), na.rm = TRUE)
  }

  rows <- lapply(GROUPS, function(g) {
    s <- meta$Status == g; n <- sum(s)
    ci_fb <- boot_ci(fb[s]); ci_lfb <- boot_ci(log2fb[s]); ci_st <- boot_ci(strep_rel[s])
    det <- sum(strep_c[s] > 0); wl <- wilson_ci(det, n)
    data.frame(Group = g, n = n,
               FB_median = median(fb[s]), FB_CI = sprintf("%.2f-%.2f", ci_fb[1], ci_fb[2]),
               Zero_Bacteroidota_n = sum(bact[s] == 0),
               log2FB_median = median(log2fb[s]), log2FB_CI = sprintf("%.2f-%.2f", ci_lfb[1], ci_lfb[2]),
               Strep_detected = sprintf("%d/%d (%.0f%%; 95%% CI %.0f-%.0f%%)", det, n, 100 * det / n, 100 * wl[1], 100 * wl[2]),
               Strep_median_relab_pct = median(strep_rel[s]),
               Strep_relab_CI = sprintf("%.3f-%.3f", ci_st[1], ci_st[2]),
               Strep_mean_relab_pct = mean(strep_rel[s]))
  })
  hcc_tab <- bind_rows(rows)
  msg("=== F/B and Streptococcus by group ==="); print(hcc_tab)

  # LOO stability: CLR of Streptococcus (pseudocount 1), HCC mean - Healthy mean
  clr_all <- { l <- log(t(cts) + 1); sweep(l, 1, rowMeans(l)) }
  strep_clr <- rowSums(clr_all[, strep_r, drop = FALSE])
  h <- strep_clr[meta$Status == "Healthy"]; hc <- strep_clr[meta$Status == "HCC"]
  full_diff <- (mean(hc) - mean(h)) / log(2)
  loo <- vapply(seq_along(hc), function(i) (mean(hc[-i]) - mean(h)) / log(2), 0)
  wt  <- wilcox.test(hc, h, exact = FALSE)
  loo_tab <- data.frame(Statistic = c("Full CLR difference (log2 units)", "LOO minimum", "LOO maximum",
                                      "Wilcoxon p (HCC vs Healthy, CLR)", "Most influential HCC sample"),
                        Value = c(round(full_diff, 2), round(min(loo), 2), round(max(loo), 2),
                                  signif(wt$p.value, 3), names(hc)[which.min(loo)]))
  msg("\n=== Streptococcus HCC vs Healthy: leave-one-out ==="); print(loo_tab)

  write_sheets(list(FB_Streptococcus_by_group = hcc_tab, Streptococcus_LOO = loo_tab,
                    Per_sample = data.frame(SampleID = names(fb), Status = meta$Status, Firmicutes_pct = firm,
                                            Bacteroidota_pct = bact, FB = fb, Streptococcus_pct = strep_rel)),
               file.path(REV_DIR, "Rev07_HCC_Uncertainty.xlsx"))
})

msg("\nAll requested sections finished. Workbooks in: %s", REV_DIR)
