#!/usr/bin/env Rscript
# Alpha diversity effect sizes vs healthy controls (Reviewer 2, comment 15)
# % difference in medians (bootstrap 95% CI), Cliff's delta, KW + BH across 5 metrics,
# and minimum detectable effect (80% power, two-sided alpha 0.05).
# Run from Desktop/Reviewers/Figure_2 (finds Desktop/16S_reprocess automatically).
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
A <- file.path(find_proj_dir(), "qiime2_analysis_v2", "Alpha_diversity")
files <- c("Observed features" = "Observed feature/Observed_feature_metadata.tsv", Chao1 = "Chao1/Chao1.tsv",
           Shannon = "Shannon/Shannon_Liver.tsv", Simpson = "Simpson/Simpson_liver_micorbiome.tsv",
           "Faith's PD" = "Faith-pd/faith-pd.tsv")
GROUPS <- c("MASLD", "Cirrhosis", "HCC"); N_BOOT <- 2000; set.seed(2024)
cliff <- function(x, y) mean(sign(outer(x, y, "-")))
res <- list(); kw <- c()
for (m in names(files)) {
  d <- read.delim(file.path(A, files[m]), check.names = FALSE); d <- d[!grepl("^#", d[[1]]), ]
  v <- as.numeric(d[[ncol(d)]]); g <- sub("^NAFLD$", "MASLD", d$Status)
  kw[m] <- kruskal.test(v ~ factor(g))$p.value
  h <- v[g == "Healthy"]
  for (G in GROUPS) {
    x <- v[g == G]
    pct <- function(a, b) 100 * (median(a) / median(b) - 1)
    bs <- replicate(N_BOOT, pct(sample(x, replace = TRUE), sample(h, replace = TRUE)))
    res[[length(res) + 1]] <- data.frame(Metric = m, Group = G, Median_group = median(x), Median_healthy = median(h),
      Pct_diff = pct(x, h), CI_low = quantile(bs, 0.025), CI_high = quantile(bs, 0.975),
      Cliffs_delta = cliff(x, h), Wilcoxon_p = wilcox.test(x, h, exact = FALSE)$p.value, row.names = NULL)
  }
}
res <- do.call(rbind, res)
kw_tab <- data.frame(Metric = names(kw), KW_p = kw, KW_q_BH = p.adjust(kw, "BH"), row.names = NULL)
mde <- sapply(c(MASLD = 59, Cirrhosis = 59, HCC = 17), function(n) {
  f <- function(dd) { se <- sqrt(1/70 + 1/n); pnorm(dd/se - qnorm(0.975)) - 0.8 }; uniroot(f, c(0.01, 3))$root })
mde_tab <- data.frame(Group = names(mde), n = c(59, 59, 17), MDE_Cohens_d = mde,
                      MDE_Cliffs_delta = 2 * pnorm(mde / sqrt(2)) - 1, row.names = NULL)
write.csv(res, "Alpha_effect_sizes_vs_Healthy.csv", row.names = FALSE)
write.csv(kw_tab, "Alpha_KW_BH.csv", row.names = FALSE)
write.csv(mde_tab, "Alpha_minimum_detectable_effect.csv", row.names = FALSE)
options(width = 200); print(kw_tab, digits = 3); print(res, digits = 3); print(mde_tab, digits = 3)
