#!/usr/bin/env Rscript
# ============================================================
# Figure 1: Study Design & Analytical Workflow
# Paper 1 - 16S rRNA Gut Microbiome in Saudi Liver Disease Cohort
# ============================================================

library(grid)

# ------------------------------------------------------------
# Working directory = folder containing this script (relative paths)
# Works with RStudio "Source", source("..."), and Rscript.
# ------------------------------------------------------------
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
OUT_NAME <- "Figure1_Study_Design"

# ============================================================
# CONSISTENT COLOR PALETTE (use across ALL figures)
# ============================================================
group_colors <- c(
  "Healthy"   = "#4DAF4A",
  "MASLD"     = "#377EB8",
  "Cirrhosis" = "#984EA3",
  "HCC"       = "#E41A1C"
)

# Supporting colors
col_dark   <- "#333333"
col_med    <- "#555555"
col_light  <- "#999999"
col_bg     <- "#F5F5F5"
col_arrow  <- "#555555"
col_white  <- "#FFFFFF"

main_font <- "sans"

# ============================================================
# HELPER FUNCTIONS
# ============================================================
draw_rounded_box <- function(x, y, w, h, fill, border, lwd = 1.5, r = unit(2, "mm")) {
  grid.roundrect(x = x, y = y, width = w, height = h, r = r,
                 gp = gpar(fill = fill, col = border, lwd = lwd), just = "centre")
}

draw_text <- function(label, x, y, size = 8, col = col_dark,
                      fontface = "plain", just = "centre", fontfamily = main_font) {
  grid.text(label, x = x, y = y,
            gp = gpar(fontsize = size, col = col, fontface = fontface, fontfamily = fontfamily),
            just = just)
}

draw_arrow_down <- function(x, y1, y2, col = col_arrow, lwd = 1.5) {
  grid.lines(x = c(x, x), y = c(y1, y2), gp = gpar(col = col, lwd = lwd),
             arrow = arrow(length = unit(2, "mm"), type = "closed"), default.units = "npc")
}

draw_branch <- function(x1, y1, x2, y2, col = col_arrow, lwd = 0.8) {
  grid.lines(x = c(x1, x2), y = c(y1, y2), gp = gpar(col = col, lwd = lwd))
}

draw_number_circle <- function(x, y, num, bg_col, r = unit(2.2, "mm")) {
  grid.circle(x = x, y = y, r = r, gp = gpar(fill = bg_col, col = NA))
  grid.text(as.character(num), x = x, y = y,
            gp = gpar(fontsize = 6.5, col = col_white, fontface = "bold", fontfamily = main_font))
}

# ============================================================
# OUTPUT: 7 x 5 inches (full-width journal figure)
# ============================================================
draw_figure1 <- function() {
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2, widths = c(0.47, 0.53))))

# ====================
# PANEL A - Study Cohort
# ====================
pushViewport(viewport(layout.pos.col = 1))

draw_text("A", x = 0.04, y = 0.96, size = 14, fontface = "bold", just = c("left", "top"))
draw_text("Study Cohort", x = 0.52, y = 0.94, size = 10, fontface = "bold")

# Top box
draw_rounded_box(x = 0.52, y = 0.82, w = 0.78, h = 0.10,
                 fill = "#F5F5F5", border = "#888888", lwd = 1.5)
draw_text("Saudi Arabian Liver Disease Cohort", x = 0.52, y = 0.835, size = 8, fontface = "bold")
draw_text("N = 205 participants", x = 0.52, y = 0.805, size = 7.5, col = col_med)

draw_arrow_down(0.52, 0.765, 0.72)

# Fecal collection
draw_rounded_box(x = 0.52, y = 0.685, w = 0.62, h = 0.06,
                 fill = col_bg, border = "#BBBBBB", lwd = 0.8)
draw_text("Fecal sample collection", x = 0.52, y = 0.685, size = 7.5)

draw_arrow_down(0.52, 0.65, 0.605)

# Sequencing
draw_rounded_box(x = 0.52, y = 0.57, w = 0.78, h = 0.06,
                 fill = col_bg, border = "#BBBBBB", lwd = 0.8)
draw_text("16S rRNA gene sequencing (V4, Illumina MiSeq)",
          x = 0.52, y = 0.57, size = 7)

draw_arrow_down(0.52, 0.535, 0.49)

# Classification label
draw_text("Disease Group Classification", x = 0.52, y = 0.47,
          size = 7.5, col = col_med, fontface = "italic")

# Branching lines
x_groups <- c(0.13, 0.37, 0.63, 0.87)
for (xg in x_groups) {
  draw_branch(0.52, 0.45, xg, 0.41)
  draw_branch(xg, 0.41, xg, 0.365)
}

# Disease group boxes
group_names <- c("Healthy", "MASLD", "Cirrhosis", "HCC")
group_n     <- c("n = 70", "n = 59", "n = 59", "n = 17")
box_w <- 0.21
box_h <- 0.12
y_center <- 0.295

for (i in 1:4) {
  draw_rounded_box(x = x_groups[i], y = y_center, w = box_w, h = box_h,
                   fill = group_colors[group_names[i]], border = "#333333", lwd = 1.2,
                   r = unit(2, "mm"))
  draw_text(group_names[i], x = x_groups[i], y = y_center + 0.02,
            size = 8.5, fontface = "bold", col = col_white)
  draw_text(group_n[i], x = x_groups[i], y = y_center - 0.025,
            size = 7.5, col = col_white)
}


popViewport()              

# ====================
# PANEL B - Analytical Workflow
# ====================
pushViewport(viewport(layout.pos.col = 2))

draw_text("B", x = 0.03, y = 0.96, size = 14, fontface = "bold", just = c("left", "top"))
draw_text("Analytical Workflow", x = 0.50, y = 0.94, size = 10, fontface = "bold")

# Colored workflow boxes (Stage 1 = teal to avoid MASLD blue clash)
stages <- list(
  list(title = "Data Processing",
       color = "#E0F2F1", edge = "#00897B",
       items = c("QIIME 2 pipeline (primer removal, cutadapt)",
                 "DADA2 denoising",
                 "SILVA 138 (99%) taxonomy",
                 "Rarefaction (13,000 reads/sample)")),
  list(title = "Diversity Analysis",
       color = "#E8F5E9", edge = "#4CAF50",
       items = c("Alpha diversity (Chao1, Observed, Shannon, Simpson, Faith's PD)",
                 "Beta diversity (weighted/unweighted UniFrac, Bray-Curtis, Jaccard)")),
  list(title = "Differential Abundance",
       color = "#FFF8E1", edge = "#FF9800",
       items = c("ANCOM-BC2 (primary)",
                 "DESeq2 (complementary)",
                 "FDR correction (q < 0.05)")),
  list(title = "Clinical Integration",
       color = "#FCE4EC", edge = "#E91E63",
       items = c("Spearman correlations (taxa x clinical markers)",
                 "Confounder-adjusted PERMANOVA (11 covariates)",
                 "Composite Dysbiosis Index"))
)

y_tops    <- c(0.88, 0.665, 0.46, 0.235)
box_h_vec <- c(0.175, 0.145, 0.155, 0.175)
box_w_all <- 0.88
x_center  <- 0.50

for (i in seq_along(stages)) {
  s  <- stages[[i]]
  yt <- y_tops[i]
  bh <- box_h_vec[i]
  yc <- yt - bh / 2

  draw_rounded_box(x = x_center, y = yc, w = box_w_all, h = bh,
                   fill = s$color, border = s$edge, lwd = 1.2)

  draw_number_circle(x = x_center - box_w_all/2 + 0.04,
                     y = yt - 0.025, num = i, bg_col = s$edge)

  draw_text(s$title, x = x_center - box_w_all/2 + 0.08,
            y = yt - 0.025, size = 8, fontface = "bold", just = "left")

  item_y <- yt - 0.055
  for (item in s$items) {
    draw_text(paste0("-  ", item),
              x = x_center - box_w_all/2 + 0.055,
              y = item_y, size = 6.5, col = col_med, just = "left")
    item_y <- item_y - 0.028
  }

  if (i < length(stages)) {
    arrow_start <- yc - bh/2 - 0.005
    arrow_end   <- y_tops[i + 1] + 0.005
    draw_arrow_down(x_center, arrow_start, arrow_end)
  }
}



popViewport()
}

# ============================================================
# SAVE: PDF (vector), PNG (600 dpi), SVG (vector)
# ============================================================
cairo_pdf(paste0(OUT_NAME, ".pdf"), width = 7, height = 5, family = main_font)
draw_figure1(); invisible(dev.off())

png(paste0(OUT_NAME, ".png"), width = 7, height = 5, units = "in", res = 600, type = "cairo")
draw_figure1(); invisible(dev.off())

svg(paste0(OUT_NAME, ".svg"), width = 7, height = 5, family = main_font)
draw_figure1(); invisible(dev.off())

cat("Figure 1 saved:", paste0(OUT_NAME, c(".pdf", ".png", ".svg")), sep = "\n  ")
cat("\nLocation:", getwd(), "\n")
