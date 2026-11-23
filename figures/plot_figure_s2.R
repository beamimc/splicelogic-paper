# Supplementary figure 2: one gene of the Glinos PTBP1 DTU transcripts drawn
# by Swan (a, its gene graph) and by splicelogic (b, plot_transcripts()), a on
# top of b.
# Run from the repo root, after benchmarks/swan/gene_match_glinos.Rmd, which
# writes the two GTFs Swan reads here. Swan runs in the conda env `swan`
# (setup at the top of gene_match_glinos.Rmd).
# Output: figures/figure_s2.png

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
  library(here)
  library(reticulate)
})
use_condaenv("swan", required = TRUE)
# splicelogic devel, loaded from the sibling repo
devtools::load_all(here("../splicelogic"), quiet = TRUE)

# the gene: on the plus strand (Swan calls nothing on the minus strand), a
# few DTU transcripts at FDR < 0.05 (ILK: 2 up, 1 down), and the two tools
# agree in place: Swan's exon skipping intron (an up transcript) runs across
# the exon splicelogic calls an SE for that same up transcript
gene_id <- "ENSG00000166333"  # ILK
fdr <- 0.05                   # as in gene_match_glinos.Rmd

no_version <- function(x) sub("\\..*", "", x)

# b: splicelogic ------------------------------------------------------------
# the gene's DTU transcripts at FDR < 0.05, their exons, and every event
sig_res <- readRDS(here("workflows/data/sig_res.rds"))
sig_res <- sig_res[sig_res$regular_FDR < fdr &
                     no_version(sig_res$gene_id) == gene_id, ]
txdb <- suppressPackageStartupMessages(
  AnnotationDbi::loadDb(here("workflows/data/gencode.v48.sqlite")))
exons <- suppressMessages(prepare_exons(txdb, sig_res, coef_col = "estimates")) |>
  preprocess(coef_col = "estimates")
events <- find_all_events(exons, verbose = FALSE)

# colours for this figure, passed in rather than the package defaults: royal
# blue up and sky blue down transcripts, and event colours chosen with them so
# the thirteen stay apart under protanopia, deuteranopia and tritanopia
# (smallest CIEDE2000 distance 10). SE, MXE, RI, A3SS and ATES are Okabe-Ito
# colours; the others keep off the two blues.
up_color   <- "#1F5FBF"
down_color <- "#7EAFD9"
event_colors <- c(se = "#D55E00", ie = "#B35693", cse = "#8F3A00",
                  cie = "#6F4177", mxe = "#CC79A7", ri = "#009E73",
                  si = "#00603D", a5ss = "#F6B5C1", a3ss = "#F0E442",
                  atss = "#A5C371", ates = "#E69F00")

# text sized for the page: label_size sets the event labels, the theme the
# legend and the x axis
panel_b <- plot_transcripts(events, exons, up_color = up_color,
                            down_color = down_color,
                            event_colors = event_colors, label_size = 4.2) +
  theme(legend.text     = element_text(size = 14),
        legend.key.size = unit(6, "mm"),
        axis.text.x     = element_text(size = 13),
        axis.title.x    = element_text(size = 14),
        # room for the last x-axis label, which sits on the panel edge
        plot.margin     = margin(5.5, 20, 5.5, 5.5))
# plot_transcripts() draws the transcript labels (id and estimate) at a fixed
# size, so that layer is resized here
for (i in seq_along(panel_b$layers)) {
  layer <- panel_b$layers[[i]]
  if (inherits(layer$geom, "GeomText") &&
      identical(rlang::as_label(layer$mapping$label), "tx_label")) {
    panel_b$layers[[i]]$aes_params$size <- 4.6
  }
}

# a: Swan -------------------------------------------------------------------
# the same inputs as gene_match_glinos.Rmd: the down transcripts are Swan's
# annotation, the up transcripts its transcriptome, so the dashed edges are
# the ones no down transcript has
gtf_dir <- here("benchmarks/swan/results/gene_match_fdr0.05")
gtfs <- c(down = file.path(gtf_dir, "down.gtf"),
          up   = file.path(gtf_dir, "up.gtf"))
if (!all(file.exists(gtfs))) {
  stop("No GTFs in ", gtf_dir, "; knit benchmarks/swan/gene_match_glinos.Rmd first.")
}
mpl <- import("matplotlib")
mpl$use("Agg")
plt <- import("matplotlib.pyplot")
swan <- import("swan_vis")

sg <- swan$SwanGraph()
invisible(sg$add_annotation(gtfs[["down"]]))
invisible(sg$add_transcriptome(gtfs[["up"]]))
sg$plot_graph(gene_id, indicate_novel = TRUE)
swan_png <- tempfile(fileext = ".png")
plt$savefig(swan_png, dpi = 300, bbox_inches = "tight")
plt$close("all")

panel_a <- wrap_elements(full = grid::rasterGrob(png::readPNG(swan_png)))

# figure ------------------------------------------------------------------
fig <- panel_a / panel_b +
  plot_layout(heights = c(1, 2.2)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 18, face = "bold"))

ggsave(here("figures/figure_s2.png"), fig, width = 12, height = 9, dpi = 300,
       bg = "white")
message("Saved -> figures/figure_s2.png")
