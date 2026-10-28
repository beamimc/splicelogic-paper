# Paper figure 1: find_se() runtime (a) next to the SE vs non-SE flank
# profiles (b), one row, each plot half the width.
# Run from the repo root, after the two scripts that save the plot objects
# read below. Output: figures/figure1.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(here)
})

# panel a: saved by benchmarks/plot_bench_find_se.R, which plots the timings in
# benchmarks/results/bench_find_se.csv (run benchmarks/bench_find_se.R first
# to create or refresh that CSV)
bench <- readRDS(here("benchmarks/results/bench_find_se.rds"))
# panel b: saved by figures/plot_flank_profiles_by_nt.R (SE vs non-SE, from
# the RBP workflow's data in workflows/data/)
flank <- readRDS(here("figures/flank_profiles_se_vs_non_se.rds"))

# both styled with 16pt tick and legend labels and 19pt titles; patchwork
# lines the panels up. Lowercase tags, so they are not read as the A
# nucleotide panel
fig <- bench + flank +
  # a and b half the width each
  plot_layout(ncol = 2, widths = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 24, face = "bold"))

# drawn on a canvas `scale` times the 18.6 x 8 in, so text and lines come out
# larger relative to the figure; dpi raised to keep the same pixel size
scale <- 0.8
ggsave(here("figures/figure1.png"), fig, width = 18.6, height = 8,
       scale = scale, dpi = 300 / scale)
message("Saved -> figures/figure1.png")
