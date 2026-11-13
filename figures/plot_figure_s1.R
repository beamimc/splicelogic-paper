# Supplementary figure 1: runtime of every find_*() function, one plot per
# finder in a 3 x 3 grid, each drawn like the find_se() panel of figure 1.
# Run from the repo root, after the two scripts that time the finders and save
# the plot objects read below:
#   Rscript benchmarks/bench_finders.R
#   Rscript benchmarks/plot_bench_finders.R
# Output: figures/figure_s1.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(here)
})

# one plot per finder, named by event; saved by benchmarks/plot_bench_finders.R
plots <- readRDS(here("benchmarks/results/finders/bench_finders_plots.rds"))

# row 1: the three finders with a `type` argument; row 2: retained and spliced
# introns, then the 5' splice site; row 3: the 3' splice site, then the
# transcription start and end sites
grid_order <- c("se", "ie", "mxe", "ri", "si", "a5ss", "a3ss", "atss", "ates")
missing <- setdiff(grid_order, names(plots))
if (length(missing) > 0) {
  stop("No plot for: ", paste(missing, collapse = ", "),
       "; run benchmarks/bench_finders.R for them, then ",
       "benchmarks/plot_bench_finders.R.")
}

fig <- wrap_plots(plots[grid_order], nrow = 3, ncol = 3) +
  # the three `type` legends are the same, so they are drawn once, below the
  # grid; the y axis title is drawn once per row
  plot_layout(guides = "collect", axis_titles = "collect_y") &
  theme(legend.position = "bottom",
        plot.margin = margin(t = 10, r = 14, b = 10, l = 6))

ggsave(here("figures/figure_s1.png"), fig, width = 19.5, height = 18.5, dpi = 200)
message("Saved -> figures/figure_s1.png")
