# Supplementary figure 1: runtime of every find_*() function, one plot per
# finder in a 4 x 3 grid, each drawn like the find_se() panel of figure 1.
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

# rows 1-2: the five finders with a `type` argument (skipped / included
# exons, mutually exclusive exons, consecutive skipped / included exons),
# then retained introns; row 3: spliced introns, then the 5' and 3' splice
# sites; row 4: the transcription start and end sites, and an empty cell
grid_order <- c("se", "ie", "mxe",
                "cse", "cie", "ri",
                "si", "a5ss", "a3ss",
                "atss", "ates")
missing <- setdiff(grid_order, names(plots))
if (length(missing) > 0) {
  stop("No plot for: ", paste(missing, collapse = ", "),
       "; run benchmarks/bench_finders.R for them, then ",
       "benchmarks/plot_bench_finders.R.")
}

fig <- wrap_plots(plots[grid_order], nrow = 4, ncol = 3) +
  # the three `type` legends are the same, so they are drawn once, below the
  # grid; the y axis title is drawn once per row
  plot_layout(guides = "collect", axis_titles = "collect_y") &
  theme(legend.position = "bottom",
        plot.margin = margin(t = 10, r = 14, b = 10, l = 6))

# a fourth row at the same row height as the 3 x 3 layout had
ggsave(here("figures/figure_s1.png"), fig, width = 19.5, height = 24.5, dpi = 200)
message("Saved -> figures/figure_s1.png")
