library(dplyr)
library(ggplot2)
library(jsonlite)

# plot_bench_find_se.R for every finder: one plot per find_*() function, from
# the timings bench_finders.R wrote to benchmarks/results/finders/. A finder
# with a `type` argument (find_se / find_ie / find_cse / find_cie / find_mxe)
# gets one line per type.
# Run from the repo root. Writes, to benchmarks/results/finders/:
#   bench_find_<event>.png   one plot per finder
#   bench_finders_plots.rds  the plot objects, named by event, for the
#                            combined figure (figures/plot_figure_s1.R)

res_dir <- "benchmarks/results/finders"

# the finders, in the order they are shown in the combined figure
finder_fns <- c(
  se   = "find_se",
  ie   = "find_ie",
  cse  = "find_cse",
  cie  = "find_cie",
  mxe  = "find_mxe",
  ri   = "find_ri",
  si   = "find_si",
  a5ss = "find_a5ss",
  a3ss = "find_a3ss",
  atss = "find_atss",
  ates = "find_ates"
)

panel_labels <- c(
  n_genes        = "Genes",
  event_ratio    = "Event ratio",
  n_pos_per_gene = "Pos tx / gene"
)

# one colour and one point shape per `type`. "over" is the blue of figure 1
# (the A nucleotide colour of slPlotFlankProfilesByNt()); a finder without a
# `type` argument is drawn in grey, so it is not read as one of the three, and
# gets its own legend entry. Every plot carries all four levels, so the eleven
# legends are identical and figure S1 collects them into one.
type_cols   <- c(over = "#3A86FF", boundary = "#FB5607", identical = "#8338EC",
                 none = "grey30")
type_shapes <- c(over = 16, boundary = 17, identical = 15, none = 16)
type_labels <- c(over = "over", boundary = "boundary", identical = "identical",
                 none = "no type argument")

# finders drawn on the same y axis, so their runtimes can be compared directly:
# the transcription start and end sites
shared_y <- list(c("atss", "ates"))

# the timings of one finder along the three sweeps, in ms, one row per point
finder_data <- function(event) {
  results <- read.csv(file.path(res_dir, sprintf("bench_find_%s.csv", event)))
  cfg <- jsonlite::read_json(
    file.path(res_dir, sprintf("grid_config_find_%s.json", event))
  )
  BL <- cfg$baseline

  bind_rows(
    results |>
      filter(event_ratio == BL$event_ratio, n_pos_per_gene == BL$n_pos_per_gene) |>
      mutate(focal_dim = "n_genes", focal_value = n_genes),
    results |>
      filter(n_genes == BL$n_genes, n_pos_per_gene == BL$n_pos_per_gene) |>
      mutate(focal_dim = "event_ratio", focal_value = event_ratio),
    results |>
      filter(n_genes == BL$n_genes, event_ratio == BL$event_ratio) |>
      mutate(focal_dim = "n_pos_per_gene", focal_value = n_pos_per_gene)
  ) |>
    mutate(
      focal_dim = factor(focal_dim, levels = names(panel_labels)),
      # NA for a finder without a `type` argument
      type      = factor(ifelse(is.na(type), "none", type),
                         levels = names(type_cols)),
      mean_ms   = mean_s * 1000,
      se_ms     = (sd_s * 1000) / sqrt(n_itr),
      lo        = mean_ms - se_ms,
      hi        = mean_ms + se_ms
    ) |>
    # thin the Pos tx / gene panel to 4 points across the swept range
    filter(!(focal_dim == "n_pos_per_gene" &
               !focal_value %in% c(1, 4, 7, 10)))
}

# `y_range`: the lowest and highest values the y axis must show; by default
# the finder's own error bars, or those of every finder sharing its axis
plot_finder <- function(event, plot_data = finder_data(event),
                        y_range = c(min(plot_data$lo), max(plot_data$hi))) {
  # y axis from just below the lowest error bar (to the step below it, so no
  # empty band) to the highest error bar. Labelled every 25 ms, or at the
  # first wider step that leaves no more than 6 labels. An error bar less than
  # a tenth of a step below a label starts the axis at that label, not a whole
  # step lower; the bar still shows, in the axis padding (coord_cartesian
  # below zooms rather than dropping data outside the limits).
  y_steps <- c(25, 50, 100, 200, 500, 1000)
  y_step  <- y_steps[diff(y_range) / y_steps <= 6][1]
  y_lim   <- c(floor(y_range[1] / y_step + 0.1) * y_step, y_range[2])
  y_breaks <- seq(y_lim[1], y_lim[2], by = y_step)

  # one line per type, or a single grey line ("none") when the finder has no
  # type. `limits` keeps all four levels in every legend, present or not, and
  # `show.legend = TRUE` draws a key for each of them (ggplot2 >= 3.5 leaves
  # out the keys of levels absent from a layer's data otherwise), so the eleven
  # legends match and figure S1 merges them into one
  ggplot(plot_data, aes(focal_value, mean_ms, colour = type, shape = type)) +
    # spread across the iterations (gc iterations kept, see bench_finders.R):
    # mean +/- SE, as error bars
    geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.6,
                   show.legend = TRUE) +
    geom_line(linewidth = 1.1, show.legend = TRUE) +
    geom_point(size = 2.6, show.legend = TRUE) +
    scale_colour_manual(values = type_cols, limits = names(type_cols),
                        labels = type_labels, name = "type") +
    scale_shape_manual(values = type_shapes, limits = names(type_cols),
                       labels = type_labels, name = "type") +
    facet_wrap(
      ~ focal_dim,
      scales = "free_x",
      strip.position = "bottom",
      labeller = labeller(focal_dim = panel_labels)
    ) +
    # 1000 -> 1k, 2000 -> 2k (the Genes panel); small values stay as they are
    scale_x_continuous(
      # about three breaks per panel, few enough that the labels do not run
      # together when eleven of these plots share the combined figure
      breaks = \(l) pretty(l, n = 3),
      labels = \(x) ifelse(x >= 1000, paste0(x / 1000, "k"), x)
    ) +
    scale_y_continuous(breaks = y_breaks) +
    coord_cartesian(ylim = y_lim) +
    labs(title = paste0(finder_fns[[event]], "()"), x = NULL,
         y = "Mean runtime (ms)") +
    theme_classic(base_size = 14) +
    theme(
      plot.title         = element_text(size = 19, margin = margin(b = 10)),
      axis.line          = element_line(colour = "black", linewidth = 0.5),
      axis.ticks         = element_line(colour = "black", linewidth = 0.5),
      axis.ticks.length  = unit(3.5, "pt"),
      axis.text          = element_text(colour = "black", size = 16),
      axis.text.x        = element_text(margin = margin(t = 6)),
      axis.text.y        = element_text(margin = margin(r = 4)),
      axis.title         = element_text(size = 19),
      axis.title.y       = element_text(margin = margin(r = 14)),
      panel.grid.major.y = element_line(colour = "grey85", linewidth = 0.3),
      strip.background   = element_blank(),
      strip.placement    = "outside",
      # the strips stand in for the x axis title, so styled like one (plain,
      # 19pt, black), as in figure 1
      strip.switch.pad.wrap = unit(8, "pt"),
      strip.text         = element_text(colour = "black", size = 19,
                                        margin = margin(t = 6, b = 6)),
      # a strip title wider than its panel ("Pos tx / gene") is not cut off
      strip.clip         = "off",
      panel.spacing.x    = unit(1.5, "lines"),
      legend.position    = "bottom",
      legend.title       = element_text(size = 19),
      legend.text        = element_text(size = 16),
      legend.key.width   = unit(2.2, "lines")
    )
}

# only the finders bench_finders.R has been run for
have <- file.exists(file.path(res_dir,
                              sprintf("bench_find_%s.csv", names(finder_fns))))
if (!any(have)) {
  stop("No bench_find_<event>.csv in ", res_dir,
       "; run benchmarks/bench_finders.R first.")
}
if (!all(have)) {
  message("No results yet for: ",
          paste(names(finder_fns)[!have], collapse = ", "))
}

events <- names(finder_fns)[have]
data <- lapply(setNames(nm = events), finder_data)

# each finder's y range, widened to cover the finders it shares an axis with
y_ranges <- lapply(data, \(d) c(min(d$lo), max(d$hi)))
for (group in shared_y) {
  group <- intersect(group, events)
  if (length(group) < 2) next
  rng <- range(unlist(y_ranges[group]))
  y_ranges[group] <- list(rng)
}

plots <- lapply(setNames(nm = events),
                \(e) plot_finder(e, data[[e]], y_ranges[[e]]))

for (event in names(plots)) {
  ggsave(file.path(res_dir, sprintf("bench_find_%s.png", event)),
         plots[[event]], width = 8, height = 7, dpi = 300)
}
# the plot objects too, for the combined figure (figures/plot_figure_s1.R)
saveRDS(plots, file.path(res_dir, "bench_finders_plots.rds"))
message("Saved -> ", res_dir, "/bench_find_<event>.png for ",
        paste(names(plots), collapse = ", "), " (+ bench_finders_plots.rds)")
