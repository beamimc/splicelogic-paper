library(dplyr)
library(ggplot2)
library(jsonlite)

results <- read.csv("benchmarks/results/bench_find_se.csv")
cfg     <- jsonlite::read_json("benchmarks/results/grid_config.json")
BL      <- cfg$baseline

panel_labels <- c(
  n_genes        = "Genes",
  event_ratio    = "Event ratio",
  n_pos_per_gene = "Pos tx / gene"
)

# same blue as slPlotFlankProfilesByNt() (the A nucleotide colour)
line_col <- "#3A86FF"

plot_data <- bind_rows(
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
    mean_ms   = mean_s * 1000,
    se_ms     = (sd_s * 1000) / sqrt(n_itr),
    lo        = mean_ms - se_ms,
    hi        = mean_ms + se_ms
  ) |>
  # thin the Pos tx / gene panel to 4 points across the swept range
  filter(!(focal_dim == "n_pos_per_gene" &
             !focal_value %in% c(1, 4, 7, 10))) |>
  # figure 1 shows the default type only
  filter(type == "over")

# y axis from just below the lowest error bar (to the 25 ms below it, so no
# empty band) to the highest error bar, labelled every 25 ms
y_step <- 25
y_lim <- c(floor(min(plot_data$lo) / 25) * 25, max(plot_data$hi))
y_breaks <- seq(ceiling(y_lim[1] / y_step) * y_step, y_lim[2], by = y_step)

p <- ggplot(plot_data, aes(focal_value, mean_ms)) +
  # spread across the iterations (gc iterations kept, see bench_find_se.R):
  # mean +/- SE, as error bars
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.6,
                 colour = line_col) +
  # as thick as the lines of the flank profiles (slPlotFlankProfilesByNt())
  geom_line(linewidth = 1.1, colour = line_col) +
  geom_point(size = 2, colour = line_col) +
  facet_wrap(
    ~ focal_dim,
    scales = "free_x",
    strip.position = "bottom",
    labeller = labeller(focal_dim = panel_labels)
  ) +
  # 1000 -> 1k, 2000 -> 2k (the Genes panel); small values stay as they are
  scale_x_continuous(
    # about five breaks per panel, few enough that the labels do not run
    # together when the plot is drawn in figure 1
    breaks = \(l) pretty(l, n = 5),
    labels = \(x) ifelse(x >= 1000, paste0(x / 1000, "k"), x)
  ) +
  scale_y_continuous(breaks = y_breaks, limits = y_lim) +
  labs(x = NULL, y = "Mean runtime (ms)") +
  theme_classic(base_size = 14) +
  theme(
    axis.line         = element_line(colour = "black", linewidth = 0.5),
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
    # 19pt, black) and 14pt below the tick labels, as in the flank figures
    strip.switch.pad.wrap = unit(8, "pt"),
    strip.text         = element_text(colour = "black", size = 19,
                                      margin = margin(t = 6, b = 6)),
    # a strip title wider than its panel ("Pos tx / gene") is not cut off
    strip.clip         = "off",
    panel.spacing.x    = unit(1.5, "lines")
  )

ggsave("benchmarks/results/bench_find_se.png", p, width = 8, height = 7, dpi = 300)
# the plot object too, for the combined figure (figures/plot_figure1.R)
saveRDS(p, "benchmarks/results/bench_find_se.rds")
message("Saved -> benchmarks/results/bench_find_se.png (+ .rds)")
