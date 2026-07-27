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
    median_ms = median * 1000
  )

p <- ggplot(plot_data, aes(x = focal_value, y = median_ms, color = type)) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(
    ~ focal_dim,
    scales = "free_x",
    labeller = labeller(focal_dim = panel_labels)
  ) +
  scale_color_brewer(palette = "Set2", name = "Match type") +
  labs(
    title    = "find_se() runtime by gene count and event density",
    subtitle = sprintf("FX: %d exons/tx; BL: %d genes, %.0f%% event ratio, %d pos tx/gene",
                       cfg$fixed_n_exons_per_tx,
                       BL$n_genes, BL$event_ratio * 100, BL$n_pos_per_gene),
    x = NULL,
    y = "Median runtime (ms)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

p

# ggsave("benchmarks/results/bench_find_se.png", p, width = 8, height = 4, dpi = 150)
# message("Saved → benchmarks/results/bench_find_se.png")
