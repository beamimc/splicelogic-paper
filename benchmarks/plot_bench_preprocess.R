library(dplyr)
library(ggplot2)
library(jsonlite)

# ---- Runs to overlay (label = file) ------------------------------------
# Both files share the same grid, so they overlay panel-for-panel.
runs <- list(
  `check_exon_rank` = "benchmarks/results/bench_preprocess.csv",
  `no check_exon_rank`         = "benchmarks/results/bench_preprocess_no_checkexonrank.csv"
)
# -------------------------------------------------------------------------

results <- bind_rows(lapply(runs, read.csv), .id = "run") |>
  mutate(run = factor(run, levels = names(runs)))

cfg <- jsonlite::read_json("benchmarks/results/grid_config_preprocess.json")
BL  <- cfg$baseline

panel_labels <- c(
  n_genes        = "Genes",
  n_tx_per_gene  = "Tx / gene",
  n_exons_per_tx = "Exons / tx"
)

plot_data <- bind_rows(
  results |>
    filter(n_tx_per_gene == BL$n_tx_per_gene, n_exons_per_tx == BL$n_exons_per_tx) |>
    mutate(focal_dim = "n_genes", focal_value = n_genes),
  results |>
    filter(n_genes == BL$n_genes, n_exons_per_tx == BL$n_exons_per_tx) |>
    mutate(focal_dim = "n_tx_per_gene", focal_value = n_tx_per_gene),
  results |>
    filter(n_genes == BL$n_genes, n_tx_per_gene == BL$n_tx_per_gene) |>
    mutate(focal_dim = "n_exons_per_tx", focal_value = n_exons_per_tx)
) |>
  mutate(
    focal_dim = factor(focal_dim, levels = names(panel_labels)),
    median_ms = median * 1000
  )

p <- ggplot(plot_data, aes(x = focal_value, y = median_ms, color = run)) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(
    ~ focal_dim,
    scales = "free_x",
    labeller = labeller(focal_dim = panel_labels)
  ) +
  scale_color_brewer(palette = "Set2", name = "Run") +
  labs(
    title    = "preprocess() runtime by input size",
    subtitle = sprintf("BL: %d genes, %d tx/gene, %d exons/tx",
                       BL$n_genes, BL$n_tx_per_gene, BL$n_exons_per_tx),
    x = NULL,
    y = "Median runtime (ms)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

p

ggsave("benchmarks/results/bench_preprocess.png", p, width = 8, height = 4, dpi = 150)
message("Saved → benchmarks/results/bench_preprocess.png")



