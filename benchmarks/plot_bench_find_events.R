library(dplyr)
library(ggplot2)
library(jsonlite)

results <- read.csv("benchmarks/results/bench_find_events.csv")
cfg <- jsonlite::read_json("benchmarks/results/grid_config_find_events.json")
BL  <- cfg$baseline

# find_all_events() is the sum of the eight finders, so it sits an order of
# magnitude above them and would flatten the panels. Plotted separately.
finders <- results |>
  filter(event != "all_events") |>
  mutate(
    event     = factor(event, levels = unlist(cfg$events)),
    median_ms = median * 1000
  )

panel_labels <- c(
  n_genes             = "Genes",
  n_pos_per_gene      = "Up-regulated tx / gene",
  n_neg_per_gene      = "Down-regulated tx / gene",
  n_pairs             = "Pairs / gene (balanced)",
  n_pairs_fixed_exons = "Pairs / gene (exons held fixed)"
)

# Each sweep varies one dimension; x is whatever that sweep moved.
focal_value <- function(df) {
  df |>
    mutate(
      focal_value = case_when(
        sweep == "n_genes"             ~ as.numeric(n_genes),
        sweep == "n_pos_per_gene"      ~ as.numeric(n_pos_per_gene),
        sweep == "n_neg_per_gene"      ~ as.numeric(n_neg_per_gene),
        sweep == "n_pairs"             ~ as.numeric(n_pairs),
        sweep == "n_pairs_fixed_exons" ~ as.numeric(n_pairs)
      ),
      sweep = factor(sweep, levels = names(panel_labels))
    )
}

plot_data <- focal_value(finders)

p <- ggplot(plot_data, aes(x = focal_value, y = median_ms, color = event)) +
  geom_line() +
  geom_point(size = 1.6) +
  facet_wrap(
    ~ sweep,
    scales = "free_x",
    labeller = labeller(sweep = panel_labels)
  ) +
  scale_color_brewer(palette = "Set2", name = "Event") +
  labs(
    title    = "find_*() runtime by input size and partition size",
    subtitle = sprintf(
      "BL: %d genes, %d up / %d down tx per gene, %d exons/tx, %g event per up tx",
      BL$n_genes, BL$n_pos_per_gene, BL$n_neg_per_gene,
      cfg$fixed_n_exons_per_tx, cfg$events_per_up_tx
    ),
    x = NULL,
    y = "Median runtime (ms)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

p

ggsave("benchmarks/results/bench_find_events.png", p, width = 10, height = 6, dpi = 150)
message("Saved -> benchmarks/results/bench_find_events.png")

# ---- Pairs, with input size held constant -------------------------------
# n_pos * n_neg pairs are compared per gene. In the "n_pairs" sweep the
# exon count grows with the pair count, so the two effects are confounded;
# "n_pairs_fixed_exons" scales n_genes down to keep the exon count at the
# baseline, leaving the pair count as the only thing that moved.
pairs_data <- plot_data |>
  filter(sweep %in% c("n_pairs", "n_pairs_fixed_exons")) |>
  mutate(
    input = if_else(sweep == "n_pairs", "exons grow with pairs",
                    "exons held fixed")
  )

p_pairs <- ggplot(
  pairs_data,
  aes(x = n_pairs, y = median_ms, color = event, linetype = input)
) +
  geom_line() +
  geom_point(size = 1.6) +
  scale_x_continuous(breaks = sort(unique(pairs_data$n_pairs))) +
  scale_color_brewer(palette = "Set2", name = "Event") +
  scale_linetype_discrete(name = "Input size") +
  labs(
    title    = "find_*() runtime by number of transcript pairs per gene",
    subtitle = sprintf("n_pos = n_neg = k, so k^2 pairs per gene; BL %d exons/tx",
                       cfg$fixed_n_exons_per_tx),
    x = "Pairs compared per gene",
    y = "Median runtime (ms)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

p_pairs

ggsave("benchmarks/results/bench_find_events_pairs.png", p_pairs,
       width = 8, height = 5, dpi = 150)
message("Saved -> benchmarks/results/bench_find_events_pairs.png")

# ---- find_all_events() --------------------------------------------------
all_data <- results |>
  filter(event == "all_events") |>
  mutate(median_ms = median * 1000) |>
  focal_value()

p_all <- ggplot(all_data, aes(x = focal_value, y = median_ms)) +
  geom_line() +
  geom_point(size = 1.6) +
  facet_wrap(
    ~ sweep,
    scales = "free_x",
    labeller = labeller(sweep = panel_labels)
  ) +
  labs(
    title    = "find_all_events() runtime",
    subtitle = "All eight finders in one pass, timed on the skipped-exon dataset",
    x = NULL,
    y = "Median runtime (ms)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

p_all

ggsave("benchmarks/results/bench_find_all_events.png", p_all,
       width = 10, height = 6, dpi = 150)
message("Saved -> benchmarks/results/bench_find_all_events.png")
