suppressPackageStartupMessages({
  devtools::load_all("../splicelogic") # or library(splicelogic) if installed
  library(bench)
  library(dplyr)
})

# ---- Baseline (matches find_ri / find_se baseline for comparability) ---
n_genes        <- 1000
event_ratio    <- 1.0
n_pos_per_gene <- 1
n_exons_per_tx <- 15
# ------------------------------------------------------------------------

n_events <- round(event_ratio * n_genes)
message(sprintf("n_genes=%d  event_ratio=%.1f  n_events=%d  n_pos=%d",
                n_genes, event_ratio, n_events, n_pos_per_gene))

set.seed(5)

gr_base <- create_mock_data(
  n_genes        = n_genes,
  n_tx_per_gene  = n_pos_per_gene + 1L,
  n_exons_per_tx = n_exons_per_tx
)

# Force exactly one negative tx per gene (min tx_id = tx_order 1).
neg_tx_ids <- as.data.frame(gr_base) |>
  group_by(gene_id) |>
  summarise(tx_id = min(tx_id), .groups = "drop") |>
  pull(tx_id)

gr_base <- gr_base |>
  mutate(estimate = if_else(estimate < 0 & !tx_id %in% neg_tx_ids,
                            abs(estimate), estimate))

# Fork: generate each event type independently from the same base object.
# Both generators shift exon boundaries in positive transcripts and call
# preprocess() internally, so they don't interfere with each other.
gr_a5 <- generate_a5ss(gr_base, n_events = n_events)
gr_a3 <- generate_a3ss(gr_base, n_events = n_events)

bm <- bench::mark(
  a5ss = find_a5ss(gr_a5),
  a3ss = find_a3ss(gr_a3),
  iterations = 20,
  time_unit  = "s",
  filter_gc  = TRUE
)

print(bm)

# Single-iteration check: verify detected events match injected events
res_a5 <- find_a5ss(gr_a5)
message("a5ss injected (sim_event in input):");  print(table(gr_a5$sim_event))
message("a5ss detected (sim_event in result):"); print(table(res_a5$sim_event))

res_a3 <- find_a3ss(gr_a3)
message("a3ss injected (sim_event in input):");  print(table(gr_a3$sim_event))
message("a3ss detected (sim_event in result):"); print(table(res_a3$sim_event))
