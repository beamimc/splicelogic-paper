suppressPackageStartupMessages({
  devtools::load_all("../splicelogic") # or library(splicelogic) if installed
  library(bench)
  library(dplyr)
})

# ---- Baseline (matches find_se baseline for comparability) -------------
n_genes        <- 1000
event_ratio    <- 1.0
n_pos_per_gene <- 1
n_exons_per_tx <- 15
# ------------------------------------------------------------------------

n_events <- round(event_ratio * n_genes)
message(sprintf("n_genes=%d  event_ratio=%.1f  n_events=%d  n_pos=%d",
                n_genes, event_ratio, n_events, n_pos_per_gene))

set.seed(5)

gr <- create_mock_data(
  n_genes        = n_genes,
  n_tx_per_gene  = n_pos_per_gene + 1L,  # 1 negative + n_pos positive
  n_exons_per_tx = n_exons_per_tx
)

# Force exactly one negative tx per gene (min tx_id = tx_order 1, always negative).
# generate_ri() modifies positive transcripts, so correct signs before calling it.
neg_tx_ids <- as.data.frame(gr) |>
  group_by(gene_id) |>
  summarise(tx_id = min(tx_id), .groups = "drop") |>
  pull(tx_id)

gr <- gr |>
  mutate(estimate = if_else(estimate < 0 & !tx_id %in% neg_tx_ids,
                            abs(estimate), estimate))

# generate_ri() merges exons 2 and 3 of positive transcripts into one large
# exon (the retained intron). find_ri() then detects introns from negative
# transcripts that fall within those merged positive exons.
gr <- generate_ri(gr, n_events = n_events)

bm <- bench::mark(
  find_ri(gr),
  iterations = 20,
  time_unit  = "s", 
  filter_gc = TRUE
)

print(bm)
