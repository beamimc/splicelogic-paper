# Mock datasets with injected events, shared by the benchmark scripts.
# Source it from the repo root, after loading splicelogic:
#   source("benchmarks/mock_events.R")

# Every event type: the generate_*() helper that injects it, the find_*()
# function that detects it, and whether that finder takes a `type` argument.
# Order as in find_all_events() and event_palette().
EVENT_TYPES <- list(
  se   = list(generate = "generate_se",   find = "find_se",   typed = TRUE),
  ie   = list(generate = "generate_ie",   find = "find_ie",   typed = TRUE),
  cse  = list(generate = "generate_cse",  find = "find_cse",  typed = TRUE),
  cie  = list(generate = "generate_cie",  find = "find_cie",  typed = TRUE),
  mxe  = list(generate = "generate_mxe",  find = "find_mxe",  typed = TRUE),
  ri   = list(generate = "generate_ri",   find = "find_ri",   typed = FALSE),
  si   = list(generate = "generate_si",   find = "find_si",   typed = FALSE),
  a5ss = list(generate = "generate_a5ss", find = "find_a5ss", typed = FALSE),
  a3ss = list(generate = "generate_a3ss", find = "find_a3ss", typed = FALSE),
  atss = list(generate = "generate_atss", find = "find_atss", typed = FALSE),
  ates = list(generate = "generate_ates", find = "find_ates", typed = FALSE)
)

# create_mock_data() guarantees only that tx_order 1 (min tx_id per gene) is
# negative and tx_order 2 positive; from the third on the sign is random.
# Force the first n_neg transcripts of each gene negative and the rest
# positive, so a dataset has exactly the partition sizes asked for.
assign_partition <- function(gr, n_neg) {
  tx <- as.data.frame(gr) |>
    distinct(gene_id, tx_id) |>
    arrange(gene_id, tx_id) |>
    group_by(gene_id) |>
    mutate(tx_order = row_number()) |>
    ungroup()
  neg_ids <- tx$tx_id[tx$tx_order <= n_neg]
  gr |>
    mutate(estimate = if_else(tx_id %in% neg_ids, -abs(estimate), abs(estimate)))
}

# n_neg down + n_pos up transcripts per gene, carrying n_events injected
# events of one type. The partition is set before the generator runs: each
# generator edits the transcripts of one sign and flags `sim_event` from the
# signs as they are when it runs, so flipping a sign afterwards would turn an
# injected event into a different one, or leave its flag on the wrong side.
make_dataset <- function(event, n_genes, n_pos, n_neg = 1L, n_exons_per_tx,
                         n_events) {
  gr <- create_mock_data(
    n_genes        = n_genes,
    n_tx_per_gene  = n_pos + n_neg,
    n_exons_per_tx = n_exons_per_tx
  )
  gr <- assign_partition(gr, n_neg)
  generate <- match.fun(EVENT_TYPES[[event]]$generate)
  generate(gr, n_events = n_events)
}
