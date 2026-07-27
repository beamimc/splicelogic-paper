suppressPackageStartupMessages({
  devtools::load_all("../splicelogic") # or library(splicelogic) if installed
  library(bench)
  library(dplyr)
  library(jsonlite)
  library(purrr)
  library(tidyr)
})

# ---- Grid configuration (edit here) ------------------------------------
sweep_n_genes        <- c(100, 500, 1000)
sweep_event_ratio    <- c(0.1, 0.5, 1.0)
sweep_n_pos_per_gene <- c(1:5)

baseline <- list(n_genes = 1000, event_ratio = 1.0, n_pos_per_gene = 1)

# Fixed parameters held constant across all cells
fixed_n_exons_per_tx <- 15
# -------------------------------------------------------------------------

size_grid <- bind_rows(
  tibble(n_genes = sweep_n_genes,     event_ratio = baseline$event_ratio,  n_pos_per_gene = baseline$n_pos_per_gene),
  tibble(n_genes = baseline$n_genes,  event_ratio = sweep_event_ratio,     n_pos_per_gene = baseline$n_pos_per_gene),
  tibble(n_genes = baseline$n_genes,  event_ratio = baseline$event_ratio,  n_pos_per_gene = sweep_n_pos_per_gene)
) |> distinct()

dir.create("benchmarks/results", recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(
  list(
    baseline             = baseline,
    fixed_n_exons_per_tx = fixed_n_exons_per_tx
  ),
  "benchmarks/results/grid_config.json",
  auto_unbox = TRUE
)

message("Benchmarking ", nrow(size_grid), " grid cells × 3 types ...")

set.seed(5)

results <- pmap_dfr(size_grid, function(n_genes, event_ratio, n_pos_per_gene) {
  n_events <- round(event_ratio * n_genes)
  message(sprintf("  n_genes=%d  event_ratio=%.2f  n_pos=%d  n_events=%d",
                  n_genes, event_ratio, n_pos_per_gene, n_events))

  gr <- create_mock_data(
    n_genes        = n_genes,
    n_tx_per_gene  = n_pos_per_gene + 1L,  # 1 negative + n_pos positive
    n_exons_per_tx = fixed_n_exons_per_tx
  )
  gr <- generate_se(gr, n_events = n_events)

  # Force all but the most-negative tx per gene to positive coefficient.
  # create_mock_data() guarantees tx_order=1 (min tx_id per gene) is negative;
  # tx_order=3+ may be random-signed, so we correct them here after generate_se.
  neg_tx_ids <- as.data.frame(gr) |>
    group_by(gene_id) |>
    summarise(tx_id = min(tx_id), .groups = "drop") |>
    pull(tx_id)

  gr <- gr |>
    mutate(estimate = if_else(estimate < 0 & !tx_id %in% neg_tx_ids,
                              abs(estimate), estimate))

  bm <- bench::mark(
    boundary = find_se(gr, type = "boundary"),
    over     = find_se(gr, type = "over"),
    `in`     = find_se(gr, type = "in"),
    iterations = 15,
    time_unit = "s",
    check = FALSE
  )

  bm |>
    select(expression, min, median, `itr/sec`, mem_alloc, n_itr, n_gc) |>
    mutate(
      min            = as.numeric(min),
      median         = as.numeric(median),
      mem_alloc      = as.numeric(mem_alloc),
      type           = as.character(expression),
      n_genes        = n_genes,
      event_ratio    = event_ratio,
      n_pos_per_gene = n_pos_per_gene,
      n_events       = n_events,
      .keep = "unused"
    )
})

total_s <- sum(results$median * results$n_itr)
message(sprintf("Total compute (sum of median × n_itr): %.2f s", total_s))

write.csv(results, "benchmarks/results/bench_find_se.csv", row.names = FALSE)
message("Saved → benchmarks/results/bench_find_se.csv")
