suppressPackageStartupMessages({
  devtools::load_all("../splicelogic") # or library(splicelogic) if installed
  library(bench)
  library(dplyr)
  library(GenomicRanges)
  library(jsonlite)
  library(purrr)
  library(tidyr)
})

# ---- Grid configuration (edit here) ------------------------------------
# preprocess() is the required gate stage: it adds key/nexons/internal and
# renames coef_col to `estimate`. Its cost has no event component, so the
# sweep dimensions are the ones that set the total number of exon rows.
sweep_n_genes        <- c(100, 500, 1000, 5000)
sweep_n_tx_per_gene  <- c(2, 3, 4, 5)
sweep_n_exons_per_tx <- c(5, 10, 15, 20)

baseline <- list(n_genes = 1000, n_tx_per_gene = 2, n_exons_per_tx = 15)

# Name of the user-supplied coefficient column in the raw input. Not "estimate",
# so the benchmark exercises preprocess()'s rename path.
coef_col <- "log2FC"
# -------------------------------------------------------------------------

size_grid <- bind_rows(
  tibble(n_genes = sweep_n_genes,    n_tx_per_gene = baseline$n_tx_per_gene, n_exons_per_tx = baseline$n_exons_per_tx),
  tibble(n_genes = baseline$n_genes, n_tx_per_gene = sweep_n_tx_per_gene,    n_exons_per_tx = baseline$n_exons_per_tx),
  tibble(n_genes = baseline$n_genes, n_tx_per_gene = baseline$n_tx_per_gene, n_exons_per_tx = sweep_n_exons_per_tx)
) |> distinct()

dir.create("benchmarks/results", recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(
  list(
    baseline = baseline,
    coef_col = coef_col
  ),
  "benchmarks/results/grid_config_preprocess.json",
  auto_unbox = TRUE
)

message("Benchmarking preprocess() over ", nrow(size_grid), " grid cells ...")

set.seed(5)

results <- pmap_dfr(size_grid, function(n_genes, n_tx_per_gene, n_exons_per_tx) {
  n_exons <- n_genes * n_tx_per_gene * n_exons_per_tx
  message(sprintf("  n_genes=%d  n_tx=%d  n_exons_per_tx=%d  n_exons=%d",
                  n_genes, n_tx_per_gene, n_exons_per_tx, n_exons))

  gr <- create_mock_data(
    n_genes        = n_genes,
    n_tx_per_gene  = n_tx_per_gene,
    n_exons_per_tx = n_exons_per_tx
  )

  # create_mock_data() calls preprocess() internally, so strip everything it
  # added (key/nexons/internal + the preprocessed flag) to recover a raw,
  # user-supplied-style input. Rename estimate -> coef_col to match the
  # documented entry point, where the coefficient column is named by the user.
  mcols(gr) <- mcols(gr)[, c("gene_id", "tx_id", "exon_rank", "estimate")]
  names(mcols(gr))[names(mcols(gr)) == "estimate"] <- coef_col
  metadata(gr) <- list()

  # filter_gc = FALSE: preprocess() materialises a full tibble copy of the
  # input, so gc fires on most iterations and filtering would drop them all.
  bm <- bench::mark(
    preprocess = preprocess(gr, coef_col = coef_col),
    iterations = 10,
    time_unit  = "s",
    filter_gc  = FALSE
  )

  bm |>
    select(expression, min, median, `itr/sec`, mem_alloc, n_itr, n_gc) |>
    mutate(
      min            = as.numeric(min),
      median         = as.numeric(median),
      mem_alloc      = as.numeric(mem_alloc),
      type           = as.character(expression),
      n_genes        = n_genes,
      n_tx_per_gene  = n_tx_per_gene,
      n_exons_per_tx = n_exons_per_tx,
      n_exons        = n_exons,
      .keep = "unused"
    )
})

total_s <- sum(results$median * results$n_itr)
message(sprintf("Total compute (sum of median × n_itr): %.2f s", total_s))

write.csv(results, "benchmarks/results/bench_preprocess.csv", row.names = FALSE)
message("Saved → benchmarks/results/bench_preprocess.csv")


