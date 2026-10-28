suppressPackageStartupMessages({
  # splicelogic devel, loaded from source: clone the repo next to this one
  #   git clone -b devel https://github.com/thelovelab/splicelogic.git ../splicelogic
  # then load it with devtools (install.packages("devtools") if needed)
  devtools::load_all("../splicelogic")
  library(bench)
  library(dplyr)
  library(jsonlite)
  library(purrr)
  library(tidyr)
})

# find_se() runtime over a grid of input sizes, one bench::mark() call per
# cell timing each of `types`. Iterations that hit a garbage collection are
# kept (filter_gc = FALSE), so its cost is part of the timing, and each
# iteration's gc level is saved. Writes:
#   bench_find_se.csv            one row per type and cell (mean, sd, ...),
#                                plotted by plot_bench_find_se.R / figure 1
#   bench_find_se_iterations.csv one row per iteration, with its gc level
#   bench_find_se_marks.rds      the bench_mark objects, for bench's autoplot()

# ---- Grid configuration (edit here) ------------------------------------
sweep_n_genes        <- c(500, 1000, 2000, 5000)
sweep_event_ratio    <- c(0.1, 0.5, 1.0)
sweep_n_pos_per_gene <- c(1:10)

baseline <- list(n_genes = 1000, event_ratio = 1.0, n_pos_per_gene = 1)

# Fixed parameters held constant across all cells
fixed_n_exons_per_tx <- 15
# timed runs of each find_se() type per cell
n_iterations <- 200
# find_se() types to time: any of "over" (the default), "boundary",
# "identical"; figure 1 shows "over" only
types <- "over"
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
    fixed_n_exons_per_tx = fixed_n_exons_per_tx,
    n_iterations         = n_iterations,
    types                = types,
    run_date             = as.character(Sys.Date()),
    splicelogic_version  = read.dcf("../splicelogic/DESCRIPTION",
                                    fields = "Version")[[1]]
  ),
  "benchmarks/results/grid_config.json",
  auto_unbox = TRUE
)

# 1 negative + n_pos positive transcripts per gene, carrying n_events injected
# events
make_dataset <- function(n_genes, n_pos_per_gene, n_events) {
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

  gr |>
    mutate(estimate = if_else(estimate < 0 & !tx_id %in% neg_tx_ids,
                              abs(estimate), estimate))
}

cell_id <- function(n_genes, event_ratio, n_pos_per_gene) {
  sprintf("genes=%d_ratio=%.2f_pos=%d", n_genes, event_ratio, n_pos_per_gene)
}

# one find_se() call per type, named by type, for bench::mark(exprs = )
type_exprs <- setNames(
  lapply(types, function(type) bquote(find_se(gr, type = .(type)))),
  types
)

message("Benchmarking ", nrow(size_grid), " grid cells × ", length(types),
        " type(s) ...")

set.seed(5)

cells <- pmap(size_grid, function(n_genes, event_ratio, n_pos_per_gene) {
  n_events <- round(event_ratio * n_genes)
  message(sprintf("  n_genes=%d  event_ratio=%.2f  n_pos=%d  n_events=%d",
                  n_genes, event_ratio, n_pos_per_gene, n_events))

  gr <- make_dataset(n_genes, n_pos_per_gene, n_events)

  bm <- bench::mark(
    exprs = type_exprs,
    iterations = n_iterations,
    # keep the iterations a garbage collection landed in, so its cost shows
    filter_gc  = FALSE,
    time_unit  = "s",
    check      = FALSE
  )

  # one row per iteration, with the highest gc level that ran during it
  # (labelled as in bench's own autoplot)
  per_itr <- map_dfr(seq_along(types), function(k) {
    tibble(type = types[[k]],
           iteration = seq_along(bm$time[[k]]),
           time_s = as.numeric(bm$time[[k]]),
           bm$gc[[k]])
  }) |>
    mutate(
      gc = case_when(level2 > 0 ~ "level2",
                     level1 > 0 ~ "level1",
                     level0 > 0 ~ "level0",
                     TRUE       ~ "none"),
      mem_alloc      = as.numeric(bm$mem_alloc)[match(type, types)],
      n_genes        = n_genes,
      event_ratio    = event_ratio,
      n_pos_per_gene = n_pos_per_gene,
      n_events       = n_events
    )

  list(per_itr = per_itr, bm = bm)
})

per_itr <- map_dfr(cells, "per_itr")
# the bench_mark objects too, for bench's autoplot()
bms <- setNames(map(cells, "bm"), pmap_chr(size_grid, cell_id))

total_s <- sum(per_itr$time_s)
message(sprintf("Total compute (sum of all iterations): %.2f s", total_s))

write.csv(per_itr, "benchmarks/results/bench_find_se_iterations.csv",
          row.names = FALSE)
saveRDS(bms, "benchmarks/results/bench_find_se_marks.rds")
message("Saved → benchmarks/results/bench_find_se_iterations.csv (+ _marks.rds)")

# one row per type and cell, the columns plot_bench_find_se.R reads
results <- per_itr |>
  group_by(type, n_genes, event_ratio, n_pos_per_gene, n_events) |>
  summarise(
    min       = min(time_s),
    median    = median(time_s),
    # per-iteration spread (seconds), for error bars in the plot
    mean_s    = mean(time_s),
    sd_s      = stats::sd(time_s),
    `itr/sec` = dplyr::n() / sum(time_s),
    mem_alloc = max(mem_alloc, na.rm = TRUE),
    n_itr     = dplyr::n(),
    # collections of any level, over all iterations
    n_gc      = sum(level0 + level1 + level2),
    .groups   = "drop"
  )

write.csv(results, "benchmarks/results/bench_find_se.csv", row.names = FALSE)
message("Saved → benchmarks/results/bench_find_se.csv")

# Single-iteration check at baseline: verify detected events match injected events
gr_check <- make_dataset(baseline$n_genes, baseline$n_pos_per_gene,
                         round(baseline$event_ratio * baseline$n_genes))
message("Injected (sim_event in input):"); print(table(gr_check$sim_event))
for (type in types) {
  res_check <- find_se(gr_check, type = type)
  message("Detected, type = \"", type, "\" (sim_event in result):")
  print(table(res_check$sim_event))
}
