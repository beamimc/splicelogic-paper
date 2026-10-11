suppressPackageStartupMessages({
  devtools::load_all("../splicelogic") # or library(splicelogic) if installed
  library(bench)
  library(dplyr)
  library(jsonlite)
  library(purrr)
  library(tidyr)
})

# Times every find_*() finder on the same grid, so the event types are
# directly comparable. Beyond the input-size dimensions swept by
# bench_preprocess.R, the grid varies the two transcript partitions
# independently: a gene with n_pos up-regulated and n_neg down-regulated
# transcripts has n_pos * n_neg pairs to compare, and every finder works
# pair by pair within a gene.

# ---- Grid configuration (edit here) ------------------------------------
sweep_n_genes        <- c(100, 200, 500, 1000, 2000)
sweep_n_pos_per_gene <- c(1, 2, 3, 5, 10)
sweep_n_neg_per_gene <- c(1, 2, 3, 5, 10)
# balanced sweep: n_pos = n_neg = k, so pairs grow as k^2
sweep_n_pairs        <- c(1, 2, 3, 5, 10)

# The baseline gene count is what the three partition sweeps sit at, so it
# also sets how big the widest cells get: at 10 up + 10 down it is already
# 200 * 21 * 15 exons. The generate_*() helpers, not the finders, are what
# slows down on those sets, and they run once per cell as untimed setup.
baseline <- list(n_genes = 200, n_pos_per_gene = 1, n_neg_per_gene = 1)

# Fixed parameters held constant across all cells
fixed_n_exons_per_tx <- 15
# One event per up-regulated transcript, so a gene carries as many injected
# events as it has transcripts in the up partition. Generators that cannot
# reach that many are capped per event type (see event_ceiling()).
events_per_up_tx     <- 1
iterations           <- 10
# -------------------------------------------------------------------------

EVENTS <- c("se", "ie", "mxe", "ri", "a5ss", "a3ss", "atss", "ates")

FINDERS <- list(
  se   = function(gr) find_se(gr),
  ie   = function(gr) find_ie(gr),
  mxe  = function(gr) find_mxe(gr),
  ri   = function(gr) find_ri(gr),
  a5ss = function(gr) find_a5ss(gr),
  a3ss = function(gr) find_a3ss(gr),
  atss = function(gr) find_atss(gr),
  ates = function(gr) find_ates(gr)
)

# ---- Data construction --------------------------------------------------

# create_mock_data() guarantees only that tx_order 1 is negative and 2 is
# positive; from tx_order 3 on the sign is random. Force the first n_neg
# transcripts of each gene negative and all the rest positive, so a cell
# has exactly the requested partition sizes.
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

# The generators do not all draw from the same pool, so "one event per
# up-regulated transcript" is not reachable for every event type:
#   se / ie  - generate_se() samples one candidate exon per *gene*, so an
#              se set carries at most one event per gene however many
#              transcripts the up partition holds
#   ri, atss, ates - one event per up-regulated transcript, which is
#              exactly the requested rate
#   mxe      - one per (down transcript, internal exon pair)
#   a5ss, a3ss - one per internal exon of an up-regulated transcript
# Asking for more only makes the generator warn and cap, so cap here and
# record the per-event count in the CSV.
event_ceiling <- function(event, n_genes, n_pos, n_neg, n_exons_per_tx) {
  switch(event,
    se   = ,
    ie   = n_genes,
    ri   = ,
    atss = ,
    ates = n_genes * n_pos,
    mxe  = n_genes * n_neg * max(n_exons_per_tx - 3L, 0L),
    a5ss = ,
    a3ss = n_genes * n_pos * max(n_exons_per_tx - 2L, 0L)
  )
}

# Build a dataset carrying n_events injected events of one type, with
# exactly n_pos positive and n_neg negative transcripts per gene.
#
# Every generate_*() helper modifies positive-estimate transcripts, so
# there is no generator for included exons: an ie is an se seen from the
# other side. Build it as an se with the partition sizes swapped, then
# negate the estimates, which turns the transcript that lost an exon into
# the down-regulated one and leaves the exon included in the up-regulated
# partition.
make_dataset <- function(event, n_genes, n_pos, n_neg, n_exons_per_tx,
                         n_events) {
  swap <- identical(event, "ie")

  gr <- create_mock_data(
    n_genes        = n_genes,
    n_tx_per_gene  = n_pos + n_neg,
    n_exons_per_tx = n_exons_per_tx
  )
  gr <- assign_partition(gr, if (swap) n_pos else n_neg)

  gr <- switch(event,
    se   = generate_se(gr, n_events = n_events),
    ie   = generate_se(gr, n_events = n_events),
    mxe  = generate_mxe(gr, n_events = n_events),
    ri   = generate_ri(gr, n_events = n_events),
    a5ss = generate_a5ss(gr, n_events = n_events),
    a3ss = generate_a3ss(gr, n_events = n_events),
    atss = generate_atss(gr, n_events = n_events),
    ates = generate_ates(gr, n_events = n_events)
  )

  if (swap) gr <- mutate(gr, estimate = -estimate)
  gr
}

# ---- Grid ---------------------------------------------------------------

# Cells are labelled by the dimension they sweep rather than deduplicated,
# so the baseline is remeasured once per sweep and each sweep plots as a
# self-contained series.
size_grid <- bind_rows(
  tibble(
    sweep = "n_genes",
    n_genes = sweep_n_genes,
    n_pos_per_gene = baseline$n_pos_per_gene,
    n_neg_per_gene = baseline$n_neg_per_gene
  ),
  tibble(
    sweep = "n_pos_per_gene",
    n_genes = baseline$n_genes,
    n_pos_per_gene = sweep_n_pos_per_gene,
    n_neg_per_gene = baseline$n_neg_per_gene
  ),
  tibble(
    sweep = "n_neg_per_gene",
    n_genes = baseline$n_genes,
    n_pos_per_gene = baseline$n_pos_per_gene,
    n_neg_per_gene = sweep_n_neg_per_gene
  ),
  tibble(
    sweep = "n_pairs",
    n_genes = baseline$n_genes,
    n_pos_per_gene = sweep_n_pairs,
    n_neg_per_gene = sweep_n_pairs
  ),
  # Same balanced pair counts, but n_genes is scaled down by k so the total
  # number of exon rows stays at the baseline. Isolates the cost of having
  # more pairs to compare from the cost of simply having more input.
  tibble(
    sweep = "n_pairs_fixed_exons",
    n_genes = round(baseline$n_genes / sweep_n_pairs),
    n_pos_per_gene = sweep_n_pairs,
    n_neg_per_gene = sweep_n_pairs
  )
)

dir.create("benchmarks/results", recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(
  list(
    baseline             = baseline,
    fixed_n_exons_per_tx = fixed_n_exons_per_tx,
    events_per_up_tx     = events_per_up_tx,
    iterations           = iterations,
    events               = EVENTS
  ),
  "benchmarks/results/grid_config_find_events.json",
  auto_unbox = TRUE
)

message("Benchmarking ", nrow(size_grid), " grid cells x ",
        length(FINDERS) + 1L, " finders ...")

set.seed(5)

results <- pmap_dfr(size_grid, function(sweep, n_genes, n_pos_per_gene,
                                        n_neg_per_gene) {
  # One event per up-regulated transcript, as far as each generator's
  # candidate pool allows (se / ie fall short: one event per gene).
  requested <- round(events_per_up_tx * n_genes * n_pos_per_gene)
  n_events <- vapply(setNames(EVENTS, EVENTS), function(ev) {
    as.integer(min(requested, event_ceiling(
      ev, n_genes, n_pos_per_gene, n_neg_per_gene, fixed_n_exons_per_tx
    )))
  }, integer(1))
  message(sprintf(
    "  [%s] n_genes=%d  n_pos=%d  n_neg=%d  pairs=%d  events requested=%d (%s)",
    sweep, n_genes, n_pos_per_gene, n_neg_per_gene,
    n_pos_per_gene * n_neg_per_gene, requested,
    paste(sprintf("%s=%d", EVENTS, n_events[EVENTS]), collapse = " ")
  ))

  build <- function(ev) {
    make_dataset(
      event          = ev,
      n_genes        = n_genes,
      n_pos          = n_pos_per_gene,
      n_neg          = n_neg_per_gene,
      n_exons_per_tx = fixed_n_exons_per_tx,
      n_events       = n_events[[ev]]
    )
  }

  # One dataset at a time: at the wide cells a single set is already tens of
  # thousands of exons, and holding all eight live at once pushed the
  # process into swap, where it spent more time paging than running.
  # filter_gc = FALSE: the finders materialise tibble copies of the input,
  # so gc fires on most iterations and filtering would drop them all.
  timed <- imap_dfr(FINDERS, function(fn, ev) {
    gr <- build(ev)
    bm <- bench::mark(
      fn(gr),
      iterations = iterations,
      time_unit  = "s",
      filter_gc  = FALSE,
      check      = FALSE
    )
    # one untimed call to record what was actually detected, so the CSV
    # doubles as a check that the finders still fire across the grid
    hits <- fn(gr)
    row <- bm |>
      select(min, median, `itr/sec`, mem_alloc, n_itr, n_gc) |>
      mutate(
        event      = ev,
        n_events   = n_events[[ev]],
        n_exons    = length(gr),
        n_hits     = length(hits),
        n_sim_hits = if (length(hits) == 0L) 0L else
          sum(hits$sim_event, na.rm = TRUE)
      )
    message(sprintf("      %-5s %7d exons  %6.3f s  hits=%d",
                    ev, length(gr), as.numeric(row$median), length(hits)))
    rm(gr, bm, hits)
    gc(verbose = FALSE)
    row
  })

  # find_all_events() runs all eight finders in turn. Timed on the se
  # dataset, so the other seven finders return few or no hits: this is the
  # cost of the sweep itself, not of a set carrying every event type.
  gr_all <- build("se")
  bm_all <- bench::mark(
    find_all_events(gr_all, verbose = FALSE),
    iterations = iterations,
    time_unit  = "s",
    filter_gc  = FALSE,
    check      = FALSE
  )
  hits_all <- find_all_events(gr_all, verbose = FALSE)
  bm_all <- bm_all |>
    select(min, median, `itr/sec`, mem_alloc, n_itr, n_gc) |>
    mutate(
      event      = "all_events",
      n_events   = n_events[["se"]],
      n_exons    = length(gr_all),
      n_hits     = length(hits_all),
      n_sim_hits = if (length(hits_all) == 0L) 0L else
        sum(hits_all$sim_event, na.rm = TRUE)
    )
  message(sprintf("      %-5s %7d exons  %6.3f s  hits=%d",
                  "all", length(gr_all), as.numeric(bm_all$median),
                  length(hits_all)))
  rm(gr_all, hits_all)
  gc(verbose = FALSE)

  bind_rows(timed, bm_all) |>
    relocate(event) |>
    mutate(
      min            = as.numeric(min),
      median         = as.numeric(median),
      mem_alloc      = as.numeric(mem_alloc),
      sweep          = sweep,
      n_genes        = n_genes,
      n_pos_per_gene = n_pos_per_gene,
      n_neg_per_gene = n_neg_per_gene,
      n_pairs        = n_pos_per_gene * n_neg_per_gene,
      n_requested    = requested
    )
})

total_s <- sum(results$median * results$n_itr)
message(sprintf("Total compute (sum of median x n_itr): %.2f s", total_s))

write.csv(results, "benchmarks/results/bench_find_events.csv", row.names = FALSE)
message("Saved -> benchmarks/results/bench_find_events.csv")

# Detection summary at the baseline cell: every finder should recover the
# events injected for its own type, with no hits outside them.
baseline_rows <- results |>
  filter(sweep == "n_genes",
         n_genes == baseline$n_genes,
         event != "all_events")
message("Baseline detection (n_hits / n_sim_hits):")
print(baseline_rows |> select(event, n_events, n_hits, n_sim_hits))
