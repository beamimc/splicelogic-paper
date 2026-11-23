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
source("benchmarks/mock_events.R")  # EVENT_TYPES, make_dataset()

# Runtime of every find_*() function over one grid (genes / event ratio / up
# transcripts per gene; the grid of bench_find_se.R), once per
# find_*() function, each on a dataset carrying injected events of its own
# type, made by that event's generate_*() helper. One bench::mark() call per
# cell, timing each `type` the finder accepts (find_se / find_ie / find_cse /
# find_cie / find_mxe; the other six have no `type` argument).
# Iterations that hit a garbage collection are kept (filter_gc = FALSE), so its
# cost is part of the timing, and each iteration's gc level is saved.
#
# Run from the repo root. All eleven finders:
#   Rscript benchmarks/bench_finders.R
# or only some of them, leaving the other finders' files as they are:
#   Rscript benchmarks/bench_finders.R ri atss
#
# Writes, per finder, to benchmarks/results/finders/:
#   bench_find_<event>.csv            one row per type and cell (mean, sd, ...),
#                                     plotted by plot_bench_finders.R
#   bench_find_<event>_iterations.csv one row per iteration, with its gc level
#   bench_find_<event>_marks.rds      the bench_mark objects, for autoplot()
#   grid_config_find_<event>.json     the grid the three files were run on

# ---- Grid configuration (edit here) ------------------------------------
sweep_n_genes        <- c(500, 1000, 2000, 5000)
sweep_event_ratio    <- c(0.1, 0.5, 1.0)
sweep_n_pos_per_gene <- c(1:10)

baseline <- list(n_genes = 1000, event_ratio = 1.0, n_pos_per_gene = 1)

# Fixed parameters held constant across all cells
fixed_n_exons_per_tx <- 15
# timed runs of each finder (and each of its types) per cell
n_iterations <- 200
# `type` values timed for the finders that have the argument
types <- c("over", "boundary", "identical")
# finders to run, by event name: every event in EVENT_TYPES
# (benchmarks/mock_events.R)
finders <- names(EVENT_TYPES)
# -------------------------------------------------------------------------

out_dir <- "benchmarks/results/finders"

# finders named on the command line replace the `finders` set above
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) finders <- args
unknown <- setdiff(finders, names(EVENT_TYPES))
if (length(unknown) > 0) {
  stop("Unknown finder(s): ", paste(unknown, collapse = ", "),
       ". Choose from: ", paste(names(EVENT_TYPES), collapse = ", "))
}

size_grid <- bind_rows(
  tibble(n_genes = sweep_n_genes,     event_ratio = baseline$event_ratio,  n_pos_per_gene = baseline$n_pos_per_gene),
  tibble(n_genes = baseline$n_genes,  event_ratio = sweep_event_ratio,     n_pos_per_gene = baseline$n_pos_per_gene),
  tibble(n_genes = baseline$n_genes,  event_ratio = baseline$event_ratio,  n_pos_per_gene = sweep_n_pos_per_gene)
) |> distinct()

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cell_id <- function(n_genes, event_ratio, n_pos_per_gene) {
  sprintf("genes=%d_ratio=%.2f_pos=%d", n_genes, event_ratio, n_pos_per_gene)
}

# ---- One finder over the whole grid -------------------------------------

bench_finder <- function(event) {
  finder <- EVENT_TYPES[[event]]
  fn     <- as.name(finder$find)
  # the type column: one value per timed expression, NA for a finder that has
  # no `type` argument
  finder_types <- if (finder$typed) types else NA_character_

  # one call per type, for bench::mark(exprs = )
  type_exprs <- if (finder$typed) {
    setNames(lapply(types, function(type) bquote(.(fn)(gr, type = .(type)))),
             types)
  } else {
    setNames(list(bquote(.(fn)(gr))), finder$find)
  }

  jsonlite::write_json(
    list(
      finder               = finder$find,
      baseline             = baseline,
      fixed_n_exons_per_tx = fixed_n_exons_per_tx,
      n_iterations         = n_iterations,
      types                = if (finder$typed) types else list(),
      run_date             = as.character(Sys.Date()),
      splicelogic_version  = read.dcf("../splicelogic/DESCRIPTION",
                                      fields = "Version")[[1]]
    ),
    file.path(out_dir, sprintf("grid_config_find_%s.json", event)),
    auto_unbox = TRUE
  )

  message(sprintf("%s(): %d grid cells x %d expression(s) x %d iterations ...",
                  finder$find, nrow(size_grid), length(type_exprs),
                  n_iterations))

  # the same seed for every finder, so a finder's datasets do not depend on
  # which other finders were run before it
  set.seed(5)

  cells <- pmap(size_grid, function(n_genes, event_ratio, n_pos_per_gene) {
    n_events <- round(event_ratio * n_genes)
    # 1 down + n_pos up transcripts per gene
    gr <- make_dataset(event, n_genes = n_genes, n_pos = n_pos_per_gene,
                       n_neg = 1L, n_exons_per_tx = fixed_n_exons_per_tx,
                       n_events = n_events)

    # one untimed call per expression: what the finder returns on this
    # dataset, so the CSV doubles as a check that the injected events are
    # found. n_injected counts the exons the generator flagged; n_sim_hits
    # the returned exons that carry the flag.
    detected <- map_dfr(seq_along(type_exprs), function(k) {
      hits <- eval(type_exprs[[k]])
      tibble(
        type       = finder_types[[k]],
        n_exons    = length(gr),
        n_injected = sum(gr$sim_event, na.rm = TRUE),
        n_hits     = length(hits),
        n_sim_hits = if (length(hits) == 0L) 0L else
          sum(hits$sim_event, na.rm = TRUE)
      )
    })
    message(sprintf(
      "  n_genes=%d  event_ratio=%.2f  n_pos=%d  n_events=%d  exons=%d  hits=%s",
      n_genes, event_ratio, n_pos_per_gene, n_events, length(gr),
      paste(detected$n_hits, collapse = "/")
    ))

    # one untimed gc() per cell, so the garbage left from building the
    # dataset is not collected inside the first timed iterations
    gc()
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
    per_itr <- map_dfr(seq_along(type_exprs), function(k) {
      tibble(type = finder_types[[k]],
             iteration = seq_along(bm$time[[k]]),
             time_s = as.numeric(bm$time[[k]]),
             bm$gc[[k]],
             mem_alloc = as.numeric(bm$mem_alloc)[[k]])
    }) |>
      mutate(
        gc = case_when(level2 > 0 ~ "level2",
                       level1 > 0 ~ "level1",
                       level0 > 0 ~ "level0",
                       TRUE       ~ "none"),
        n_genes        = n_genes,
        event_ratio    = event_ratio,
        n_pos_per_gene = n_pos_per_gene,
        n_events       = n_events
      )

    detected <- detected |>
      mutate(n_genes = n_genes, event_ratio = event_ratio,
             n_pos_per_gene = n_pos_per_gene)

    list(per_itr = per_itr, detected = detected, bm = bm)
  })

  per_itr <- map_dfr(cells, "per_itr") |>
    mutate(finder = finder$find, .before = 1)
  detected <- map_dfr(cells, "detected")
  # the bench_mark objects too, for bench's autoplot()
  bms <- setNames(map(cells, "bm"), pmap_chr(size_grid, cell_id))

  write.csv(per_itr,
            file.path(out_dir, sprintf("bench_find_%s_iterations.csv", event)),
            row.names = FALSE)
  saveRDS(bms, file.path(out_dir, sprintf("bench_find_%s_marks.rds", event)))

  # one row per type and cell, the columns plot_bench_finders.R reads
  results <- per_itr |>
    group_by(finder, type, n_genes, event_ratio, n_pos_per_gene, n_events) |>
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
    ) |>
    left_join(detected,
              by = c("type", "n_genes", "event_ratio", "n_pos_per_gene"))

  write.csv(results, file.path(out_dir, sprintf("bench_find_%s.csv", event)),
            row.names = FALSE)
  message(sprintf(
    "  Saved -> %s/bench_find_%s.csv (+ _iterations.csv, _marks.rds); %.1f s timed",
    out_dir, event, sum(per_itr$time_s)
  ))

  results
}

# ---- Run ----------------------------------------------------------------

message("Benchmarking ", length(finders), " finder(s): ",
        paste(finders, collapse = ", "))

started <- Sys.time()
all_results <- map_dfr(finders, bench_finder)

message(sprintf(
  "Total compute (sum of all timed iterations): %.1f s; wall clock %.1f min",
  sum(all_results$mean_s * all_results$n_itr),
  as.numeric(difftime(Sys.time(), started, units = "mins"))
))

# Detection at the baseline cell: exons flagged by the generator (n_injected),
# exons returned (n_hits) and returned exons carrying the flag (n_sim_hits)
message("Baseline detection:")
all_results |>
  filter(n_genes == baseline$n_genes,
         event_ratio == baseline$event_ratio,
         n_pos_per_gene == baseline$n_pos_per_gene) |>
  select(finder, type, n_events, n_injected, n_hits, n_sim_hits) |>
  as.data.frame() |>
  print()
