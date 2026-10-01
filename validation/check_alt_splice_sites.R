library(here)
devtools::load_all(here("../splicelogic"))
library(plyranges)
library(readr)
library(dplyr)

txps_rds <- here("validation/biosurfer/txps.rds")
cbt_rds  <- here("validation/biosurfer/cbt.rds")

txps <- readRDS(txps_rds)
cbt  <- readRDS(cbt_rds)

pblocks <- read_tsv(
  here("validation/biosurfer/biosurfer_gencode_v42_output/pblocks.tsv"),
  show_col_types = FALSE
)

source(here("validation/utils.R"))

# --- A5SS cases (d: alternative donor in other) ------------------------------

a5ss_cases <- pblocks |>
  filter(
    events == "frozenset({'d'})",
    internal,
    !compound_splicing,
    !frameshift
  ) |>
  group_by(anchor, other) |>
  filter(dplyr::n() == 1) |>
  ungroup()

set.seed(5)
sample_a5ss <- slice_sample(a5ss_cases, n = 2000)

# --- A3SS cases (a: alternative acceptor in other) ---------------------------

a3ss_cases <- pblocks |>
  filter(
    events == "frozenset({'a'})",
    internal,
    !compound_splicing,
    !frameshift
  ) |>
  group_by(anchor, other) |>
  filter(dplyr::n() == 1) |>
  ungroup()

set.seed(5)
sample_a3ss <- slice_sample(a3ss_cases, n = 2000)

# --- A5SS detection ----------------------------------------------------------
# find_a5ss output: tx_id = other (positive), event_tx_id = anchor (negative)

a5ss_anchor_gr <- make_cds_gr(unique(sample_a5ss$anchor), txps, cbt)
a5ss_other_gr  <- make_cds_gr(unique(sample_a5ss$other),  txps, cbt)

a5ss_gr <- prepare_exons_by_partition(up = a5ss_other_gr, down = a5ss_anchor_gr) |>
  preprocess(coef_col = "estimate")

a5ss_all <- find_a5ss(a5ss_gr)
GenomeInfoDb::seqlevelsStyle(a5ss_all) <- "NCBI"

a5ss_counts <- as_tibble(a5ss_all) |>
  inner_join(
    sample_a5ss |> select(anchor, other),
    by = c("tx_id" = "other", "event_tx_id" = "anchor")
  ) |>
  count(anchor = event_tx_id, other = tx_id, name = "n_detections")

n_a5ss_detected   <- nrow(a5ss_counts)
n_a5ss_once       <- sum(a5ss_counts$n_detections == 1L)
n_a5ss_multi      <- sum(a5ss_counts$n_detections >  1L)
n_a5ss_undetected <- nrow(sample_a5ss) - n_a5ss_detected

# --- A3SS detection ----------------------------------------------------------
# find_a3ss output: tx_id = other (positive), event_tx_id = anchor (negative)

a3ss_anchor_gr <- make_cds_gr(unique(sample_a3ss$anchor), txps, cbt)
a3ss_other_gr  <- make_cds_gr(unique(sample_a3ss$other),  txps, cbt)

a3ss_gr <- prepare_exons_by_partition(up = a3ss_other_gr, down = a3ss_anchor_gr) |>
  preprocess(coef_col = "estimate")

a3ss_all <- find_a3ss(a3ss_gr)
GenomeInfoDb::seqlevelsStyle(a3ss_all) <- "NCBI"

a3ss_counts <- as_tibble(a3ss_all) |>
  inner_join(
    sample_a3ss |> select(anchor, other),
    by = c("tx_id" = "other", "event_tx_id" = "anchor")
  ) |>
  count(anchor = event_tx_id, other = tx_id, name = "n_detections")

n_a3ss_detected   <- nrow(a3ss_counts)
n_a3ss_once       <- sum(a3ss_counts$n_detections == 1L)
n_a3ss_multi      <- sum(a3ss_counts$n_detections >  1L)
n_a3ss_undetected <- nrow(sample_a3ss) - n_a3ss_detected

# --- Summary -----------------------------------------------------------------

message(
  "A5SS detection: ", n_a5ss_detected, " / ", nrow(sample_a5ss), " cases\n",
  "  detected exactly once: ", n_a5ss_once, "\n",
  "  detected > once:       ", n_a5ss_multi, "\n",
  "  not detected:          ", n_a5ss_undetected
)

if (n_a5ss_multi > 0L) {
  message("A5SS multiply-detected cases:")
  print(a5ss_counts |> filter(n_detections > 1L))
}

message(
  "A3SS detection: ", n_a3ss_detected, " / ", nrow(sample_a3ss), " cases\n",
  "  detected exactly once: ", n_a3ss_once, "\n",
  "  detected > once:       ", n_a3ss_multi, "\n",
  "  not detected:          ", n_a3ss_undetected
)

if (n_a3ss_multi > 0L) {
  message("A3SS multiply-detected cases:")
  print(a3ss_counts |> filter(n_detections > 1L))
}
