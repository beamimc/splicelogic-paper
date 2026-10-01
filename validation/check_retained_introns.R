library(here)
devtools::load_all(here("../splicelogic"))
library(plyranges)
library(readr)
library(dplyr)

txps_rds <- here("validation/biosurfer/txps.rds")
ebt_rds  <- here("validation/biosurfer/ebt.rds")

txps <- readRDS(txps_rds)
ebt  <- readRDS(ebt_rds)

pblocks <- read_tsv(
  here("validation/biosurfer/biosurfer_gencode_v42_output/pblocks.tsv"),
  show_col_types = FALSE
)

# --- RI cases (I: anchor splices out intron, other retains it) --------------
# other=up (has the big retaining exon), anchor=down (has the split exons)

ri_cases <- pblocks |>
  filter(
    events == "frozenset({'I'})",
    !compound_splicing,
    !frameshift,
    !split_codons,
    aa_gain > 0
  ) |>
  group_by(anchor, other) |>
  filter(dplyr::n() == 1) |>
  ungroup()

set.seed(5)
sample_ri <- slice_sample(ri_cases, n = 2000)

source(here("validation/utils.R"))

# --- RI detection ------------------------------------------------------------
# use exon-level (not CDS) data: the retaining exon in other spans the intron
# fully at the exon boundary, whereas cdsBy can truncate at an early stop codon

ri_anchor_gr <- make_cds_gr(unique(sample_ri$anchor), txps, ebt)
ri_other_gr  <- make_cds_gr(unique(sample_ri$other),  txps, ebt)

ri_gr <- prepare_exons_by_partition(up = ri_other_gr, down = ri_anchor_gr) |>
  preprocess(coef_col = "estimate")

ri_all <- find_ri(ri_gr)
GenomeInfoDb::seqlevelsStyle(ri_all) <- "NCBI"

# --- RI summary --------------------------------------------------------------
# Count detections per case before deduplication to catch multiply-detected cases.
# ri_all: tx_id = other (pos, retaining), event_tx_id = anchor (neg, splicing)

ri_counts <- as_tibble(ri_all) |>
  inner_join(
    sample_ri |> select(anchor, other),
    by = c("tx_id" = "other", "event_tx_id" = "anchor")
  ) |>
  count(anchor = event_tx_id, other = tx_id, name = "n_detections")

n_detected   <- nrow(ri_counts)
n_once       <- sum(ri_counts$n_detections == 1L)
n_multi      <- sum(ri_counts$n_detections >  1L)
n_undetected <- nrow(sample_ri) - n_detected

message(
  "RI detection: ", n_detected, " / ", nrow(sample_ri), " cases\n",
  "  detected exactly once: ", n_once, "\n",
  "  detected > once:       ", n_multi, "\n",
  "  not detected:          ", n_undetected
)

if (n_multi > 0L) {
  message("Multiply-detected cases:")
  print(ri_counts |> filter(n_detections > 1L))
}
