library(here)
devtools::load_all(here("../splicelogic"))
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(plyranges)
library(Biostrings)
library(readr)
library(dplyr)

txps_rds <- here("validation/biosurfer/txps.rds")
cbt_rds  <- here("validation/biosurfer/cbt.rds")

txps <- readRDS(txps_rds)
cbt  <- readRDS(cbt_rds)

bsg <- BSgenome.Hsapiens.NCBI.GRCh38

pblocks <- read_tsv(
  here("validation/biosurfer/biosurfer_gencode_v42_output/pblocks.tsv"),
  show_col_types = FALSE
)

# --- IE cases (E: exon included in other, absent from anchor) ---------------

ie_cases <- pblocks |>
  filter(
    events == "frozenset({'E'})",
    internal,
    !compound_splicing,
    !frameshift,
    !split_codons,
    aa_gain > 0
  ) |>
  group_by(anchor, other) |>
  filter(dplyr::n() == 1) |>
  ungroup()

set.seed(5)
sample_ie <- slice_sample(ie_cases, n = 2000)

source(here("validation/utils.R"))

# --- IE detection ------------------------------------------------------------

ie_anchor_gr <- make_cds_gr(unique(sample_ie$anchor), txps, cbt)
ie_other_gr  <- make_cds_gr(unique(sample_ie$other),  txps, cbt)

ie_gr <- prepare_exons_by_partition(up = ie_other_gr, down = ie_anchor_gr) |>
  preprocess(coef_col = "estimate")

ie_all <- find_ie(ie_gr, type = "boundary")
GenomeInfoDb::seqlevelsStyle(ie_all) <- "NCBI"

# --- IE coordinate validation ------------------------------------------------
# Mirror of SE validation: included exon is in the other (up) transcript.
# ie_all: tx_id = other (has the exon), event_tx_id = anchor (missing it)
# Phase is cumulative CDS in other before the included exon.

other_txids <- txps$tx_id[match(ie_all$tx_id, txps$transcript_name)]
cbt_per_ie <- cbt[as.character(other_txids)]

preceding_ie <- mendoapply(function(cds, str, ie_start, ie_end) {
  if (str == "+") cds[end(cds) < ie_start] else cds[start(cds) > ie_end]
}, cbt_per_ie, as.list(as.character(strand(ie_all))),
   as.list(start(ie_all)), as.list(end(ie_all)))

ie_all$phase_vec <- sum(width(preceding_ie)) %% 3L

ie_phase0 <- ie_all[
  !is.na(ie_all$phase_vec) & ie_all$phase_vec == 0L & width(ie_all) %% 3L == 0L
]

ie_phase0_pairs <- as_tibble(ie_phase0) |>
  inner_join(sample_ie |> select(anchor, other),
             by = c("tx_id" = "other", "event_tx_id" = "anchor")) |>
  distinct(tx_id, event_tx_id)

n_ie_phase0 <- nrow(ie_phase0_pairs)

ie_phase0 <- ie_phase0 %>%
  mutate(
    dna      = get_seq(., bsg),
    aastring = translate(dna, no.init.codon = TRUE),
    aa       = as.character(aastring)
  )

ie_results <- as_tibble(ie_phase0) |>
  inner_join(
    sample_ie |> select(anchor, other, other_seq, aa_gain),
    by = c("tx_id" = "other", "event_tx_id" = "anchor")
  ) |>
  filter(width == aa_gain * 3)

n_ie_width_mismatch <- n_ie_phase0 - nrow(ie_results)
n_ie_success <- sum(ie_results$aa == ie_results$other_seq, na.rm = TRUE)

# --- Non-phase-0 IE validation -----------------------------------------------

ie_nonphase0 <- ie_all[
  !is.na(ie_all$phase_vec) & ie_all$phase_vec != 0L & width(ie_all) %% 3L == 0L
]
GenomeInfoDb::seqlevelsStyle(ie_nonphase0) <- "NCBI"

ie_nonphase0 <- ie_nonphase0 %>%
  mutate(
    dna_full   = get_seq(., bsg),
    trim_start = (3L - phase_vec) %% 3L,
    dna_inner  = subseq(dna_full,
                        start = trim_start + 1L,
                        width = ((width - trim_start) %/% 3L) * 3L),
    aa_inner   = as.character(translate(dna_inner, no.init.codon = TRUE))
  )

ie_results_nonphase0 <- as_tibble(ie_nonphase0) |>
  inner_join(
    sample_ie |> select(anchor, other, other_seq, aa_gain),
    by = c("tx_id" = "other", "event_tx_id" = "anchor")
  ) |>
  # exclude selenoproteins: biosurfer records selenocysteine as U, but
  # Biostrings::translate() uses the standard genetic code and returns * for UGA
  filter(width == aa_gain * 3, !grepl("U", other_seq)) |>
  anti_join(ie_phase0_pairs, by = c("tx_id", "event_tx_id")) |>
  group_by(tx_id, event_tx_id) |>
  filter(n() == 1) |>
  ungroup() |>
  mutate(
    other_inner = if_else(
      phase_vec == 1L,
      substr(other_seq, 2L, nchar(other_seq)),
      substr(other_seq, 1L, nchar(other_seq) - 1L)
    )
  )

n_ie_nonphase0_success <- sum(
  ie_results_nonphase0$aa_inner == ie_results_nonphase0$other_inner,
  na.rm = TRUE
)

# --- Diagnostics -------------------------------------------------------------

failures_ie_nonphase0 <- ie_results_nonphase0 |>
  filter(aa_inner != other_inner | is.na(aa_inner)) |>
  select(tx_id, event_tx_id, phase_vec, width, aa_gain, aa_inner, other_inner, other_seq)

# --- Summary -----------------------------------------------------------------

message(
  "IE Phase-0 validation\n",
  "  ", n_ie_phase0, " / ", nrow(sample_ie), " cases: phase-0 and testable\n",
  "  ", n_ie_width_mismatch, " cases: phase-0 but width != aa_gain * 3\n",
  "  ", n_ie_success, " / ", nrow(ie_results), " testable cases matched biosurfer other_seq\n",
  "\n",
  "IE Non-phase-0 validation (inner sequence, split-codon AAs trimmed)\n",
  "  ", nrow(ie_results_nonphase0), " / ", nrow(sample_ie),
      " cases: non-phase-0, width divisible by 3, single inclusion, testable\n",
  "  ", n_ie_nonphase0_success, " / ", nrow(ie_results_nonphase0),
      " testable cases matched biosurfer other_seq (inner)\n",
  "\nIE Non-phase-0 failures (", nrow(failures_ie_nonphase0), "):"
)
print(failures_ie_nonphase0)
