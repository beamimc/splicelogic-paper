library(here)
# splicelogic devel, loaded from source: clone the repo next to this one
#   git clone -b devel https://github.com/thelovelab/splicelogic.git ../splicelogic
# then load it with devtools (install.packages("devtools") if needed)
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

source(here("validation/utils.R"))

n_sample <- 2000

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
sample_a5ss <- slice_sample(a5ss_cases, n = min(n_sample, nrow(a5ss_cases)))

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
sample_a3ss <- slice_sample(a3ss_cases, n = min(n_sample, nrow(a3ss_cases)))

# --- A5SS detection ----------------------------------------------------------
# find_a5ss output: tx_id = other (positive), event_tx_id = anchor (negative)

a5ss_anchor_gr <- make_cds_gr(unique(sample_a5ss$anchor), txps, cbt)
a5ss_other_gr  <- make_cds_gr(unique(sample_a5ss$other),  txps, cbt)

a5ss_gr <- prepare_exons_by_partition(up = a5ss_other_gr, down = a5ss_anchor_gr) |>
  preprocess(coef_col = "estimate")

a5ss_all <- find_a5ss(a5ss_gr)

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

# --- Sequence validation -----------------------------------------------------
# find_a5ss()/find_a3ss() return the other (positive) exon. Recover the
# overlapping anchor exon; the two share one boundary, and the bases between
# the two alternative splice sites are the delta. The delta lies in whichever
# exon is longer, so biosurfer reports it as anchor_seq (aa_loss) when the
# anchor exon is longer and other_seq (aa_gain) when the other exon is.
# Phase is the CDS position of the delta within that transcript. Non-phase-0
# deltas trim the split-codon AA at the leading boundary, as for SE/IE.

check_alt_ss_seq <- function(ss_all, anchor_gr, sample_cases) {
  # keep pairs detected exactly once: multi-detections are separate exons
  # that biosurfer collapses into one pblock
  single <- as_tibble(ss_all) |>
    count(tx_id, event_tx_id) |>
    filter(n == 1L)
  ss_all <- ss_all |>
    filter(paste(tx_id, event_tx_id) %in% paste(single$tx_id, single$event_tx_id))

  hits <- findOverlaps(ss_all, anchor_gr)
  hits <- hits[ss_all$event_tx_id[queryHits(hits)] == anchor_gr$tx_id[subjectHits(hits)]]
  other_ex  <- ss_all[queryHits(hits)]
  anchor_ex <- anchor_gr[subjectHits(hits)]

  delta_tbl <- tibble(
    other     = other_ex$tx_id,
    anchor    = other_ex$event_tx_id,
    seqnames  = as.character(seqnames(other_ex)),
    strand    = as.character(strand(other_ex)),
    o_start   = start(other_ex),
    o_end     = end(other_ex),
    a_start   = start(anchor_ex),
    a_end     = end(anchor_ex)
  ) |>
    inner_join(
      sample_cases |>
        select(anchor, other, anchor_seq, other_seq, aa_loss, aa_gain, split_codons),
      by = c("anchor", "other")
    ) |>
    mutate(
      one_boundary  = xor(o_start == a_start, o_end == a_end),
      anchor_longer = (a_end - a_start) > (o_end - o_start),
      d_start       = if_else(o_start == a_start, pmin(o_end, a_end) + 1L, pmin(o_start, a_start)),
      d_end         = if_else(o_start == a_start, pmax(o_end, a_end), pmax(o_start, a_start) - 1L),
      width         = d_end - d_start + 1L,
      seq_tx        = if_else(anchor_longer, anchor, other),
      exp_seq       = if_else(anchor_longer, anchor_seq, other_seq),
      exp_n_aa      = if_else(anchor_longer, aa_loss, aa_gain)
    )

  testable <- delta_tbl |>
    filter(
      one_boundary,
      width %% 3L == 0L,
      width == exp_n_aa * 3,
      # selenoproteins: biosurfer records U, Biostrings::translate() gives *
      !grepl("U", exp_seq)
    )

  delta <- GRanges(
    testable$seqnames,
    IRanges(testable$d_start, testable$d_end),
    strand = testable$strand
  )
  seq_txids <- txps$tx_id[match(testable$seq_tx, txps$transcript_name)]
  delta_tx <- GenomicFeatures::pmapToTranscripts(delta, cbt[as.character(seq_txids)])
  testable$phase <- (start(delta_tx) - 1L) %% 3L

  GenomeInfoDb::seqlevelsStyle(delta) <- "NCBI"
  delta <- delta %>%
    mutate(
      phase      = testable$phase,
      dna_full   = get_seq(., bsg),
      trim_start = (3L - phase) %% 3L,
      dna_inner  = subseq(dna_full,
                          start = trim_start + 1L,
                          width = ((width - trim_start) %/% 3L) * 3L),
      aa_inner   = as.character(translate(dna_inner, no.init.codon = TRUE))
    )

  results <- testable |>
    mutate(
      aa_inner  = delta$aa_inner,
      exp_inner = case_when(
        phase == 0L ~ exp_seq,
        phase == 1L ~ substr(exp_seq, 2L, nchar(exp_seq)),
        phase == 2L ~ substr(exp_seq, 1L, nchar(exp_seq) - 1L)
      ),
      match = aa_inner == exp_inner
    )

  list(delta_tbl = delta_tbl, results = results)
}

a5ss_seq <- check_alt_ss_seq(a5ss_all, a5ss_anchor_gr, sample_a5ss)
a3ss_seq <- check_alt_ss_seq(a3ss_all, a3ss_anchor_gr, sample_a3ss)

# --- Summary -----------------------------------------------------------------

summarize_seq <- function(label, ss_seq, n_cases) {
  d <- ss_seq$delta_tbl
  r <- ss_seq$results
  message(
    label, " sequence validation\n",
    "  ", nrow(d), " / ", n_cases, " cases: detected exactly once\n",
    "  ", sum(!d$one_boundary), " cases: both boundaries differ (excluded)\n",
    "  ", sum(d$one_boundary & d$width != d$exp_n_aa * 3), " cases: delta width != n_aa * 3 (",
        sum(d$one_boundary & d$width != d$exp_n_aa * 3 & d$split_codons), " of them split_codons)\n",
    "  ", nrow(r), " cases testable (",
        sum(r$phase == 0L), " phase-0, ", sum(r$phase != 0L), " non-phase-0; ",
        sum(r$anchor_longer), " anchor-longer, ", sum(!r$anchor_longer), " other-longer)\n",
    "  ", sum(r$match, na.rm = TRUE), " / ", nrow(r), " matched biosurfer sequence\n",
    "    phase-0:     ", sum(r$match[r$phase == 0L], na.rm = TRUE), " / ", sum(r$phase == 0L), "\n",
    "    non-phase-0: ", sum(r$match[r$phase != 0L], na.rm = TRUE), " / ", sum(r$phase != 0L)
  )
  failures <- r |> filter(!match | is.na(match))
  if (nrow(failures) > 0L) {
    message(label, " failures:")
    print(failures |> select(anchor, other, phase, width, exp_n_aa, aa_inner, exp_inner))
  }
}

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

summarize_seq("A5SS", a5ss_seq, nrow(sample_a5ss))

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

summarize_seq("A3SS", a3ss_seq, nrow(sample_a3ss))
