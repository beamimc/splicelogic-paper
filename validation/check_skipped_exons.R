library(here)
devtools::load_all(here("../splicelogic"))
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(plyranges)
library(Biostrings)
library(readr)
library(dplyr)

gtf <- here("validation/biosurfer/gencode.v42.annotation.gtf.gz")
db_path <- here("validation/biosurfer/gencode.v42.TxDb.sqlite")
txps_rds <- here("validation/biosurfer/txps.rds")
cbt_rds <- here("validation/biosurfer/cbt.rds")

if (!file.exists(gtf)) {
  download.file(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz",
    destfile = gtf
  )
}

if (!file.exists(db_path)) {
  library(txdbmaker)
  txdb <- txdbmaker::makeTxDbFromGFF(gtf, format = "gtf")
  AnnotationDbi::saveDb(txdb, db_path)
}

if (file.exists(txps_rds) && file.exists(cbt_rds)) {
  txps <- readRDS(txps_rds)
  cbt <- readRDS(cbt_rds)
} else {
  library(GenomicFeatures)
  txdb <- AnnotationDbi::loadDb(db_path)
  txps <- GenomicFeatures::transcripts(txdb)
  # transcript_name (e.g. "DDX11L2-202") is not stored in TxDb; recover from GTF
  gtf_gr <- rtracklayer::import(gtf, feature.type = "transcript")
  tx_name_map <- setNames(gtf_gr$transcript_name, gtf_gr$transcript_id)
  txps$transcript_name <- tx_name_map[txps$tx_name]
  cbt <- GenomicFeatures::cdsBy(txdb, by = "tx")
  saveRDS(txps, txps_rds)
  saveRDS(cbt, cbt_rds)
}

# load the H.s. genome
bsg <- BSgenome.Hsapiens.NCBI.GRCh38

# --- 500 simple SE cases ---------------------------------------------------

pblocks <- read_tsv(
  here("validation/biosurfer/biosurfer_gencode_v42_output/pblocks.tsv"),
  show_col_types = FALSE
)

# filter to clean internal SE cases: single event, no frameshift, no split codons
se_cases <- pblocks |>
  filter(
    events == "frozenset({'e'})",
    internal,
    !compound_splicing,
    !frameshift,
    !split_codons,
    aa_loss > 0
  ) |>
  group_by(anchor, other) |>
  filter(dplyr::n() == 1) |>
  ungroup()

set.seed(5)
sample_cases <- slice_sample(se_cases, n = 100)

# build CDS GRanges for each unique anchor (down) and other (up) transcript;
# re-number exon_rank to 1:n because cdsBy ranks reflect the full transcript
# (may skip non-coding exons), which would break preprocess()'s internal flag
# and rank ± 1 neighbor key lookup
make_cds_gr <- function(tx_names, txps, cbt) {
  exon_list <- vector("list", length(tx_names))
  for (i in seq_along(tx_names)) {
    tx_name <- tx_names[i]
    txid <- txps$tx_id[txps$transcript_name == tx_name]
    if (length(txid) == 0) next
    cds <- cbt[[as.character(txid)]]
    if (is.null(cds)) next
    cds$exon_rank <- seq_along(cds)
    mcols(cds)$tx_id <- tx_name
    mcols(cds)$gene_id <- sub("-\\d+$", "", tx_name)
    exon_list[[i]] <- cds
  }
  bind_ranges(exon_list)
}

anchor_gr <- make_cds_gr(unique(sample_cases$anchor), txps, cbt)
other_gr <- make_cds_gr(unique(sample_cases$other), txps, cbt)

# other = up (positive), anchor = down (negative)
gr <- prepare_exons_by_partition(up = other_gr, down = anchor_gr) |>
  preprocess(coef_col = "estimate")

# find_se output carries tx_id (anchor) and event_tx_id (other), so we can
# join directly back to sample_cases without any case_id tracking
se_all <- find_se(gr, type = "boundary")
GenomeInfoDb::seqlevelsStyle(se_all) <- "NCBI"

# filter to phase-0 exons (cumulative CDS before the exon is divisible by 3)
anchor_txids <- txps$tx_id[match(se_all$tx_id, txps$transcript_name)]
cbt_per_se <- cbt[as.character(anchor_txids)]

preceding <- mendoapply(function(cds, str, se_start, se_end) {
  if (str == "+") cds[end(cds) < se_start] else cds[start(cds) > se_end]
}, cbt_per_se, as.list(as.character(strand(se_all))),
   as.list(start(se_all)), as.list(end(se_all)))

se_all$phase_vec <- sum(width(preceding)) %% 3L

# keep only exons where (1) cumulative CDS before the exon is divisible by 3
# (phase-0: exon starts at a codon boundary) and (2) exon width is divisible
# by 3 (no split codon at the right boundary); both required for clean
# in-frame translation comparable to biosurfer anchor_seq
se_phase0 <- se_all[!is.na(se_all$phase_vec) & se_all$phase_vec == 0L & width(se_all) %% 3L == 0L]

n_phase0 <- as_tibble(se_phase0) |>
  inner_join(sample_cases |> select(anchor, other),
             by = c("tx_id" = "anchor", "event_tx_id" = "other")) |>
  distinct(tx_id, event_tx_id) |>
  nrow()

#n_non_phase0 <- n_found - n_phase0

se_phase0 <- se_phase0 %>%
  mutate(
    dna = get_seq(., bsg),
    aastring = translate(dna, no.init.codon = TRUE),
    aa = as.character(aastring)
  )

results <- as_tibble(se_phase0) |>
  inner_join(
    sample_cases |> select(anchor, other, anchor_seq, aa_loss),
    by = c("tx_id" = "anchor", "event_tx_id" = "other")
  ) |>
  filter(width == aa_loss * 3)

n_width_mismatch <- n_phase0 - nrow(results)
n_success <- sum(results$aa == results$anchor_seq, na.rm = TRUE)

# --- Non-phase-0 SE validation (trimmed inner AA) ---------------------------
# For exons where the cumulative CDS at the exon start is not divisible by 3,
# both boundaries split a codon. Trim the leading (3 - phase) nt and enough
# trailing nt to reach the next codon boundary; the remaining inner sequence
# is fully in-frame and comparable to biosurfer's anchor_seq minus the one
# split-codon AA that biosurfer includes at the leading boundary.
# phase 1 → biosurfer anchor_seq[1] is the split AA → compare to anchor_seq[2:]
# phase 2 → biosurfer anchor_seq[end] is the split AA → compare to anchor_seq[:-1]

se_nonphase0 <- se_all[
  !is.na(se_all$phase_vec) & se_all$phase_vec != 0L & width(se_all) %% 3L == 0L
]
GenomeInfoDb::seqlevelsStyle(se_nonphase0) <- "NCBI"

se_nonphase0 <- se_nonphase0 %>%
  mutate(
    dna_full   = get_seq(., bsg),
    trim_start = (3L - phase_vec) %% 3L,
    dna_inner  = subseq(dna_full,
                        start = trim_start + 1L,
                        width = ((width - trim_start) %/% 3L) * 3L),
    aa_inner   = as.character(translate(dna_inner, no.init.codon = TRUE))
  )

results_nonphase0 <- as_tibble(se_nonphase0) |>
  inner_join(
    sample_cases |> select(anchor, other, anchor_seq, aa_loss),
    by = c("tx_id" = "anchor", "event_tx_id" = "other")
  ) |>
  filter(width == aa_loss * 3) |>
  distinct(tx_id, event_tx_id, .keep_all = TRUE) |>
  mutate(
    anchor_inner = if_else(
      phase_vec == 1L,
      substr(anchor_seq, 2L, nchar(anchor_seq)),
      substr(anchor_seq, 1L, nchar(anchor_seq) - 1L)
    )
  )

n_nonphase0_success <- sum(
  results_nonphase0$aa_inner == results_nonphase0$anchor_inner,
  na.rm = TRUE
)

# --- Summary ----------------------------------------------------------------

message(
  "Phase-0 validation\n",
  "  ", n_phase0, " / ", nrow(sample_cases), " cases: phase-0 and testable\n",
  "  ", n_width_mismatch, " / ", nrow(sample_cases), " cases: phase-0 but width != aa_loss * 3\n",
  "  ", n_success, " / ", nrow(results), " testable cases matched biosurfer anchor_seq\n",
  "\n",
  "Non-phase-0 validation (inner sequence, split-codon AAs trimmed)\n",
  "  ", nrow(results_nonphase0), " / ", nrow(sample_cases),
      " cases: non-phase-0, width divisible by 3, testable\n",
  "  ", n_nonphase0_success, " / ", nrow(results_nonphase0),
      " testable cases matched biosurfer anchor_seq (inner)"
)
