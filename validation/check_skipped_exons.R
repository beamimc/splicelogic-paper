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
ebt_rds <- here("validation/biosurfer/ebt.rds")
cbt_rds <- here("validation/biosurfer/cbt.rds")

if (file.exists(txps_rds) && file.exists(ebt_rds) && file.exists(cbt_rds)) {
  txps <- readRDS(txps_rds)
  ebt <- readRDS(ebt_rds)
  cbt <- readRDS(cbt_rds)
} else {
  library(txdbmaker)
  library(GenomicFeatures)
  txdb <- AnnotationDbi::loadDb(db_path)
  txps <- GenomicFeatures::transcripts(txdb)
  # transcript_name (e.g. "DDX11L2-202") is not stored in TxDb; recover from GTF
  gtf_gr <- rtracklayer::import(gtf, feature.type = "transcript")
  tx_name_map <- setNames(gtf_gr$transcript_name, gtf_gr$transcript_id)
  txps$transcript_name <- tx_name_map[txps$tx_name]
  ebt <- GenomicFeatures::exonsBy(txdb, by = "tx")
  cbt <- GenomicFeatures::cdsBy(txdb, by = "tx")
  saveRDS(txps, txps_rds)
  saveRDS(ebt, ebt_rds)
  saveRDS(cbt, cbt_rds)
}

# load the H.s. genome
bsg <- BSgenome.Hsapiens.NCBI.GRCh38

# --- 100 simple SE cases ---------------------------------------------------

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
sample_cases <- slice_sample(se_cases, n = 500)

# build combined GRanges across all cases; gene_id is unique per case to
# prevent find_se from pairing exons across different cases
exon_list <- vector("list", nrow(sample_cases))

for (i in seq_len(nrow(sample_cases))) {
  row <- sample_cases[i, ]
  anchor_name <- row$anchor
  other_name <- row$other
  gene_id_case <- paste0(sub("-\\d+$", "", anchor_name), "_case", i)

  anchor_txid <- txps$tx_id[txps$transcript_name == anchor_name]
  other_txid <- txps$tx_id[txps$transcript_name == other_name]

  if (length(anchor_txid) == 0 || length(other_txid) == 0) next

  anchor_exons <- cbt[[as.character(anchor_txid)]]
  other_exons <- cbt[[as.character(other_txid)]]

  if (is.null(anchor_exons) || is.null(other_exons)) next

  # cdsBy exon_rank reflects transcript-level rank (may skip non-coding exons);
  # re-number to 1:n so find_se's rank ± 1 neighbor lookup doesn't produce NAs
  anchor_exons$exon_rank <- seq_along(anchor_exons)
  other_exons$exon_rank <- seq_along(other_exons)

  # tx_id must be unique per case: if the same transcript appears in two cases,
  # preprocess() groups all copies under the same tx_id and computes nexons as a
  # multiple of the true count, causing internal exon mislabeling and NA rank lookups
  mcols(anchor_exons) <- cbind(mcols(anchor_exons), DataFrame(
    tx_id = paste0(anchor_name, "_case", i), gene_id = gene_id_case,
    estimate = -1L, case_id = i
  ))
  mcols(other_exons) <- cbind(mcols(other_exons), DataFrame(
    tx_id = paste0(other_name, "_case", i), gene_id = gene_id_case,
    estimate = 1L, case_id = i
  ))

  exon_list[[i]] <- c(anchor_exons, other_exons)
}

gr <- bind_ranges(exon_list) |>
  preprocess(coef_col = "estimate")

# note: this will drop when the skipped exon lives in UTR
se_all <- find_se(gr, type = "boundary")

GenomeInfoDb::seqlevelsStyle(se_all) <- "NCBI"

# filter to phase-0 exons (cumulative CDS before the exon is divisible by 3)
# strip the per-case suffix to recover the original transcript name for cbt lookup
anchor_txids <- txps$tx_id[match(sub("_case\\d+$", "", se_all$tx_id), txps$transcript_name)]
cbt_per_se <- cbt[as.character(anchor_txids)]

# complex code that gets the preceding CDS exons and returns a GRangesList
preceding <- mendoapply(function(cds, str, se_start, se_end) {
  if (str == "+") cds[end(cds) < se_start] else cds[start(cds) > se_end]
}, cbt_per_se, as.list(as.character(strand(se_all))),
   as.list(start(se_all)), as.list(end(se_all)))

se_all$phase_vec <- sum(width(preceding)) %% 3L
n_found <- length(unique(se_all$case_id))
n_not_in_cds <- nrow(sample_cases) - n_found

# keep only exons where (1) cumulative CDS before the exon is divisible by 3
# (phase-0: exon starts at a codon boundary) and (2) exon width is divisible
# by 3 (no split codon at the right boundary); both are required for a clean
# in-frame translation that can be compared directly to biosurfer anchor_seq
se_phase0 <- se_all[!is.na(se_all$phase_vec) & se_all$phase_vec == 0L & width(se_all) %% 3L == 0L]
n_phase0 <- length(unique(se_phase0$case_id))
n_non_phase0 <- n_found - n_phase0

se_phase0 <- se_phase0 %>%
  mutate(
    dna = get_seq(., bsg),
    aastring = translate(dna, no.init.codon = TRUE),
    aa = as.character(aastring)
  )

results <- as_tibble(se_phase0) |>
  left_join(
    mutate(sample_cases, case_id = row_number()) |>
      select(case_id, anchor_seq, aa_loss),
    by = "case_id"
  ) |>
  filter(width == aa_loss * 3)

n_width_mismatch <- n_phase0 - nrow(results)
n_success <- sum(results$aa == results$anchor_seq, na.rm = TRUE)
message(n_not_in_cds, " / ", nrow(sample_cases), " cases: skipped exon not in CDS")
message(n_non_phase0, " / ", nrow(sample_cases), " cases: non-phase-0, excluded")
message(n_width_mismatch, " / ", nrow(sample_cases), " cases: phase-0 but width != aa_loss * 3")
message(n_success, " / ", nrow(results), " testable cases matched biosurfer anchor_seq")
