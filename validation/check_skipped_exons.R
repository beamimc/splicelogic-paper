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
sample_cases <- slice_sample(se_cases, n = 10)

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

  mcols(anchor_exons) <- cbind(mcols(anchor_exons), DataFrame(
    tx_id = anchor_name, gene_id = gene_id_case, estimate = -1L, case_id = i
  ))
  mcols(other_exons) <- cbind(mcols(other_exons), DataFrame(
    tx_id = other_name, gene_id = gene_id_case, estimate = 1L, case_id = i
  ))

  exon_list[[i]] <- c(anchor_exons, other_exons)
}

gr <- bind_ranges(exon_list) |>
  preprocess(coef_col = "estimate")

se_all <- find_se(gr, type = "boundary")

GenomeInfoDb::seqlevelsStyle(se_all) <- "NCBI"

# filter to phase-0 exons (cumulative CDS before the exon is divisible by 3)
phase_vec <- vapply(seq_along(se_all), function(i) {
  se <- se_all[i]
  anchor_txid <- txps$tx_id[txps$transcript_name == se$tx_id]
  if (length(anchor_txid) == 0) return(NA_integer_)
  cds_exons <- cbt[[as.character(anchor_txid)]]
  if (is.null(cds_exons)) return(NA_integer_)
  cds_df <- as.data.frame(cds_exons)
  cumcds <- if (as.character(strand(se)) == "+") {
    sum(cds_df$width[cds_df$end < start(se)])
  } else {
    sum(cds_df$width[cds_df$start > end(se)])
  }
  cumcds %% 3L
}, integer(1))

se_all$phase_vec <- phase_vec
n_found <- length(se_all)
n_not_in_cds <- nrow(sample_cases) - n_found

se_all <- se_all[!is.na(se_all$phase_vec) & se_all$phase_vec == 0L]
n_phase0 <- length(se_all)
n_non_phase0 <- n_found - n_phase0

se_all <- se_all %>%
  mutate(
    dna = get_seq(., bsg),
    aastring = translate(dna),
    aa = as.character(aastring)
  )

results <- as_tibble(se_all) |>
  left_join(
    mutate(sample_cases, case_id = row_number()) |>
      select(case_id, anchor_seq, aa_loss),
    by = "case_id"
  ) |>
  filter(width == aa_loss * 3)

n_success <- sum(results$aa == results$anchor_seq, na.rm = TRUE)
message(n_not_in_cds, " / ", nrow(sample_cases), " cases: skipped exon not in CDS")
message(n_non_phase0, " / ", nrow(sample_cases), " cases: non-phase-0, excluded")
message(n_success, " / ", n_phase0, " testable cases matched biosurfer anchor_seq")
