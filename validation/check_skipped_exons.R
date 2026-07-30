library(here)
library(txdbmaker)
library(GenomicFeatures)
devtools::load_all(here("../splicelogic"))
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(plyranges)
library(Biostrings)
library(readr)
library(dplyr)

gtf <- here("validation/biosurfer/gencode.v42.annotation.gtf.gz")
db_path <- here("validation/biosurfer/gencode.v42.TxDb.sqlite")

# txdb <- makeTxDbFromGFF(gtf, format = "gtf")
# saveDb(txdb, file = db_path)

txdb <- loadDb(db_path)

txps <- transcripts(txdb)

# transcript_name (e.g. "DDX11L2-202") is not stored in TxDb; recover from GTF
gtf_gr <- rtracklayer::import(gtf, feature.type = "transcript")
tx_name_map <- setNames(gtf_gr$transcript_name, gtf_gr$transcript_id)
txps$transcript_name <- tx_name_map[txps$tx_name]

bsg <- BSgenome.Hsapiens.NCBI.GRCh38
ebt <- exonsBy(txdb, by = "tx")

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

  anchor_exons <- ebt[[as.character(anchor_txid)]]
  other_exons <- ebt[[as.character(other_txid)]]

  if (is.null(anchor_exons) || is.null(other_exons)) next

  anchor_exons$tx_id <- anchor_name
  anchor_exons$gene_id <- gene_id_case
  anchor_exons$estimate <- -1L
  anchor_exons$case_id <- i

  other_exons$tx_id <- other_name
  other_exons$gene_id <- gene_id_case
  other_exons$estimate <- 1L
  other_exons$case_id <- i

  exon_list[[i]] <- c(anchor_exons, other_exons)
}

gr <- do.call(c, Filter(Negate(is.null), exon_list)) |>
  preprocess(coef_col = "estimate")

se_all <- find_se(gr, type = "boundary")

# check results per case
n_success <- 0L

for (i in seq_len(nrow(sample_cases))) {
  row <- sample_cases[i, ]

  se <- se_all[se_all$case_id == i]
  if (length(se) == 0) next
  if (BiocGenerics::width(se) != row$aa_loss * 3) next

  skipped <- se
  GenomeInfoDb::seqlevelsStyle(skipped) <- "NCBI"
  dna <- get_seq(skipped, bsg)
  aa <- as.character(translate(dna))

  if (aa == row$anchor_seq) n_success <- n_success + 1L
}

message(n_success, " / ", nrow(sample_cases), " SE cases matched biosurfer anchor_seq")
