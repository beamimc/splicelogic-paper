library(here)
devtools::load_all(here("../splicelogic"))
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(plyranges)
library(Biostrings)
library(readr)
library(dplyr)

txps_rds <- here("validation/biosurfer/txps.rds")
cbt_rds <- here("validation/biosurfer/cbt.rds")

txps <- readRDS(txps_rds)
cbt <- readRDS(cbt_rds)

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
sample_ie <- slice_sample(ie_cases, n = 100)

# --- RI cases (B: intron retained in other, spliced out in anchor) ----------

ri_cases <- pblocks |>
  filter(
    events == "frozenset({'B'})",
    !compound_splicing,
    !frameshift,
    !split_codons,
    aa_gain > 0
  ) |>
  group_by(anchor, other) |>
  filter(dplyr::n() == 1) |>
  ungroup()

set.seed(5)
sample_ri <- slice_sample(ri_cases, n = 100)

source(here("validation/utils.R"))

# --- IE detection ------------------------------------------------------------

ie_anchor_gr <- make_cds_gr(unique(sample_ie$anchor), txps, cbt)
ie_other_gr <- make_cds_gr(unique(sample_ie$other), txps, cbt)

ie_gr <- prepare_exons_by_partition(up = ie_other_gr, down = ie_anchor_gr) |>
  preprocess(coef_col = "estimate")

ie_all <- find_ie(ie_gr, type = "boundary")
GenomeInfoDb::seqlevelsStyle(ie_all) <- "NCBI"

# --- RI detection ------------------------------------------------------------

ri_anchor_gr <- make_cds_gr(unique(sample_ri$anchor), txps, cbt)
ri_other_gr <- make_cds_gr(unique(sample_ri$other), txps, cbt)

ri_gr <- prepare_exons_by_partition(up = ri_other_gr, down = ri_anchor_gr) |>
  preprocess(coef_col = "estimate")

ri_all <- find_ri(ri_gr)
GenomeInfoDb::seqlevelsStyle(ri_all) <- "NCBI"
