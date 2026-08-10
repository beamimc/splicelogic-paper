make_cds_gr <- function(tx_names, txps, grl) {
  exon_list <- vector("list", length(tx_names))
  for (i in seq_along(tx_names)) {
    tx_name <- tx_names[i]
    txid <- txps$tx_id[txps$transcript_name == tx_name]
    if (length(txid) == 0) next
    gr <- grl[[as.character(txid)]]
    if (is.null(gr)) next
    mcols(gr)$tx_id <- tx_name
    mcols(gr)$gene_id <- sub("-\\d+$", "", tx_name)
    exon_list[[i]] <- gr
  }
  bind_ranges(exon_list)
}

message(
  "validation/utils.R loaded:\n",
  "  make_cds_gr(tx_names, txps, cbt) -- build a CDS GRanges for a vector of transcript names"
)
