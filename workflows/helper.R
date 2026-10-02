# ---------------------------------------------------------------------------
# Sequence-context helpers for RBP_workflow.Rmd
# Claude Code reused SPLain code and rewrote the functions here to be
# more general and reusable.
#
# splicelogic package itself only exports get_seq() for sequence work, so there
# is no package equivalent to call directly here.
#
# Needs: Biostrings, plyranges, GenomicRanges, GenomeInfoDb, BSgenome, dplyr,
# tidyr, ggplot2, ggseqlogo, patchwork. The shuffled control additionally needs
# universalmotif.
# ---------------------------------------------------------------------------


# ===========================================================================
# Genome + flanking regions
# ===========================================================================

.slGenome <- function(genome) {
    if (methods::is(genome, "BSgenome")) return(genome)
    if (is.character(genome) && length(genome) == 1L) {
        return(BSgenome::getBSgenome(genome))
    }
    stop("`genome` must be a BSgenome object or an assembly name, ",
         "e.g. \"hg38\".")
}

# TRUE for ranges that fall entirely within their sequence, so that flanking
# regions running off a chromosome end can be dropped rather than trimmed
# (a trimmed flank would give an uneven number of windows).
.slInBounds <- function(gr, genome) {
    lens <- GenomeInfoDb::seqlengths(genome)[
        as.character(GenomicRanges::seqnames(gr))
    ]
    GenomicRanges::start(gr) >= 1L &
        (is.na(lens) | GenomicRanges::end(gr) <= lens)
}

# The two flanking regions of each exon, strand-aware, optionally extended
# across the splice site into the exon itself. Exons on sequences absent from
# the genome, or whose flank runs off a chromosome end, are dropped from both
# sides at once so the two sets stay row-aligned.
.slFlankRegions <- function(events, genome, width, exon_bases = 0L) {

    stopifnot(methods::is(events, "GRanges"))
    if (length(events) == 0L) stop("`events` has no rows.")

    known <- as.character(GenomicRanges::seqnames(events)) %in%
        GenomeInfoDb::seqnames(genome)
    if (!all(known)) {
        warning("Dropping ", sum(!known),
                " exon(s) on sequences absent from `genome`.")
    }
    gr <- events[known]
    if (length(gr) == 0L) stop("No exons left on sequences present in `genome`.")
    # keeps getSeq() from warning about seqlevels that are no longer used
    GenomeInfoDb::seqlevels(gr) <- GenomeInfoDb::seqlevelsInUse(gr)

    up   <- plyranges::flank_upstream(gr, width = width)
    down <- plyranges::flank_downstream(gr, width = width)

    if (exon_bases > 0L) {
        # anchoring the intronic end means the region grows across the splice
        # site into the exon, on either strand
        up   <- plyranges::stretch(plyranges::anchor_5p(up), exon_bases)
        down <- plyranges::stretch(plyranges::anchor_3p(down), exon_bases)
    }

    keep <- .slInBounds(up, genome) & .slInBounds(down, genome)
    if (!all(keep)) {
        message("Dropped ", sum(!keep),
                " exon(s) whose flanking region runs off a chromosome end.")
    }
    if (!any(keep)) stop("No exon has both flanking regions in bounds.")

    list(upstream = up[keep], downstream = down[keep])
}


# ===========================================================================
# Nucleotide composition in sliding windows across the flanks
# ===========================================================================

# Per-sequence nucleotide fractions, one row per sequence.
.slBpPercent <- function(seqs) {
    counts <- Biostrings::oligonucleotideFrequency(seqs, width = 1)
    counts / rowSums(counts)
}

# Slide overlapping windows along each range, fetch its sequence, and return a
# matrix of nucleotide fractions: one row per input range, one column per
# window/nucleotide pair ("w1_A", "w1_C", ...). Windows are numbered in
# transcription order, so w1 is the 5'-most window on either strand.
.slWindowProfile <- function(gr, genome, window_width = 10, overlap = 5) {
    step <- window_width - overlap
    if (step <= 0L) stop("`overlap` must be smaller than `window_width`.")

    windows <- plyranges::slide_ranges(gr, width = window_width, step = step)
    seqs <- Biostrings::RNAStringSet(Biostrings::getSeq(genome, windows))

    df <- as.data.frame(.slBpPercent(seqs))
    df$partition <- windows$partition
    df$strand    <- as.character(GenomicRanges::strand(windows))

    nts <- intersect(c("A", "C", "G", "U"), names(df))
    long <- df |>
        dplyr::group_by(.data$partition) |>
        dplyr::mutate(
            window_order = if (dplyr::first(.data$strand) == "-") {
                dplyr::n() - dplyr::row_number() + 1L
            } else {
                dplyr::row_number()
            }
        ) |>
        dplyr::ungroup() |>
        tidyr::pivot_longer(cols = dplyr::all_of(nts),
                            names_to = "nt", values_to = "value") |>
        dplyr::mutate(col = paste0("w", .data$window_order, "_", .data$nt))

    wide <- long |>
        dplyr::select("partition", "col", "value") |>
        tidyr::pivot_wider(id_cols = "partition", names_from = "col",
                           values_from = "value") |>
        dplyr::arrange(.data$partition)

    mat <- as.matrix(dplyr::select(wide, -"partition"))
    rownames(mat) <- as.character(wide$partition)
    # Columns in window order, nucleotides alphabetical within each window.
    win <- as.integer(sub("^w(\\d+)_.*$", "\\1", colnames(mat)))
    mat[, order(win, colnames(mat)), drop = FALSE]
}

#' Nucleotide composition upstream and downstream of event exons
#'
#' Walks overlapping windows across the flanking region on each side of a set of
#' exons and returns the per-window nucleotide fractions. Flanks are
#' strand-aware, so "upstream" is 5' of the exon in transcription order.
#'
#' @param events A `GRanges` of exons (a `find_*()` result works directly).
#' @param genome A `BSgenome` object, or an assembly name such as `"hg38"`.
#' @param width Width (bp) of the flanking region on each side.
#' @param window_width Window width (bp).
#' @param overlap Overlap (bp) between consecutive windows.
#' @return A list with elements `upstream` and `downstream`, each a matrix with
#'   one row per retained exon and one column per window/nucleotide pair.
slEventFlankProfiles <- function(events, genome, width = 100,
                                 window_width = 10, overlap = 5) {

    genome  <- .slGenome(genome)
    regions <- .slFlankRegions(events, genome, width)

    structure(
        list(
            upstream   = .slWindowProfile(regions$upstream, genome,
                                          window_width, overlap),
            downstream = .slWindowProfile(regions$downstream, genome,
                                          window_width, overlap)
        ),
        width        = width,
        window_width = window_width,
        overlap      = overlap
    )
}

#' Plot the nucleotide composition flanking a set of exons
#'
#' Mean nucleotide fraction per window, upstream and downstream panels split by
#' a gap standing in for the exon itself, with standard-error bars.
#'
#' @param profiles Either the list returned by [slEventFlankProfiles()] or the
#'   upstream matrix, with the downstream matrix passed as `downstream`.
#' @param downstream Downstream matrix, when `profiles` is a matrix.
#' @param exon_label Label drawn in the middle box, e.g. `"SE"`.
#' @param width,step Flank width (bp) and distance (bp) between consecutive
#'   window starts. Taken from the attributes of `profiles` when available.
#' @param gap Width of the exon box, in window units.
#' @param ylim y-axis range (the panel is clipped to it, means use every window).
#' @param palette Named colour vector for the four nucleotides.
#' @return A `ggplot` object.
slPlotFlankProfiles <- function(profiles, downstream = NULL,
                                exon_label = "exon",
                                width   = NULL,
                                step    = NULL,
                                gap     = 8,
                                ylim    = c(0, 0.7),
                                palette = c(A = "#F84040", C = "skyblue",
                                            G = "#FFB400", U = "#06D6A0")) {

    if (is.list(profiles) && !is.matrix(profiles)) {
        if (is.null(width)) width <- attr(profiles, "width")
        if (is.null(step)) {
            step <- attr(profiles, "window_width") - attr(profiles, "overlap")
        }
        downstream <- profiles$downstream
        profiles   <- profiles$upstream
    }
    if (is.null(downstream)) {
        stop("Pass either the list from slEventFlankProfiles() or both ",
             "matrices.")
    }
    if (is.null(width)) width <- 100
    if (is.null(step))  step  <- 5

    to_long <- function(mat, set) {
        as.data.frame(mat) |>
            tibble::rownames_to_column("rep") |>
            tidyr::pivot_longer(
                -"rep",
                names_to  = c("window", "nt"),
                names_sep = "_",
                values_to = "value"
            ) |>
            dplyr::mutate(
                set     = set,
                win_idx = as.integer(sub("w", "", .data$window))
            )
    }

    up   <- to_long(profiles,   "upstream")
    down <- to_long(downstream, "downstream")
    nwin <- max(up$win_idx)

    up$x   <- up$win_idx
    down$x <- down$win_idx + nwin + gap
    df     <- dplyr::bind_rows(up, down)

    breaks <- c(seq_len(nwin), seq_len(nwin) + nwin + gap)
    labels <- c(-width + seq_len(nwin) * step, seq_len(nwin) * step)
    labels[seq(2, length(labels), by = 2)] <- ""

    xmin <- nwin + 0.5
    xmax <- nwin + gap + 0.5
    ymid <- mean(ylim)
    yoff <- diff(ylim) * 0.07

    ggplot2::ggplot(df, ggplot2::aes(
        x = .data$x, y = .data$value, color = .data$nt,
        group = interaction(.data$set, .data$nt)
    )) +
        ggplot2::annotate("rect",
                          xmin = xmin, xmax = xmax,
                          ymin = ymid - yoff, ymax = ymid + yoff,
                          fill = "grey70") +
        ggplot2::annotate("text",
                          x = (xmin + xmax) / 2, y = ymid,
                          label = exon_label, fontface = "bold") +
        ggplot2::stat_summary(fun      = mean,    geom = "line",
                              linewidth = 1) +
        ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar",
                              width = 0.2) +
        ggplot2::scale_color_manual(values = palette) +
        ggplot2::scale_x_continuous(breaks = breaks, labels = labels) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        # clipped rather than filtered, so the means use every window
        ggplot2::coord_cartesian(ylim = ylim) +
        ggplot2::labs(x = "Relative location (bp)",
                      y = "Nucleotide percentage (%)",
                      color = "Nucleotide") +
        ggplot2::theme_classic() +
        ggplot2::theme(
            legend.position = "top",
            axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1)
        )
}


# ===========================================================================
# Background (non-event) exon sets
# ===========================================================================

#' Random non-event exons from transcripts that failed the DTU test
#'
#' Draws a background set the same size as an event set, sampled uniformly from
#' the exons of transcripts whose DTU p-value is above `pval_min`.
#'
#' @param dtu_table The unfiltered transcript-level DTU results.
#' @param exons A `GRangesList` of exons-by-transcript, named by transcript id.
#' @param n Number of exons to draw, or a `GRanges` whose length sets it.
#' @param tx_id_col,pval_col Column names in `dtu_table`.
#' @param pval_min Sample from transcripts with a p-value strictly above this.
#' @param internal_only Exclude first and last exons.
#' @param seed Optional integer for a reproducible draw.
#' @return A `GRanges` of `n` exons.
slSampleNonEventExons <- function(dtu_table, exons, n,
                                  tx_id_col     = "tx_name",
                                  pval_col      = "pval",
                                  pval_min      = 0.5,
                                  internal_only = TRUE,
                                  seed          = NULL) {

    if (methods::is(n, "GRanges")) n <- length(n)
    n <- as.integer(n)
    stopifnot(length(n) == 1L, !is.na(n), n > 0L)

    tab <- as.data.frame(dtu_table)
    miss <- setdiff(c(tx_id_col, pval_col), names(tab))
    if (length(miss) > 0L) {
        stop("`dtu_table` is missing column(s): ",
             paste(miss, collapse = ", "))
    }

    pvals <- tab[[pval_col]]
    tab   <- tab[!is.na(pvals) & pvals > pval_min, , drop = FALSE]
    tx    <- intersect(unique(as.character(tab[[tx_id_col]])), names(exons))
    if (length(tx) == 0L) {
        stop("No transcript with ", pval_col, " > ", pval_min,
             " is present in `exons`.")
    }

    sub  <- exons[tx]
    nex  <- S4Vectors::elementNROWS(sub)
    cand <- unlist(sub, use.names = FALSE)
    cand$tx_id <- rep(tx, nex)

    if (internal_only) {
        if (is.null(cand$exon_rank)) {
            stop("`exons` has no exon_rank mcol; set internal_only = FALSE.")
        }
        cand <- cand[cand$exon_rank > 1L & cand$exon_rank < rep(nex, nex)]
        if (length(cand) == 0L) stop("No internal exon among the candidates.")
    }

    if (n > length(cand)) {
        warning("Only ", length(cand), " candidate exon(s) available; ",
                "returning all of them.")
        n <- length(cand)
    }

    if (!is.null(seed)) set.seed(seed)
    out <- cand[sort(sample.int(length(cand), n))]
    out$pval <- tab[[pval_col]][match(out$tx_id, tab[[tx_id_col]])]
    out
}

#' Matched non-event exons from the event transcripts themselves
#'
#' For each event, the internal exon of the same transcript that is shared with
#' its event partner and closest in rank to the event exon. A within-transcript
#' control, as opposed to the unrelated-transcript background from
#' [slSampleNonEventExons()].
#'
#' @param events A `GRanges` from a `find_*()`, with `tx_id`, `event_tx_id` and
#'   `exon_rank`.
#' @param exons A `GRangesList` of exons-by-transcript, named by transcript id.
#' @return A `GRanges` of at most `length(events)` exons.
slMatchedNonEventExons <- function(events, exons) {
    stopifnot(methods::is(events, "GRanges"))

    pick_one <- function(i) {
        tx_id    <- events$tx_id[i]
        ev_tx_id <- events$event_tx_id[i]
        if (!all(c(tx_id, ev_tx_id) %in% names(exons))) return(NULL)

        ex_tx <- exons[[tx_id]]
        ex_ev <- exons[[ev_tx_id]]

        # internal exons of the reference transcript
        internal <- ex_tx$exon_rank > 1 & ex_tx$exon_rank < length(ex_tx)
        ex_int   <- ex_tx[internal]

        # shared with the partner (matched on exon_id)
        shared <- ex_int[ex_int$exon_id %in% ex_ev$exon_id]
        if (length(shared) == 0L) return(NULL)

        # closest in rank to the event exon
        target <- events$exon_rank[i]
        shared[which.min(abs(shared$exon_rank - target))]
    }

    picks    <- lapply(seq_along(events), pick_one)
    keep_idx <- !vapply(picks, is.null, logical(1))

    if (!any(keep_idx)) {
        warning("No shared internal exons found for any event.")
        return(events[FALSE])
    }
    if (!all(keep_idx)) {
        message("Dropped ", sum(!keep_idx),
                " event(s) with no shared internal exon.")
    }

    new_gr  <- do.call(c, picks[keep_idx])
    base_mc <- S4Vectors::mcols(events[keep_idx])
    base_mc$exon_id   <- new_gr$exon_id
    base_mc$exon_name <- new_gr$exon_name
    base_mc$exon_rank <- new_gr$exon_rank
    S4Vectors::mcols(new_gr) <- base_mc
    new_gr
}


# ===========================================================================
# Per-position sequence logos of the splice sites flanking an exon
#
# Additionally needs ggseqlogo and patchwork.
# ===========================================================================

# Per-position nucleotide probabilities, one column per base.
.slPositionMatrix <- function(gr, genome) {
    seqs <- Biostrings::RNAStringSet(Biostrings::getSeq(genome, gr))
    cm   <- Biostrings::consensusMatrix(seqs)[c("A", "C", "G", "U"), ,
                                              drop = FALSE]
    tot  <- colSums(cm)
    if (any(tot == 0)) {
        warning(sum(tot == 0), " position(s) had no called base (all N).")
        tot[tot == 0] <- NA_real_
    }
    sweep(cm, 2, tot, "/")
}

# Position labels relative to the splice site: negative before the exon/intron
# boundary, positive after it, with no zero.
.slFlankPositions <- function(width, exon_bases, side) {
    ex <- if (exon_bases > 0L) seq_len(exon_bases) else integer(0)
    if (side == "upstream") {
        c(seq(-width, -1L), ex)          # intron ... | exon
    } else {
        c(if (exon_bases > 0L) seq(-exon_bases, -1L) else integer(0),
          seq_len(width))                # exon | ... intron
    }
}

#' Position frequency matrices for the splice sites flanking a set of exons
#'
#' Where [slEventFlankProfiles()] averages over sliding windows, this keeps
#' every base position separate, which is what a sequence logo needs. Each
#' region runs `width` bases into the flanking intron and `exon_bases` into the
#' exon, so the splice site itself sits inside the matrix.
#'
#' The region upstream of an exon holds the 3' splice site (acceptor); the
#' region downstream holds the 5' splice site (donor).
#'
#' @param events A `GRanges` of exons (a `find_*()` result works directly).
#' @param genome A `BSgenome` object, or an assembly name such as `"hg38"`.
#' @param width Bases into the intron on each side.
#' @param exon_bases Bases into the exon on each side.
#' @return A list with elements `upstream` and `downstream`, each a 4 x N matrix
#'   of probabilities with rows `A`, `C`, `G`, `U`.
slFlankLogoMatrices <- function(events, genome, width = 25, exon_bases = 3) {

    genome <- .slGenome(genome)
    regions <- .slFlankRegions(events, genome, width, exon_bases)

    out <- lapply(c(upstream = "upstream", downstream = "downstream"),
                  function(side) {
        m <- .slPositionMatrix(regions[[side]], genome)
        colnames(m) <- .slFlankPositions(width, exon_bases, side)
        m
    })

    structure(out, width = width, exon_bases = exon_bases,
              n_exons = length(regions$upstream))
}

# One logo panel: letters scaled by probability, with the exon/intron boundary
# marked and the exonic side shaded.
.slOneLogo <- function(m, method, col_scheme, title, exon_bases, side,
                       exon_fill, break_every) {

    n     <- ncol(m)
    bound <- if (side == "upstream") n - exon_bases + 0.5 else exon_bases + 0.5
    pos   <- as.integer(colnames(m))

    # always label the two bases either side of the splice site
    keep <- abs(pos) <= 2 | seq_len(n) %% break_every == 0

    p <- ggseqlogo::ggseqlogo(m, method = method, col_scheme = col_scheme)

    if (exon_bases > 0L) {
        rng <- if (side == "upstream") c(bound, n + 0.5) else c(0.5, bound)
        p <- p + ggplot2::annotate("rect", xmin = rng[1], xmax = rng[2],
                                   ymin = -Inf, ymax = Inf,
                                   fill = exon_fill, alpha = 0.35)
        p <- p + ggplot2::geom_vline(xintercept = bound, linetype = "dashed",
                                     linewidth = 0.4, colour = "grey30")
    }

    # ggseqlogo sets its own x scale; replacing it is intended, so the
    # "scale already present" message is not worth passing on
    suppressMessages(
        p +
            ggplot2::scale_x_continuous(breaks = seq_len(n)[keep],
                                        labels = pos[keep],
                                        expand = c(0.01, 0)) +
            ggplot2::labs(title = title, x = NULL,
                          y = if (method == "bits") "bits" else "probability") +
            ggplot2::theme(
                plot.title  = ggplot2::element_text(face = "bold", size = 10),
                axis.text.x = ggplot2::element_text(size = 7)
            )
    )
}

#' Sequence logos of the splice sites flanking a set of exons
#'
#' Two logo panels, one per side, with each base drawn at a height equal to its
#' frequency (`method = "prob"`) or its information content (`method = "bits"`).
#'
#' @param mats The list returned by [slFlankLogoMatrices()].
#' @param method `"prob"` or `"bits"`.
#' @param titles Panel titles.
#' @param palette Named colour vector for the four nucleotides.
#' @param exon_fill Shading colour for the exonic bases.
#' @param break_every Label every nth position.
#' @return A `patchwork` object.
slPlotFlankLogo <- function(mats,
                            method  = c("prob", "bits"),
                            titles  = c(upstream   = "Upstream: 3' splice site (acceptor)",
                                        downstream = "Downstream: 5' splice site (donor)"),
                            palette = c(A = "#F84040", C = "skyblue",
                                        G = "#FFB400", U = "#06D6A0"),
                            exon_fill   = "grey70",
                            break_every = 5) {

    method <- match.arg(method)
    exon_bases <- attr(mats, "exon_bases")
    if (is.null(exon_bases)) exon_bases <- 0L

    cs <- ggseqlogo::make_col_scheme(chars = names(palette),
                                     cols  = unname(palette))

    panels <- lapply(c("upstream", "downstream"), function(side) {
        .slOneLogo(mats[[side]], method, cs, titles[[side]], exon_bases, side,
                   exon_fill, break_every)
    })

    patchwork::wrap_plots(panels, nrow = 1) +
        patchwork::plot_annotation(
            caption = paste0("n = ", attr(mats, "n_exons"), " exons")
        )
}


# ===========================================================================
# RBP motif scanning in the flanks
# ===========================================================================

# Motifs may be given with U or T; scanning happens in DNA space because that
# is what getSeq() returns.
.slAsDna <- function(motifs) {
    stopifnot(is.character(motifs), length(motifs) > 0L)
    if (is.null(names(motifs))) names(motifs) <- motifs
    setNames(toupper(gsub("U", "T", motifs)), names(motifs))
}

# Occurrences of one IUPAC pattern per sequence. Ambiguity codes are expanded
# in the pattern but not in the subject.
.slCountMotif <- function(seqs, pattern) {
    Biostrings::vcountPattern(
        pattern, seqs, fixed = c(pattern = FALSE, subject = TRUE)
    )
}

# One row per distinct exon. A finder returns one row per (exon, partner
# transcript) pair, so an exon pairing with several partners repeats; their
# flanking sequence is identical and would otherwise be counted repeatedly.
.slDedupExons <- function(gr) {
    key <- paste(GenomicRanges::seqnames(gr), GenomicRanges::start(gr),
                 GenomicRanges::end(gr), GenomicRanges::strand(gr))
    dup <- duplicated(key)
    if (any(dup)) {
        message("Collapsed ", sum(dup), " duplicate exon range(s) of ",
                length(gr), ".")
    }
    gr[!dup]
}

# Match start positions per sequence, as a list of integer vectors.
.slMotifStarts <- function(seqs, pattern) {
    hits <- Biostrings::vmatchPattern(
        pattern, seqs, fixed = c(pattern = FALSE, subject = TRUE)
    )
    lapply(BiocGenerics::start(hits), as.integer)
}

#' Distinct exons from a finder result
#'
#' A public wrapper around the deduplication [slFlankMotifHits()] applies
#' internally, for sizing a background set against the number of exons actually
#' scanned rather than the number of event rows.
#'
#' @param gr A `GRanges`, typically a `find_*()` result.
#' @return `gr` with repeated coordinates collapsed to one row each.
slUniqueExons <- function(gr) .slDedupExons(gr)

#' Motif hits in the flanks of a set of exons
#'
#' Flanks each exon on both sides (strand-aware) then counts occurrences of each
#' motif in each flank. Works on any `GRanges`.
#'
#' @param events A `GRanges` of exons.
#' @param genome A `BSgenome` object, or an assembly name such as `"hg38"`.
#' @param motifs Character vector of motifs, IUPAC codes allowed, U or T.
#' @param width Width (bp) of the flanking region scanned on each side.
#' @param exon_bases Bases of the exon to include across the splice site.
#' @param set Label carried into the output, e.g. `"SE"` or `"background"`.
#' @param dedup Collapse exons appearing at the same coordinates more than once.
#' @return A data frame, one row per exon / side / motif.
slFlankMotifHits <- function(events, genome, motifs,
                             width      = 250,
                             exon_bases = 0L,
                             set        = "events",
                             dedup      = TRUE) {

    genome  <- .slGenome(genome)
    motifs  <- .slAsDna(motifs)
    if (dedup) events <- .slDedupExons(events)
    regions <- .slFlankRegions(events, genome, width, exon_bases)

    one_side <- function(side) {
        gr   <- regions[[side]]
        seqs <- Biostrings::getSeq(genome, gr)
        # positional ids, so that n_exons never disagrees with the number
        # of scanned regions even if exon_id is absent or repeated
        ids  <- as.character(seq_along(gr))

        do.call(rbind, lapply(names(motifs), function(nm) {
            n <- .slCountMotif(seqs, motifs[[nm]])
            df <- data.frame(
                set     = set,
                side    = side,
                motif   = nm,
                exon    = ids,
                count   = n,
                present = n > 0L,
                width      = GenomicRanges::width(gr),
                exon_bases = as.integer(exon_bases),
                stringsAsFactors = FALSE
            )
            df$starts <- .slMotifStarts(seqs, motifs[[nm]])
            df
        }))
    }

    rbind(one_side("upstream"), one_side("downstream"))
}

#' Compare motif occurrence between an event set and a background set
#'
#' Two readouts per side and motif: presence (Fisher exact test) and density
#' (hits per kb, quasi-Poisson GLM offset by region width).
#'
#' @param event_hits,bg_hits Outputs of [slFlankMotifHits()].
#' @param p_adjust Method passed to [stats::p.adjust()].
#' @return A data frame, one row per side / motif.
slMotifTest <- function(event_hits, bg_hits, p_adjust = "BH") {

    hits <- rbind(event_hits, bg_hits)
    hits$set <- factor(hits$set, levels = c(unique(bg_hits$set)[1],
                                            unique(event_hits$set)[1]))

    keys <- unique(hits[, c("side", "motif")])

    res <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
        d  <- hits[hits$side == keys$side[i] & hits$motif == keys$motif[i], ]
        ev <- d[d$set == levels(hits$set)[2], ]
        bg <- d[d$set == levels(hits$set)[1], ]

        tab <- matrix(c(sum(ev$present),  sum(!ev$present),
                        sum(bg$present),  sum(!bg$present)),
                      nrow = 2, byrow = TRUE)
        ft <- stats::fisher.test(tab)

        fit <- stats::glm(
            count ~ set + offset(log(width)),
            data = d, family = stats::quasipoisson()
        )
        co <- summary(fit)$coefficients

        data.frame(
            side          = keys$side[i],
            motif         = keys$motif[i],
            n_event       = nrow(ev),
            n_bg          = nrow(bg),
            frac_event    = mean(ev$present),
            frac_bg       = mean(bg$present),
            odds_ratio    = unname(ft$estimate),
            or_lo         = ft$conf.int[1],
            or_hi         = ft$conf.int[2],
            p_presence    = ft$p.value,
            per_kb_event  = 1000 * sum(ev$count) / sum(ev$width),
            per_kb_bg     = 1000 * sum(bg$count) / sum(bg$width),
            rate_ratio    = exp(co[2, 1]),
            p_density     = co[2, 4],
            stringsAsFactors = FALSE
        )
    }))

    res$padj_presence <- stats::p.adjust(res$p_presence, method = p_adjust)
    res$padj_density  <- stats::p.adjust(res$p_density,  method = p_adjust)
    res[order(res$p_presence), ]
}

#' Motif hits binned by distance from the splice site
#'
#' Converts the match positions from [slFlankMotifHits()] into distances from
#' the exon boundary and bins them. Upstream distances are measured from the 3'
#' splice site, downstream from the 5' splice site; both count outward into the
#' intron, so 1 is the base adjacent to the exon.
#'
#' @param hits One or more outputs of [slFlankMotifHits()], row-bound.
#' @param width Flank width used when the hits were generated.
#' @param breaks Bin edges in bp.
#' @return A data frame of hits per kb, per set / side / motif / distance bin.
slMotifDistanceProfile <- function(hits, width = 250,
                                   breaks = c(0, 25, 50, 100, width)) {

    stopifnot(all(c("starts", "exon_bases") %in% names(hits)))

    per_row <- lapply(seq_len(nrow(hits)), function(i) {
        st <- hits$starts[[i]]
        if (length(st) == 0L) return(NULL)
        eb <- hits$exon_bases[i]
        # position 1 of each sequence is the 5'-most base in transcription
        # order, on either strand. An upstream region runs intron-then-exon, a
        # downstream one exon-then-intron, so subtracting exon_bases puts
        # distance 1 at the base adjacent to the splice site in both cases.
        dist <- if (hits$side[i] == "upstream") {
            (hits$width[i] - eb) - st + 1L
        } else {
            st - eb
        }
        data.frame(
            set   = hits$set[i],
            side  = hits$side[i],
            motif = hits$motif[i],
            exon  = hits$exon[i],
            dist  = dist,
            stringsAsFactors = FALSE
        )
    })

    rows <- do.call(rbind, per_row)
    if (is.null(rows)) {
        warning("No motif hits in `hits`.")
        return(NULL)
    }

    # matches starting inside the exon, possible only when exon_bases > 0
    inside <- rows$dist < 1L
    if (any(inside)) {
        message("Dropped ", sum(inside),
                " match(es) starting inside the exon.")
        rows <- rows[!inside, , drop = FALSE]
    }

    rows$bin <- cut(rows$dist, breaks = breaks, include.lowest = TRUE)

    n_exons <- stats::aggregate(
        exon ~ set + side + motif, data = hits,
        FUN = function(x) length(unique(x))
    )
    names(n_exons)[4] <- "n_exons"

    # every set/side/motif crossed with every bin, so a bin with no hits is
    # reported as zero rather than dropping out of the profile entirely
    grid <- merge(n_exons, data.frame(
        bin = factor(levels(rows$bin), levels = levels(rows$bin)),
        bp  = diff(breaks)
    ))

    key <- c("set", "side", "motif", "bin")
    n_hits <- stats::aggregate(exon ~ set + side + motif + bin, data = rows,
                               FUN = length)
    names(n_hits)[5] <- "n_hits"
    # exons contributing at least one hit, as opposed to hits
    n_hit_ex <- stats::aggregate(exon ~ set + side + motif + bin, data = rows,
                                 FUN = function(x) length(unique(x)))
    names(n_hit_ex)[5] <- "n_exons_hit"

    out <- merge(grid, merge(n_hits, n_hit_ex, by = key), by = key, all.x = TRUE)
    out$n_hits[is.na(out$n_hits)] <- 0L
    out$n_exons_hit[is.na(out$n_exons_hit)] <- 0L

    # hits per kb of scanned sequence
    out$per_kb    <- 1000 * out$n_hits / (out$n_exons * out$bp)
    # fraction of exons contributing at least one hit
    out$frac_exon <- out$n_exons_hit / out$n_exons

    out <- out[, c("set", "side", "motif", "bin",
                   "n_exons", "n_hits", "per_kb", "frac_exon",
                   "n_exons_hit", "bp")]
    out[order(out$set, out$side, out$motif, out$bin), ]
}

#' Dinucleotide-shuffled control for a set of flanks
#'
#' Shuffling each flank while preserving its dinucleotide composition keeps base
#' composition and destroys the motif, isolating motif structure from base
#' composition. Needs the universalmotif package.
#'
#' @param events A `GRanges` of exons, as passed to [slFlankMotifHits()].
#' @param genome,motifs,width,exon_bases As in [slFlankMotifHits()].
#' @param n_shuffles Shuffled replicates per exon.
#' @param seed Optional RNG seed.
#' @param dedup Collapse duplicate exon ranges first.
#' @return A data frame in the format of [slFlankMotifHits()], `set = "shuffled"`.
slShuffledFlankHits <- function(events, genome, motifs,
                                width      = 250,
                                exon_bases = 0L,
                                n_shuffles = 1L,
                                seed       = NULL,
                                dedup      = TRUE) {

    if (!requireNamespace("universalmotif", quietly = TRUE)) {
        stop("Package 'universalmotif' is required for the shuffled control. ",
             "Install with: BiocManager::install('universalmotif')")
    }
    if (!is.null(seed)) set.seed(seed)

    genome  <- .slGenome(genome)
    motifs  <- .slAsDna(motifs)
    if (dedup) events <- .slDedupExons(events)
    regions <- .slFlankRegions(events, genome, width, exon_bases)

    one_side <- function(side) {
        gr   <- regions[[side]]
        seqs <- Biostrings::getSeq(genome, gr)
        sh   <- universalmotif::shuffle_sequences(
            rep(seqs, n_shuffles), k = 2, method = "euler"
        )
        do.call(rbind, lapply(names(motifs), function(nm) {
            n  <- .slCountMotif(sh, motifs[[nm]])
            df <- data.frame(
                set        = "shuffled",
                side       = side,
                motif      = nm,
                exon       = paste0("shuf", seq_along(sh)),
                count      = n,
                present    = n > 0L,
                width      = Biostrings::width(sh),
                exon_bases = as.integer(exon_bases),
                stringsAsFactors = FALSE
            )
            df$starts <- .slMotifStarts(sh, motifs[[nm]])
            df
        }))
    }

    rbind(one_side("upstream"), one_side("downstream"))
}
