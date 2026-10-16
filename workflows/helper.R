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

# x-axis of the flank plots: bp from the exon every `break_every`, the exon
# bar, and the splice site at the inner end of each flank. Upstream window i
# sits at -width + i * step bp and is drawn at x = i, downstream window i at
# +i * step bp, drawn at x = i + nwin + gap.
.slFlankXAxis <- function(nwin, gap, width, step, break_every, exon_label) {
    exon_x  <- nwin + (gap + 1) / 2
    cand    <- seq(break_every, width, by = break_every)
    up_bp   <- rev(-cand[-cand >= -width + step &
                             -cand <= -width + nwin * step])
    down_bp <- cand[cand >= step & cand <= nwin * step]
    list(
        exon_x = exon_x,
        breaks = c((up_bp + width) / step, nwin, exon_x, nwin + gap + 1,
                   down_bp / step + nwin + gap),
        # the splice-site labels go on a second line, under the bp labels
        labels = c(up_bp, "\n3'ss", exon_label, "\n5'ss",
                   paste0("+", down_bp))
    )
}

# The exon between the two flanks: the gridlines are blanked out across the
# gap, so the exon stands on its own, and a thin bar stands in for it,
# spanning the panel whatever its y range.
.slFlankExonBar <- function(nwin, gap, exon_x) {
    list(
        ggplot2::annotate("rect",
                          xmin = nwin + 0.5, xmax = nwin + gap + 0.5,
                          ymin = -Inf, ymax = Inf, fill = "white"),
        ggplot2::annotate("rect",
                          xmin = exon_x - 0.2, xmax = exon_x + 0.2,
                          ymin = -Inf, ymax = Inf, fill = "grey75")
    )
}

# Shared look of the flank plots: no axis lines or ticks, only light
# horizontal gridlines.
.slFlankTheme <- function() {
    ggplot2::theme_classic() +
        ggplot2::theme(
            legend.position    = "top",
            axis.line          = ggplot2::element_blank(),
            axis.ticks         = ggplot2::element_blank(),
            axis.text          = ggplot2::element_text(colour = "grey40"),
            panel.grid.major.y = ggplot2::element_line(colour = "grey85",
                                                       linewidth = 0.3)
        )
}

#' Plot the nucleotide composition flanking a set of exons
#'
#' Mean nucleotide fraction per window, upstream and downstream panels split by
#' a thin bar standing in for the exon itself, with standard-error bars.
#'
#' @param profiles Either the list returned by [slEventFlankProfiles()] or the
#'   upstream matrix, with the downstream matrix passed as `downstream`.
#' @param downstream Downstream matrix, when `profiles` is a matrix.
#' @param exon_label Label under the exon bar, e.g. `"SE"`.
#' @param exon_label_size Text size of `exon_label`, in mm as for
#'   [ggplot2::geom_text()].
#' @param width,step Flank width (bp) and distance (bp) between consecutive
#'   window starts. Taken from the attributes of `profiles` when available.
#' @param gap Space between the two flanks, in window units. The exon bar sits
#'   in the middle of it.
#' @param break_every Spacing (bp) of the labelled x-axis positions.
#' @param ylim y-axis range. The means use every window, but the panel is not
#'   clipped (the exon label is drawn below it), so a mean outside `ylim` is
#'   drawn past the panel edge.
#' @param hlines y-axis breaks, each with a light gridline. `NULL` for the
#'   default breaks.
#' @param palette Named colour vector for the four nucleotides.
#' @return A `ggplot` object.
slPlotFlankProfiles <- function(profiles, downstream = NULL,
                                exon_label      = "exon",
                                exon_label_size = 5,
                                width           = NULL,
                                step            = NULL,
                                gap             = 4,
                                break_every     = 30,
                                ylim            = c(0, 0.7),
                                hlines          = seq(ylim[1], ylim[2],
                                                      by = 0.1),
                                palette         = c(A = "#3A86FF",
                                                    C = "#F0545C",
                                                    G = "#2EB67D",
                                                    U = "#F2B705")) {

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
    # the exon label is left off the axis, which has one text size for every
    # label, and drawn under the exon bar as its own, larger text
    ax     <- .slFlankXAxis(nwin, gap, width, step, break_every, "")

    ggplot2::ggplot(df, ggplot2::aes(
        x = .data$x, y = .data$value, color = .data$nt,
        group = interaction(.data$set, .data$nt)
    )) +
        ggplot2::stat_summary(fun      = mean,    geom = "line",
                              linewidth = 0.8) +
        ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar",
                              width = 0.2) +
        .slFlankExonBar(nwin, gap, ax$exon_x) +
        # placed in points below the panel, under the two lines of axis text
        # (bp, then splice sites), so it clears the splice-site labels
        ggplot2::annotation_custom(
            grid::textGrob(
                exon_label, y = grid::unit(-22, "pt"), vjust = 1,
                gp = grid::gpar(fontsize = exon_label_size * ggplot2::.pt,
                                fontface = "bold", col = "grey20")
            ),
            xmin = ax$exon_x, xmax = ax$exon_x, ymin = -Inf, ymax = -Inf
        ) +
        ggplot2::scale_color_manual(values = palette) +
        ggplot2::scale_x_continuous(breaks = ax$breaks, labels = ax$labels) +
        ggplot2::scale_y_continuous(
            breaks = if (length(hlines) > 0L) hlines else ggplot2::waiver(),
            expand = ggplot2::expansion(mult = 0.02)
        ) +
        # limited rather than filtered, so the means use every window; not
        # clipped, so the exon label can sit below the panel
        ggplot2::coord_cartesian(ylim = ylim, clip = "off") +
        ggplot2::labs(x = "Relative position in intronic sequence (nt)",
                      y = "Nucleotide fraction",
                      color = "Nucleotide") +
        .slFlankTheme() +
        # room for the exon label between the axis text and the axis title
        ggplot2::theme(axis.title.x = ggplot2::element_text(
            margin = ggplot2::margin(t = exon_label_size * ggplot2::.pt + 6)
        ))
}

#' Compare the flanking nucleotide composition of several exon sets
#'
#' One panel per nucleotide, with the mean fraction per window of each exon set
#' overlaid as a line and a standard-error ribbon. The first set (e.g. SE) is
#' drawn solid in the nucleotide's colour, the others (e.g. non-SE) in
#' `background` with their own line type. With `difference = TRUE`, plots the
#' first set minus the second instead, shaded to zero.
#'
#' The exon is drawn as a thin bar between the two flanks, so the x-axis runs
#' through the intronic sequence on either side: the 3' splice site ends the
#' upstream flank and the 5' splice site starts the downstream one.
#'
#' @param profiles A named list of [slEventFlankProfiles()] results, all built
#'   with the same `width`, `window_width` and `overlap`. The names label the
#'   sets in the legend.
#' @param difference If `TRUE`, plot the difference between the mean fractions
#'   of the two sets (first minus second). Needs exactly two sets.
#' @param exon_label x-axis label under the exon bar.
#' @param gap Space between the two flanks, in window units. The exon bar sits
#'   in the middle of it.
#' @param break_every Spacing (bp) of the labelled x-axis positions.
#' @param ylim y-axis range shared by every panel. `NULL` (the default) gives
#'   each nucleotide its own range, so small differences between sets stay
#'   visible, or with `difference = TRUE` one range symmetric around zero.
#' @param ncol Number of panel columns.
#' @param palette Named colour vector for the four nucleotides, used for the
#'   first set.
#' @param background Colour for every set after the first.
#' @param linetypes Line types for the sets, named by set or in set order.
#' @return A `ggplot` object.
slPlotFlankProfilesByNt <- function(profiles,
                                    difference  = FALSE,
                                    exon_label  = "exon",
                                    gap         = 4,
                                    break_every = 30,
                                    ylim        = NULL,
                                    ncol        = 1,
                                    palette     = c(A = "#3A86FF",
                                                    C = "#F0545C",
                                                    G = "#2EB67D",
                                                    U = "#F2B705"),
                                    background  = "grey50",
                                    linetypes   = c("solid", "dashed",
                                                    "dotted", "dotdash")) {

    if (!is.list(profiles) || is.matrix(profiles) ||
        is.null(names(profiles)) || any(names(profiles) == "")) {
        stop("`profiles` must be a named list of slEventFlankProfiles() ",
             "results, e.g. list(SE = se_profiles, `non-SE` = bg_profiles).")
    }
    params <- lapply(profiles, function(p) {
        c(attr(p, "width"), attr(p, "window_width"), attr(p, "overlap"))
    })
    if (length(unique(params)) != 1L) {
        stop("All profiles must share the same width, window_width and ",
             "overlap.")
    }
    width <- params[[1]][1]
    step  <- params[[1]][2] - params[[1]][3]

    sets <- names(profiles)
    if (difference && length(sets) != 2L) {
        stop("`difference = TRUE` needs exactly two sets, got ",
             length(sets), ".")
    }
    if (length(linetypes) < length(sets)) {
        stop("`linetypes` has ", length(linetypes), " line types for ",
             length(sets), " sets.")
    }
    if (is.null(names(linetypes))) {
        linetypes <- stats::setNames(linetypes[seq_along(sets)], sets)
    }

    to_long <- function(mat, set, side) {
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
                side    = side,
                win_idx = as.integer(sub("w", "", .data$window))
            )
    }

    up <- dplyr::bind_rows(lapply(sets, function(s) {
        to_long(profiles[[s]]$upstream, s, "upstream")
    }))
    down <- dplyr::bind_rows(lapply(sets, function(s) {
        to_long(profiles[[s]]$downstream, s, "downstream")
    }))
    nwin <- max(up$win_idx)

    up$x   <- up$win_idx
    down$x <- down$win_idx + nwin + gap
    # the first set takes the nucleotide's colour, the rest share `background`
    df     <- dplyr::bind_rows(up, down) |>
        dplyr::mutate(
            hue = ifelse(.data$set == sets[1], .data$nt, .data$set),
            set = factor(.data$set, levels = sets)
        )
    hues <- c(palette, stats::setNames(rep(background, length(sets) - 1L),
                                       sets[-1]))
    ax   <- .slFlankXAxis(nwin, gap, width, step, break_every, exon_label)

    if (difference) {
        df <- df |>
            dplyr::group_by(.data$nt, .data$side, .data$x, .data$set) |>
            dplyr::summarise(value = mean(.data$value), .groups = "drop") |>
            tidyr::pivot_wider(names_from = "set", values_from = "value") |>
            dplyr::mutate(
                value = .data[[sets[1]]] - .data[[sets[2]]]
            )
        if (is.null(ylim)) {
            ylim <- c(-1, 1) * max(abs(df$value), na.rm = TRUE) * 1.1
        }
        # gridlines at zero and one round value either side of it
        top <- pretty(c(0, ylim[2]), n = 2)
        top <- max(top[top <= ylim[2]])

        p <- ggplot2::ggplot(df, ggplot2::aes(
            x = .data$x, y = .data$value, group = .data$side
        )) +
            ggplot2::geom_hline(yintercept = 0, colour = "grey40",
                                linewidth = 0.4) +
            ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = .data$value,
                                              fill = .data$nt),
                                 alpha = 0.3) +
            ggplot2::geom_line(ggplot2::aes(colour = .data$nt),
                               linewidth = 0.8) +
            ggplot2::scale_colour_manual(values = palette, guide = "none") +
            ggplot2::scale_fill_manual(values = palette, guide = "none") +
            ggplot2::scale_y_continuous(
                breaks = c(-top, 0, top),
                labels = scales::label_number(style_positive = "plus"),
                expand = ggplot2::expansion(mult = 0.02)
            ) +
            ggplot2::labs(
                y        = paste0("\u0394 nucleotide fraction (", sets[1],
                                  " \u2212 ", sets[2], ")"),
                subtitle = paste0("Above the line: higher in ", sets[1],
                                  " \u00b7 below: lower in ", sets[1])
            )
    } else {
        p <- ggplot2::ggplot(df, ggplot2::aes(
            x = .data$x, y = .data$value,
            group = interaction(.data$side, .data$set)
        )) +
            ggplot2::stat_summary(ggplot2::aes(fill = .data$hue),
                                  fun.data = ggplot2::mean_se,
                                  geom = "ribbon", alpha = 0.25, colour = NA) +
            ggplot2::stat_summary(ggplot2::aes(colour   = .data$hue,
                                               linetype = .data$set),
                                  fun = mean, geom = "line",
                                  linewidth = 0.8) +
            # the panel titles already name the nucleotide, so only the sets
            # get a legend
            ggplot2::scale_colour_manual(values = hues, guide = "none") +
            ggplot2::scale_fill_manual(values = hues, guide = "none") +
            ggplot2::scale_linetype_manual(
                values = linetypes,
                guide  = ggplot2::guide_legend(override.aes = list(
                    colour = c("grey20", rep(background, length(sets) - 1L))
                ))
            ) +
            # little padding above and below, so the profiles fill the panel,
            # and two or three gridlines per panel at round values
            ggplot2::scale_y_continuous(
                breaks = function(l) {
                    b <- pretty(l, n = 4)
                    b <- b[b >= l[1] & b <= l[2]]
                    if (length(b) > 3L) b <- b[seq(1L, length(b), by = 2L)]
                    b
                },
                expand = ggplot2::expansion(mult = 0.02)
            ) +
            ggplot2::labs(y = "Nucleotide fraction", linetype = NULL)
    }

    p +
        .slFlankExonBar(nwin, gap, ax$exon_x) +
        ggplot2::facet_wrap(
            ~ nt, ncol = ncol, strip.position = "left",
            scales = if (is.null(ylim)) "free_y" else "fixed"
        ) +
        ggplot2::scale_x_continuous(breaks = ax$breaks, labels = ax$labels) +
        # clipped rather than filtered, so the means use every window
        ggplot2::coord_cartesian(ylim = ylim) +
        ggplot2::labs(x = "Relative position in intronic sequence (nt)") +
        .slFlankTheme() +
        ggplot2::theme(
            legend.key.width   = ggplot2::unit(1.5, "cm"),
            strip.background   = ggplot2::element_blank(),
            strip.placement    = "outside",
            strip.text.y.left  = ggplot2::element_text(angle = 0,
                                                       face = "bold",
                                                       size = 14),
            # room between the panels so each exon bar is separate
            panel.spacing.y    = ggplot2::unit(1, "lines")
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
                            palette = c(A = "#3A86FF", C = "#F0545C",
                                        G = "#2EB67D", U = "#F2B705"),
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
