# Paper figure: nucleotide composition of the flanks of SE exons vs a non-SE
# background, one panel per nucleotide (RBP_workflow.Rmd, GTEx v9 ONT RBP
# knockdown).
# Run from the repo root. Output: figures/flank_profiles_se_vs_non_se.png

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(dplyr)
  library(ggplot2)
  library(here)
})
devtools::load_all(here("../splicelogic"))
source(here("workflows/helper.R"))

my_paste <- \(x) paste(unique(x), collapse = ",")

# 10 nt windows, each starting 5 nt after the previous one (overlap = 10 - 5),
# over 100 nt flanks: (100 - 10) / 5 = 18 steps, so 19 windows tile each flank
# exactly. With a remainder the last window comes out shorter, at the exon end
# of the flank on + strand and the far end on -
flank_width <- 100
window_width <- 10
window_overlap <- 5

# one row per exon (gene, chromosome, start, end, strand), as in the workflows
collapse_events <- function(events) {
  events |>
    dplyr::mutate(exon_key = paste(gene_id, seqnames, start, end, strand, sep = ":")) |>
    dplyr::group_by(exon_key, event_type) |>
    plyranges::reduce_ranges_directed(
      n = plyranges::n(),
      gene_id = my_paste(gene_id),
      exon_id = my_paste(exon_id),
      tx_id = my_paste(tx_id),
      exon_rank = my_paste(exon_rank),
      event_tx_id = my_paste(event_tx_id)
    )
}

# event exons vs sampled non-event exons -> the two named profile sets
make_profiles <- function(event_exons, label, dtu_res, exons) {
  event_profiles <- slEventFlankProfiles(
    event_exons, genome = "hg38", width = flank_width,
    window_width = window_width, overlap = window_overlap
  )
  non_event <- slSampleNonEventExons(
    dtu_res, exons, n = event_exons, pval_min = 0.5, seed = 1
  )
  non_event_profiles <- slEventFlankProfiles(
    non_event, genome = "hg38", width = flank_width,
    window_width = window_width, overlap = window_overlap
  )
  setNames(
    list(event_profiles, non_event_profiles),
    c(label, paste0("non-", label))
  )
}

plot_flank_figure <- function(profiles) {
  p <- slPlotFlankProfilesByNt(profiles)
  # each panel spans its own data, so drop the expand_limits(y = 0) that
  # slPlotFlankProfilesByNt() adds to start every panel at 0
  p$layers <- Filter(\(l) !inherits(l$geom, "GeomBlank"), p$layers)
  # x positions as slPlotFlankProfilesByNt() lays them out (its default
  # gap = 4): upstream windows at x = 1..nwin, the exon at exon_x, downstream
  # windows at x = nwin + 5..2 * nwin + 4
  step <- attr(profiles[[1]], "window_width") - attr(profiles[[1]], "overlap")
  xs <- p$scales$get_scales("x")
  exon_x <- xs$breaks[xs$labels == "exon"]
  nwin <- exon_x - 2.5
  up0 <- nwin            # the window next to the 3' splice site is position 0
  down0 <- nwin + 5      # the window next to the 5' splice site is position 0
  bp <- seq(30, (nwin - 1) * step, by = 30)
  # one tick per bp break, none at 0 or under the exon; "exon" is an axis
  # label too, so it lines up with the bp labels, but its tick is transparent
  breaks <- c(up0 - rev(bp) / step, exon_x, down0 + bp / step)
  labels <- c(-rev(bp), "exon", paste0("+", bp))
  tick_col <- ifelse(labels == "exon", "transparent", "black")
  # a single x axis line in the bottom panel, broken only across the exon gap:
  # from the y axis to position 0 upstream, and from position 0 downstream on
  edge <- data.frame(nt = tail(levels(factor(p$data$nt)), 1))
  seg <- function(x, xend) {
    ggplot2::geom_segment(
      data = edge, inherit.aes = FALSE, colour = "black", linewidth = 0.5,
      ggplot2::aes(x = x, xend = xend, y = -Inf, yend = -Inf)
    )
  }
  # blank the gridlines from position 0 upstream to position 0 downstream; put
  # first, so it sits under the profiles and covers only the gridlines
  p <- p + ggplot2::annotate("rect", xmin = up0, xmax = down0,
                             ymin = -Inf, ymax = Inf, fill = "white")
  p$layers <- c(utils::tail(p$layers, 1), utils::head(p$layers, -1))
  p +
    # the nucleotide to the left of each panel, as a plain label; each panel
    # with its own y range
    ggplot2::facet_wrap(~ nt, ncol = 1, strip.position = "left",
                        scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = breaks, labels = labels) +
    # a little padding at both ends, and two or three round breaks per panel
    ggplot2::scale_y_continuous(
      breaks = function(l) {
        b <- pretty(l, n = 3)
        b[b >= l[1] & b <= l[2]]
      },
      expand = ggplot2::expansion(mult = 0.05)
    ) +
    # not clipped, so the axis segments are drawn at full width like the y axis
    ggplot2::coord_cartesian(clip = "off") +
    seg(-Inf, up0) +
    seg(down0, Inf) +
    ggplot2::theme(
      axis.line.x = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_line(colour = "black", linewidth = 0.5),
      axis.ticks.y = ggplot2::element_line(colour = "black", linewidth = 0.5),
      axis.ticks.x = ggplot2::element_line(colour = tick_col, linewidth = 0.5),
      axis.ticks.length = ggplot2::unit(3.5, "pt"),
      panel.spacing.y = ggplot2::unit(0.5, "lines"),
      axis.title = ggplot2::element_text(size = 19),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 14)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 14)),
      axis.text = ggplot2::element_text(colour = "black", size = 16),
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 4)),
      strip.placement = "outside",
      # 19pt, as the axis titles and the strips of the benchmark figure
      strip.text.y.left = ggplot2::element_text(angle = 0, face = "plain",
                                                size = 19),
      legend.text = ggplot2::element_text(size = 16)
    )
}

# ---------------------------------------------------------------------------
# SE vs non-SE: RBP workflow
# ---------------------------------------------------------------------------
rbp_dtu_res_file <- here("workflows/data/dtu_res.rds")

if (!file.exists(rbp_dtu_res_file)) {
  # the RBP workflow only saves the significant transcripts; the background
  # needs every tested one, so rerun the DTU fit once and cache the table
  suppressPackageStartupMessages({
    library(satuRn)
    library(limma)
  })
  se <- readRDS(here("workflows/data/GTEx_v9_ONT_SummarizedExperiment.rds"))
  se <- se[rowSums(assay(se) >= 10) >= 6, ]
  colData(se)$condition <- dplyr::case_when(
    colData(se)$condition == "ctrl" ~ "WT",
    colData(se)$condition == "exp" ~ "KD"
  )
  colData(se)$condition <- factor(colData(se)$condition, levels = c("KD", "WT"))
  condition <- colData(se)$condition
  donor <- factor(colData(se)$donor)
  design <- model.matrix(~ 0 + condition + donor)
  L <- limma::makeContrasts(WT_vs_KD = conditionWT - conditionKD, levels = design)

  rowData(se)$isoform_id <- rowData(se)$tx_name
  rowData(se)$gene_id <- vapply(rowData(se)$gene_id, `[`, character(1), 1)
  se <- se[
    duplicated(rowData(se)$gene_id) |
      duplicated(rowData(se)$gene_id, fromLast = TRUE),
  ]
  se <- satuRn::fitDTU(se, formula = ~ 0 + condition + donor)
  se <- satuRn::testDTU(object = se, contrasts = L, diagplot1 = FALSE, diagplot2 = FALSE)

  gene_ids <- rowData(se)$gene_id
  rbp_dtu_res <- rowData(se)[["fitDTUResult_WT_vs_KD"]] |>
    as.data.frame() |>
    tibble::rownames_to_column("tx_name") |>
    tibble::as_tibble() |>
    dplyr::mutate(gene_id = gene_ids) |>
    dplyr::select(tx_name, gene_id, estimates, pval, regular_FDR) |>
    dplyr::arrange(pval)
  saveRDS(rbp_dtu_res, rbp_dtu_res_file)
}
rbp_dtu_res <- readRDS(rbp_dtu_res_file)

rbp_exons <- readRDS(here("workflows/data/glinos_exons.rds"))
sig_exons <- readRDS(here("workflows/data/sl_exons_glinos_0.1.rds"))

se_events <- sig_exons |> preprocess(coef_col = "estimates") |> find_se()
se_exons <- collapse_events(se_events)

se_profiles <- make_profiles(se_exons, "SE", rbp_dtu_res, rbp_exons)

# ---------------------------------------------------------------------------
# Save, each nucleotide panel with the y range of its own data
# ---------------------------------------------------------------------------
se_plot <- plot_flank_figure(se_profiles)

ggsave(here("figures/flank_profiles_se_vs_non_se.png"), se_plot,
       width = 8, height = 7, dpi = 300)
# the SE plot object too, for the combined figure (figures/plot_figure1.R)
saveRDS(se_plot, here("figures/flank_profiles_se_vs_non_se.rds"))
