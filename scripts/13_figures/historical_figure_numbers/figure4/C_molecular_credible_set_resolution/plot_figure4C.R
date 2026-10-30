#!/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/Rscript


Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           NUMEXPR_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
ROOT <- "/path/to/project/analysis_20260817"
OUT_ROOT <- file.path(ROOT, "article/figure/figure4/C_molecular_credible_set_resolution")
INPUT_ROOT <- file.path(ROOT, "5-dogs/10.finemaping/formal/8.5_no_gwas_build_expression_ase_regulatory_catalogue/20260827_201055")
VALID <- file.path(INPUT_ROOT, "no_gwas_expression_ase_valid_block_catalogue.tsv")
INVALID <- file.path(INPUT_ROOT, "no_gwas_expression_ase_invalid_block_qc.tsv")


WIDTH_MM <- 80
HEIGHT_MM <- 54
FONT <- "Arial"
ORIGINS <- c("dog5_manual_new_gene", "reference_annotation_gene")
TISSUES <- c("cerebellum", "cortex_frontal", "cortex_occipital", "pituitary_gland")
BIN_LABELS <- c("1-10", "11-50", "51-100", "101-500", ">500")

BIN_COLOURS <- setNames(c("#08519C", "#3182BD", "#6BAED6", "#BDD7E7", "#E8EEF3"), BIN_LABELS)

main <- function() {
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out <- file.path(OUT_ROOT, stamp)
  if (!dir.create(out, recursive = TRUE, showWarnings = FALSE)) stop("output directory already exists or cannot create: ", out)
  log_dir <- file.path(OUT_ROOT, "log")
  dir.create(log_dir, showWarnings = FALSE)
  log_path <- file.path(log_dir, paste0("plot_figure4C.", stamp, ".", Sys.getpid(), ".log"))
  con <- file(log_path, "wt")
  sink(con, type = "output")
  sink(con, type = "message")
  on.exit({sink(type = "message"); sink(type = "output"); close(con)}, add = TRUE)
  log_stage <- function(s) {cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S ]"), s, "\n"); flush(con)}
  tryCatch({
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(ggplot2))
    setDTthreads(1)
    font_paths <- systemfonts::match_fonts(c(FONT, FONT), weight = c("normal", "bold"))$path
    expected_fonts <- c("/path/to/fonts/Arial/Arial.TTF",
                        "/path/to/fonts/Arial/Arialbd.TTF")
    if (!all(file.exists(expected_fonts)) ||
        !identical(normalizePath(font_paths), normalizePath(expected_fonts)))
      stop("Arial regular or bold font is unavailable; please install the specified font, do not silently substitute.")
    write_tsv <- function(x, name) fwrite(x, file.path(out, name), sep = "\t", na = "NA")
    log_stage("[1/4] read complete valid and invalid block directory; concurrent_tasks=1 threads_per_task=1")
    cat("VALID=", VALID, "\nINVALID=", INVALID, "\nOUTPUT=", out, "\n", sep = "")
    d <- fread(file = VALID, sep = "\t", na.strings = "NA", showProgress = FALSE)
    excluded <- fread(file = INVALID, sep = "\t", na.strings = "NA", showProgress = FALSE)
    required <- c("block_id", "feature_origin", "tissue", "expression_gene_id", "result_interpretable",
                  "credible_set_size", "observed_marker_count", "selected_sample_count", "max_ppa", "top_ppa_tie_count")
    if (!all(required %in% names(d))) stop("the valid catalog is missing required fields")
    if (anyDuplicated(d$block_id)) stop("block_id is not unique, cannot duplicate count")
    if (anyNA(d[, ..required]) || any(d$result_interpretable != 1L) ||
        any(!d$feature_origin %in% ORIGINS) || any(!d$tissue %in% TISSUES) ||
        any(!is.finite(d$credible_set_size) | !is.finite(d$observed_marker_count)) ||
        any(d$credible_set_size < 1 | d$credible_set_size != floor(d$credible_set_size) |
            d$credible_set_size > d$observed_marker_count)) stop("classification or credible-set fields in the valid catalog do not satisfy the plotting input contract")
    if (!all(c("block_id", "feature_origin", "tissue") %in% names(excluded)) ||
        any(!excluded$feature_origin %in% ORIGINS) || any(!excluded$tissue %in% TISSUES) ||
        anyDuplicated(excluded$block_id) || any(excluded$block_id %in% d$block_id)) stop("exclude directory identity or group mismatch")

    log_stage("[2/4] summarize mutually exclusive bins and cumulative thresholds by tissue and source; no additional association filtering")
    d[, cs_bin := as.character(cut(credible_set_size, breaks = c(0, 10, 50, 100, 500, Inf),
                                    labels = BIN_LABELS, right = TRUE))]
    d[, cs_fraction := credible_set_size / observed_marker_count]
    groups <- CJ(feature_origin = ORIGINS, tissue = TISSUES, sorted = FALSE)

    g <- d[, .(n_blocks = .N, n_expression_genes = uniqueN(expression_gene_id),
               min_dogs = min(selected_sample_count), max_dogs = max(selected_sample_count),
               median_cs_size = as.numeric(median(credible_set_size)),
               median_input_markers = as.numeric(median(observed_marker_count)),
               median_cs_fraction = as.numeric(median(cs_fraction)), tied_top_blocks = sum(top_ppa_tie_count > 1)),
           by = .(feature_origin, tissue)]
    g <- merge(groups, g, by = c("feature_origin", "tissue"), all.x = TRUE)
    if (anyNA(g$n_blocks)) stop("predefined group no valid block, cannot plot percentage records")
    bad <- excluded[, .(n_excluded_blocks = .N), by = .(feature_origin, tissue)]
    g <- merge(g, bad, by = c("feature_origin", "tissue"), all.x = TRUE)
    g[is.na(n_excluded_blocks), n_excluded_blocks := 0L]

    g[, y := 10L - match(tissue, TISSUES) - 5L * (match(feature_origin, ORIGINS) - 1L)]
    counts <- d[, .(n_bin = .N), by = .(feature_origin, tissue, cs_bin)]
    bins <- merge(CJ(feature_origin = ORIGINS, tissue = TISSUES, cs_bin = BIN_LABELS, sorted = FALSE),
                  counts, by = c("feature_origin", "tissue", "cs_bin"), all.x = TRUE)
    bins[is.na(n_bin), n_bin := 0L]
    bins <- merge(bins, g, by = c("feature_origin", "tissue"))
    bins[, bin_order := match(cs_bin, BIN_LABELS)]
    setorder(bins, feature_origin, tissue, bin_order)
    bins[, fraction := n_bin / n_blocks]
    bins[, xmax := cumsum(fraction), by = .(feature_origin, tissue)]
    bins[, xmin := xmax - fraction]
    bins[, cumulative_blocks := cumsum(n_bin), by = .(feature_origin, tissue)]
    bins[, upper_bound := c(10, 50, 100, 500, Inf)[bin_order]]
    source_fields <- unique(c(required, "cs_bin", "cs_fraction", "credible_set_ppa_sum", "input_dir", "result_dir"))
    write_tsv(d[, ..source_fields], "Figure4C_block_source.tsv")
    write_tsv(g, "Figure4C_group_summary.tsv")
    write_tsv(bins, "Figure4C_size_bins.tsv")
    write_tsv(bins[, .(feature_origin, tissue, upper_bound, cumulative_blocks, n_blocks,
                       cumulative_fraction = xmax)], "Figure4C_cumulative_thresholds.tsv")
    write_tsv(excluded[, .(block_id, feature_origin, tissue, qc_category, qc_reason)], "Figure4C_excluded_blocks.tsv")

    log_stage("[3/4] export 80 × 54 mm of shared coordinates bar figure and traceable source table")
    headings <- data.table(x = 0, y = c(9.85, 4.85),
                           label = c("Candidate loci", "Reference genes"))
    p <- ggplot() +
      geom_rect(data = bins, aes(xmin = xmin, xmax = xmax, ymin = y - 0.36, ymax = y + 0.36, fill = cs_bin),
                colour = NA, show.legend = TRUE) +
      geom_text(data = g, aes(x = 1.025, y = y, label = format(n_blocks, big.mark = ",", trim = TRUE)),
                hjust = 0, size = 6.5, size.unit = "pt", family = FONT) +
      geom_text(data = headings, aes(x = x, y = y, label = label),
                hjust = 0, fontface = "bold", size = 8, size.unit = "pt", family = FONT) +
      annotate("text", x = 1.025, y = 9.85, label = "Blocks", hjust = 0, size = 6.5, size.unit = "pt", family = FONT) +
      scale_fill_manual(values = BIN_COLOURS, breaks = BIN_LABELS, limits = BIN_LABELS, drop = FALSE,
                        name = "Variants in the 95% credible set", guide = guide_legend(nrow = 1)) +
      scale_x_continuous(breaks = seq(0, 1, 0.25), labels = c("0", "25", "50", "75", "100"),
                         expand = c(0, 0)) +
      scale_y_continuous(breaks = c(1:4, 6:9), labels = rep(c("Pit", "OC", "FC", "CB"), 2),
                         expand = c(0, 0)) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0.45, 10.25), clip = "off") +

      labs(x = "Interpretable blocks (%)", y = NULL, tag = "C") +
      theme_classic(base_size = 7, base_family = FONT) +
      theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
            axis.line.x = element_line(linewidth = 0.3), axis.ticks.x = element_line(linewidth = 0.3),
            axis.text = element_text(colour = "black", size = 7),
            axis.title.x = element_text(size = 7.5, margin = margin(t = 2)),
            legend.position = "bottom", legend.justification = "left", legend.title.position = "top",
            legend.title = element_text(size = 7), legend.text = element_text(size = 6.5),
            legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(1, "mm"),
            legend.key.spacing.x = grid::unit(0.8, "mm"),
            legend.key.width = grid::unit(2.5, "mm"), legend.key.height = grid::unit(2.5, "mm"),
            plot.tag = element_text(face = "bold", size = 10), plot.tag.position = "topleft",

            plot.margin = margin(5, 35, 5, 5))
    prefix <- file.path(out, "Figure4C_molecular_credible_set_resolution")
    ggsave(paste0(prefix, ".pdf"), p, device = grDevices::cairo_pdf, width = WIDTH_MM, height = HEIGHT_MM, units = "mm", bg = "white")
    ggsave(paste0(prefix, ".svg"), p, device = svglite::svglite, width = WIDTH_MM, height = HEIGHT_MM, units = "mm", bg = "white")
    ggsave(paste0(prefix, ".png"), p, device = "png", type = "cairo", width = WIDTH_MM, height = HEIGHT_MM, units = "mm", dpi = 600, bg = "white")
    ggsave(paste0(prefix, ".tiff"), p, device = "tiff", type = "cairo", compression = "lzw", width = WIDTH_MM, height = HEIGHT_MM, units = "mm", dpi = 600, bg = "white")
    writeLines(c("width_mm=80; height_mm=54; font=Arial",
      "font_pt: panel=10; group_heading=8; axis_title=7.5; ticks/legend_title=7; block_counts/legend_labels=6.5",
      "Figure4 planned left column: A=80x40 mm; gap=4 mm; C=80x54 mm; right B=94x98 mm; column_gap=4 mm",
      "intended_scale=1; assembly_and_caption_page_fit=NOT_CHECKED",
      "x_axis=Interpretable blocks (%); bar_height_y_units=0.72; row_step_y_units=1",
      paste("font_paths", paste(font_paths, collapse=";"))),
      file.path(out, "Figure4C_layout.txt"))
    writeLines(c(
      "Figure 4C. Credible-set sizes in exploratory molecular fine-mapping.",
      "Each bar represents all numerically interpretable expression+ASE blocks in one tissue and feature-origin group. Blocks are equally weighted; no PPA, GWAS, or LINK-positive preselection is applied.",
      "The axis label Interpretable blocks (%) refers to numerically interpretable blocks. Numerically interpretable denotes the upstream numerical QC classification, not evidence of a significant molecular association or validation of the model. A block is an upstream gene- and tissue-specific analysis unit; the same gene may contribute blocks in multiple tissues.",
      "The five mutually exclusive bins are 1-10, 11-50, 51-100, 101-500 and >500 members. Segment widths sum to 100%; boundaries correspond to cumulative thresholds. Numbers to the right are numerically interpretable block counts used as each bar's denominator, not dogs.",
      "Candidate loci and reference genes follow the upstream feature_origin annotation. CB: cerebellum; FC: frontal cortex; OC: occipital cortex; Pit: pituitary gland.",
      "The upstream model searched exactly one causal variant with target credible-set coverage 0.95, no null-causal configuration, no association/FDR screening, and global overdispersion set to zero. These are exploratory conditional posterior summaries, not significant eQTL counts or demonstrated causal accuracy.",
      "The 95% label denotes target cumulative posterior coverage under that model, not a 95% causal probability for each member. This panel summarizes the original upstream outputs without refitting the model.",
      "Effective sample counts and marker sets vary across blocks and tissues. Tissue distributions are descriptive and not paired tests; blocks from the same dog cohort are not independent biological replicates. Excluded blocks are separately listed, never treated as zero-sized sets.",
      "Absolute credible-set size depends on the input marker set and does not by itself establish better resolution in one group. The source tables also report each block's credible-set size divided by its input marker count, group medians of that fraction, and counts of blocks with tied highest PPA values; ties do not identify a unique top variant.",
      "Methods and visual reference: Wang et al., AJHG (2020), doi:10.1016/j.ajhg.2019.12.011, Figures 5 and 7. Structural adaptation of plot_thresh in austintwang/plasma run/interpret.py; original R implementation using disjoint bins instead of overlaid cumulative bars.",
      "Source: https://github.com/austintwang/plasma/blob/c5fe317d2a59814c0a9275c2b3216e3c978f7cc9/run/interpret.py"
    ), file.path(out, "Figure4C_caption.txt"))

    write_tsv(data.table(parameter = c("valid_input", "invalid_input", "bin_breaks", "width_mm", "height_mm", "threads", "rendered_QA"),
                         value = c(VALID, INVALID, "0,10,50,100,500,Inf", as.character(WIDTH_MM), as.character(HEIGHT_MM), "1", "pending_separate_review")), "Figure4C_parameters.tsv")
    capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
    log_stage("[4/4] completed normally; one shared coordinate panel, no independent panel scaling")
    writeLines(c("execution_status=success", paste0("output_dir=", out), paste0("master_log=", log_path)),
               file.path(out, "execution_completion.marker"))
    log_stage("task completed / done; execution_status=success; please check separately final dimensions, font, collision, source table with biological interpretation.")
  }, error = function(e) {
    log_stage(paste("execution_status=failed;", conditionMessage(e)))
    stop(e)
  })
}

main()
