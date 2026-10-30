#!/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/Rscript


Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           NUMEXPR_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
ROOT <- "/path/to/project/analysis_20260817/article/figure/figure5/B_candidate_partner_support"
INPUT_DIR <- file.path(ROOT, "preparation/20260929_231104")
OUT_ROOT <- file.path(ROOT, "plotting")
CELL_INPUT <- file.path(INPUT_DIR, "Figure5B_pair_sample_evidence.tsv")
PAIR_INPUT <- file.path(INPUT_DIR, "Figure5B_pair_catalogue.tsv")
CANDIDATE_INPUT <- file.path(INPUT_DIR, "Figure5B_selected_candidates.tsv")
DOGS <- c("Ace_396226", "CodRou_346627", "Crak_358336", "Dozer_386519", "SanMcg_349891")
DOG_LABELS <- c("Ace", "CodRou", "Crak", "Dozer", "SanMcg")
TISSUES <- c("cerebellum", "cortex_frontal", "cortex_occipital", "pituitary_gland")
TISSUE_LABELS <- c("CB", "FC", "OC", "Pit")
FONT <- "Arial"
WIDTH_MM <- 178

MAIN_ROW_MM <- 3.6
MAIN_MATRIX_LEFT_MM <- 74
MAIN_PARTNER_LEFT_MM <- 45.5
LEFT_MM <- 8
RIGHT_MM <- 2
LABEL_WIDTH_MM <- 46
DOG_GAP_MM <- 3
ROW_MM <- 4.2
GROUP_HEADER_MM <- 7
GROUP_GAP_MM <- 1
LEGEND_COLUMNS <- 4L
LEGEND_ROW_MM <- 5.5
LEGEND_ITEM_GAP_MM <- 6

MATRIX_LEFT_MM <- LEFT_MM + LABEL_WIDTH_MM
MATRIX_RIGHT_MM <- WIDTH_MM - RIGHT_MM
COLOURS <- c(link_detected = "#0072B2", A_link = "#0072B2", detected = "#777777", not_detected = "#FFFFFF",
             not_evaluable = "#EBEBEB",
             unknown = "#FFFFFF", outside_scope = "#FCFCFC")
STATE_LABELS <- c(link_detected = "Detected LINK", A_link = "Detected LINK",
                  detected = "Within-scope LINK not selected in A",
                  not_detected = "Evaluable; not detected", not_evaluable = "Not evaluable",
                  unknown = "Unresolved",
                  outside_scope = "Outside A tissue scope")

read_source <- function(path) {


  read.delim(path, sep = "\t", quote = "\"", comment.char = "", colClasses = "character",
             na.strings = "NA", check.names = FALSE, encoding = "UTF-8")
}

write_source <- function(d, path) {
  write.table(d, path, sep = "\t", row.names = FALSE, quote = TRUE, qmethod = "double",
              na = "NA", fileEncoding = "UTF-8")
}

log_stage <- function(message) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S ]"), message, "\n")
  flush.console()
}

make_cells <- function(cells, background = FALSE) {
  d <- cells
  d$plot_state <- d$pair_state
  if (background) {

    d$plot_state[d$pair_state == "detected"] <- "link_detected"
  } else {
    d$plot_state[!d$within_A_tissue_scope] <- "outside_scope"
    d$plot_state[d$A_selected_link] <- "A_link"
  }
  d$display_scope <- if (background) "cross_tissue_background" else "A_integration_scope"
  d$dog_index <- match(d$animal_id, DOGS)
  d$tissue_index <- match(d$tissue, TISSUES)
  d
}

layout_rows <- function(pairs, candidates, legend_count, compact = FALSE) {
  p <- pairs[order(match(pairs$new_gene_id, candidates$new_gene_id), pairs$partner_symbol,
                   pairs$partner_gene_id), , drop = FALSE]
  if (compact) {

    cursor <- 14
    p$row_top_mm <- NA_real_
    p$group_top_mm <- NA_real_
    for (gene in unique(p$new_gene_id)) {
      ii <- which(p$new_gene_id == gene)
      p$group_top_mm[ii] <- cursor
      p$row_top_mm[ii] <- cursor + (seq_along(ii) - 1) * MAIN_ROW_MM
      cursor <- cursor + length(ii) * MAIN_ROW_MM + GROUP_GAP_MM
    }
    footer_top <- max(p$row_top_mm) + 4.7
    return(list(pairs = p, height = footer_top + 3.5,
                footer_top = footer_top, compact = TRUE, row_mm = MAIN_ROW_MM))
  }
  cursor <- 23
  p$row_top_mm <- NA_real_
  p$group_top_mm <- NA_real_
  for (gene in unique(p$new_gene_id)) {
    ii <- which(p$new_gene_id == gene)
    p$group_top_mm[ii] <- cursor
    cursor <- cursor + if (compact) 0 else GROUP_HEADER_MM
    p$row_top_mm[ii] <- cursor + (seq_along(ii) - 1) * ROW_MM
    cursor <- cursor + length(ii) * ROW_MM + GROUP_GAP_MM
  }

  footer_top <- cursor + 2
  legend_rows <- ceiling(legend_count / LEGEND_COLUMNS)
  note_top <- footer_top + (legend_rows - 1) * LEGEND_ROW_MM + 7
  list(pairs = p, height = note_top - 2, footer_top = footer_top, note_top = note_top,
       compact = compact, row_mm = ROW_MM)
}

draw_matrix <- function(d, layout, columns, candidates, background) {
  grid::grid.newpage()
  h <- layout$height
  p <- layout$pairs
  cell_step <- layout$cell_step
  text <- function(label, x, y, size = 7, bold = FALSE, colour = "#252525", just = "left") {
    grid::grid.text(label, x = grid::unit(x, "mm"), y = grid::unit(h - y, "mm"), just = just,
      gp = grid::gpar(fontfamily = FONT, fontsize = size, fontface = if (bold) "bold" else "plain", col = colour))
  }
  line <- function(x, y, colour = "#D0D0D0", width = 0.5) {
    grid::grid.lines(x = grid::unit(x, "mm"), y = grid::unit(h - y, "mm"),
                     gp = grid::gpar(col = colour, lwd = width))
  }
  glyph <- function(state, x, y, width = cell_step - 0.35, height = layout$row_mm - 0.55) {
    border <- if (state %in% c("A_link", "link_detected")) COLOURS[[state]] else
      if (state == "outside_scope") "#EEEEEE" else "#D6D6D6"
    grid::grid.rect(x = grid::unit(x, "mm"), y = grid::unit(h - y, "mm"),
      width = grid::unit(width, "mm"), height = grid::unit(height, "mm"),
      gp = grid::gpar(fill = COLOURS[[state]], col = border, lwd = 0.4))
    if (state == "outside_scope") line(c(x - 0.9, x + 0.9), c(y, y), "#AAAAAA")
    if (state == "unknown") text("?", x, y, size = 7, colour = "#777777", just = "centre")
  }

  if (layout$compact) {
    text("B", 2, 4, 10, TRUE)
    text("Candidate", LEFT_MM, 8.5, 7, TRUE)
    text("Coding partner", MAIN_PARTNER_LEFT_MM, 8.5, 7, TRUE)
    for (gene in unique(p$new_gene_id)) {
      ii <- which(p$new_gene_id == gene)
      text(gene, LEFT_MM, mean(p$row_top_mm[ii]), 6.5, TRUE)
      if (min(ii) > 1) {
        yline <- mean(p$row_top_mm[c(min(ii) - 1, min(ii))])
        line(c(LEFT_MM, MATRIX_RIGHT_MM), c(yline, yline), "#E4E4E4")
      }
    }
  } else {
    text("Candidate / coding partner", LEFT_MM, 16, 7, TRUE)
  }
  for (dog in unique(columns$dog_index)) {
    cc <- columns[columns$dog_index == dog, , drop = FALSE]
    xs <- cc$plot_x_mm
    text(DOG_LABELS[dog], mean(xs), if (layout$compact) 4 else 11, 7, TRUE, just = "centre")
    for (j in seq_len(nrow(cc))) text(TISSUE_LABELS[cc$tissue_index[j]], xs[j], if (layout$compact) 8.5 else 16, 7, just = "centre")
    line(c(min(xs) - cell_step / 2, max(xs) + cell_step / 2), rep(if (layout$compact) 11 else 20, 2))
  }
  for (gene in if (layout$compact) character(0) else unique(p$new_gene_id)) {
    first <- which(p$new_gene_id == gene)[1]
    a <- candidates[match(gene, candidates$new_gene_id), , drop = FALSE]
    tissue <- if (a$evidence_route == "Molecular") TISSUE_LABELS[match(a$eqtl_tissue, TISSUES)] else NULL
    top <- p$group_top_mm[first]

    text(gene, LEFT_MM, top + 1.5, 7, TRUE)

    if (!background) {
      description <- paste(c(a$evidence_route, a$inclusion_status, tissue), collapse = " | ")
      text(description, MATRIX_LEFT_MM, top + 1.5, 7, colour = "#555555")
    }
    line(c(LEFT_MM, MATRIX_RIGHT_MM), c(top + 4, top + 4), "#E4E4E4")
  }
  for (i in seq_len(nrow(p))) {
    yy <- p$row_top_mm[i]
    text(p$partner_symbol[i], if (layout$compact) MAIN_PARTNER_LEFT_MM else LEFT_MM + 3, yy, 7)
    z <- d[d$pair_id == p$pair_id[i], , drop = FALSE]
    for (j in seq_len(nrow(z))) glyph(z$plot_state[j], z$plot_x_mm[j], yy)
  }
  legend_states <- names(STATE_LABELS)[names(STATE_LABELS) %in% d$plot_state]
  legend_labels <- STATE_LABELS

  label_widths <- vapply(legend_labels[legend_states], function(label) {
    grid::convertWidth(grid::grobWidth(grid::textGrob(label,
      gp = grid::gpar(fontfamily = FONT, fontsize = 7))), "mm", valueOnly = TRUE)
  }, numeric(1))
  legend_left <- LEFT_MM
  for (i in seq_along(legend_states)) {
    if ((i - 1) %% LEGEND_COLUMNS == 0) legend_left <- LEFT_MM
    xx <- legend_left + 2
    yy <- layout$footer_top + ((i - 1) %/% LEGEND_COLUMNS) * LEGEND_ROW_MM
    glyph(legend_states[i], xx, yy, 4, 3)
    text(legend_labels[[legend_states[i]]], xx + 4, yy, 7)
    legend_left <- legend_left + 4 + 2 + label_widths[i] + LEGEND_ITEM_GAP_MM
  }

}

export_matrix <- function(d, pairs, candidates, out, stem, background = FALSE) {
  layout <- layout_rows(pairs, candidates, sum(names(STATE_LABELS) %in% d$plot_state),
                        compact = background)

  columns <- unique(d[, c("animal_id", "tissue", "dog_index", "tissue_index"), drop = FALSE])
  columns <- columns[order(columns$dog_index, columns$tissue_index), , drop = FALSE]
  group_index <- match(columns$dog_index, unique(columns$dog_index))

  matrix_left <- if (layout$compact) MAIN_MATRIX_LEFT_MM else MATRIX_LEFT_MM
  layout$cell_step <- (MATRIX_RIGHT_MM - matrix_left -
    (length(unique(group_index)) - 1) * DOG_GAP_MM) / nrow(columns)
  columns$plot_x_mm <- matrix_left + (seq_len(nrow(columns)) - 0.5) * layout$cell_step +
    (group_index - 1) * DOG_GAP_MM
  d$plot_x_mm <- columns$plot_x_mm[match(paste(d$animal_id, d$tissue), paste(columns$animal_id, columns$tissue))]
  d$row_index <- match(d$pair_id, layout$pairs$pair_id)
  d <- d[order(d$row_index, d$dog_index, d$tissue_index), , drop = FALSE]
  prefix <- file.path(out, stem)

  for (format in c("pdf", "svg", "png", "tiff")) {
    if (format == "pdf") {
      grDevices::cairo_pdf(paste0(prefix, ".pdf"), width = WIDTH_MM / 25.4,
                          height = layout$height / 25.4, family = FONT, bg = "white")
    } else if (format == "svg") {
      svglite::svglite(paste0(prefix, ".svg"), width = WIDTH_MM / 25.4,
                      height = layout$height / 25.4, bg = "white")
    } else if (format == "png") {
      grDevices::png(paste0(prefix, ".png"), width = WIDTH_MM, height = layout$height,
                    units = "mm", res = 600, type = "cairo", bg = "white")
    } else {
      grDevices::tiff(paste0(prefix, ".tiff"), width = WIDTH_MM, height = layout$height,
                     units = "mm", res = 600, type = "cairo", compression = "lzw", bg = "white")
    }
    tryCatch(draw_matrix(d, layout, columns, candidates, background), finally = grDevices::dev.off())
  }
  writeLines(c(paste0("width_mm=", WIDTH_MM, "; height_mm=", layout$height),
    "font=Arial; panel_pt=10; ordinary_pt=7; full_candidate_id_pt=6.5",
    "intended_placement_width_mm=178; intended_scale=1; raster_ppi=600; TIFF=LZW",
    "PDF_font_embedding=NOT_CHECKED; visual_layout=NOT_CHECKED; assembled_figure_caption=NOT_CHECKED",
    "Dimensions and font sizes are specified by the script; after execution, check PDF font embedding, text boundaries, and the page space occupied by the full Figure 5 caption."),
    paste0(prefix, "_layout.txt"))
  write_source(d, paste0(prefix, "_cells.tsv"))
  write_source(layout$pairs, paste0(prefix, "_rows.tsv"))
  write_source(columns, paste0(prefix, "_columns.tsv"))
  log_stage(paste(stem, "pair count=", nrow(pairs), "display sample count=", nrow(columns),
                  "number of displayed cells=", nrow(d), "width mm=", WIDTH_MM, "height mm=", layout$height))
}

main <- function() {
  out <- file.path(OUT_ROOT, format(Sys.time(), "%Y%m%d_%H%M%S"))
  if (!dir.create(out, showWarnings = FALSE)) stop("cannot create run directory, or a directory with the same name already exists: ", out)
  log_path <- file.path(out, "plot_figure5B.log")
  con <- file(log_path, "wt")
  sink(con, type = "output")
  sink(con, type = "message")
  on.exit({sink(type = "message"); sink(type = "output"); close(con)}, add = TRUE)
  tryCatch({
    log_stage("[1/4] reading existing prepared source tables; concurrent_tasks=1 threads_per_task=1")
    cat("LC_CTYPE=", Sys.getlocale("LC_CTYPE"), "\n", sep = "")
    cat("CELL_INPUT=", CELL_INPUT, "\nPAIR_INPUT=", PAIR_INPUT,
        "\nCANDIDATE_INPUT=", CANDIDATE_INPUT, "\nOUTPUT=", out, "\n", sep = "")

    font_paths <- systemfonts::match_fonts(c(FONT, FONT), weight = c("normal", "bold"))$path
    expected_fonts <- c("/path/to/fonts/Arial/Arial.TTF",
                        "/path/to/fonts/Arial/Arialbd.TTF")
    if (!all(file.exists(expected_fonts)) ||
        !identical(normalizePath(font_paths), normalizePath(expected_fonts)))
      stop("Arial regular or bold font is unavailable; do not silently substitute fonts.")
    log_stage(paste("Arial font: ", paste(font_paths, collapse = "; ")))
    cells <- read_source(CELL_INPUT)
    pairs <- read_source(PAIR_INPUT)
    candidates <- read_source(CANDIDATE_INPUT)
    for (field in c("within_A_tissue_scope", "A_selected_link", "sample_present", "background_only_link")) {
      cells[[field]] <- cells[[field]] == "True"
    }

    sample_key <- paste(cells$animal_id, cells$tissue, sep = "|")
    missing_keys <- setdiff(unique(sample_key), unique(sample_key[cells$sample_present %in% TRUE]))
    omitted <- unique(cells[sample_key %in% missing_keys, c("animal_id", "tissue"), drop = FALSE])
    omitted$omission_reason <- rep("sample_unavailable", nrow(omitted))
    write_source(omitted, file.path(out, "Figure5B_omitted_samples.tsv"))
    cells <- cells[!sample_key %in% missing_keys, , drop = FALSE]
    omission_note <- if (nrow(omitted)) paste0("Unavailable samples omitted from both displays: ",
      paste(paste(omitted$animal_id, omitted$tissue, sep = " / "), collapse = "; "),
      ". Complete preparation tables retain these missing-sample records.") else "No unavailable sample columns were omitted."
    log_stage(omission_note)
    log_stage("[2/4] supplement matrix: all pair, by A integrate span display; supplementary figure numbering remains to be finalized")

    scope_cells <- make_cells(cells)
    page_genes <- list()
    current <- character(0)
    for (gene in unique(candidates$new_gene_id)) {
      if (!gene %in% pairs$new_gene_id) next
      proposed <- c(current, gene)
      pp <- pairs[pairs$new_gene_id %in% proposed, , drop = FALSE]
      dd <- scope_cells[scope_cells$new_gene_id %in% proposed, , drop = FALSE]
      proposed_height <- layout_rows(pp, candidates,
        sum(names(STATE_LABELS) %in% dd$plot_state))$height
      if (length(current) && proposed_height > 200) {
        page_genes[[length(page_genes)+1L]] <- current
        current <- gene
      } else current <- proposed
    }
    if (length(current)) page_genes[[length(page_genes)+1L]] <- current
    for (page in seq_along(page_genes)) {
      genes <- page_genes[[page]]
      export_matrix(scope_cells[scope_cells$new_gene_id %in% genes, , drop = FALSE],
        pairs[pairs$new_gene_id %in% genes, , drop = FALSE], candidates, out,
        sprintf("Figure5B_all_pairs_A_scope_part%02d", page))
    }
    log_stage("[3/4] main figure B: display all coding partners and available samples for the primary molecular integration candidates")

    main_genes <- candidates$new_gene_id[candidates$evidence_route %in% "Molecular" &
                                         candidates$inclusion_status %in% "Primary"]
    main_pairs <- pairs$pair_id[pairs$new_gene_id %in% main_genes]
    if (length(main_pairs)) {
      export_matrix(make_cells(cells[cells$pair_id %in% main_pairs, , drop = FALSE], TRUE),
        pairs[pairs$pair_id %in% main_pairs, , drop = FALSE], candidates, out,
        "Figure5B_molecular_partner_support", TRUE)
    } else {
      log_stage("the preparation table has no coding partners for primary molecular integration candidates; do not generate an empty main figure B.")
    }
    log_stage("[4/4] save figure legend, input record and execution status")
    write_source(data.frame(role = c("pair_sample_evidence", "pair_catalogue", "selected_candidates"),
      absolute_path = c(CELL_INPUT, PAIR_INPUT, CANDIDATE_INPUT)), file.path(out, "Figure5B_input_paths.tsv"))

    common_caption <- c(
      "Full DOG5_NEWGENE_ identifiers identify candidate groups. The main panel uses separate candidate and coding-partner columns; supplementary matrices use candidate group headings. CB, cerebellum; FC, frontal cortex; OC, occipital cortex; Pit, pituitary gland.",
      "Each cell represents a pair in one dog-tissue sample. White denotes evaluable but undetected; light grey denotes non-evaluable. An unresolved status, if present, is marked with a question mark.",
      omission_note,
      "Eligibility is inherited from preparation: both specific endpoints require at least 10 site-level allele supports and a same-autosome 500-kb flanking window. Allelic bias is not an eligibility filter. The plotting script does not recompute counts or LINK calls.",
      "Colours do not encode linkage strength or regulatory direction. Repeated tissues from one dog are not independent animals. Reported mechanisms, multigene attribution review flags, and phase-comparability status are retained in the source tables; these records do not establish causal regulation.",
      "Existing molecular implementation limitations are not resolved by this visualization.")
    writeLines(c(
      "Supplementary matrix: complete candidate-protein-coding partner support within the Figure 5A integration scope. The supplementary figure number will be assigned during manuscript assembly.",
      "All selected pairs are retained across the numbered part files, split between candidates to keep each current-data page below 200 mm high. Each part repeats sample headers and its applicable legend. Pairs are grouped in A candidate order and sorted alphabetically by partner within each candidate. Blue denotes LINK selected in A. A dash denotes outside A tissue scope, not absence of LINK. Dark grey, if present, denotes detected within-scope LINK not selected in A.",
      "Molecular candidates retain their A tissue scope; positional candidates retain their multi-tissue scope. Primary and Exploratory labels inherit A inclusion status. This matrix provides the complete context for the molecular candidates displayed in the main panel.",
      common_caption
    ), file.path(out, "Figure5B_all_pairs_A_scope_caption.txt"))
    if (length(main_pairs)) {
      display_pairs <- pairs[pairs$pair_id %in% main_pairs, , drop = FALSE]

      main_candidates <- candidates[candidates$new_gene_id %in% main_genes, , drop = FALSE]
      a_tissue_note <- paste0("The molecular-integration tissue in panel A is: ",
        paste(paste(main_candidates$new_gene_id,
                    TISSUE_LABELS[match(main_candidates$eqtl_tissue, TISSUES)], sep = " / "),
              collapse = "; "), ".")
      writeLines(c(
        "Figure 5B. Protein-coding partner LINK support for primary molecular-integration candidates.",
        paste0("Displayed pair(s): ", paste(paste(display_pairs$new_gene_id,
          display_pairs$partner_symbol, sep = " - "), collapse = "; "), "."),
        "All protein-coding partner pairs of candidates labelled Molecular and Primary in panel A are displayed, irrespective of additional-tissue LINK detection. Rows are grouped in A candidate order and sorted alphabetically by partner; columns show all available dog-tissue samples.",
        "Blue denotes detected LINK in any displayed tissue, regardless of whether that record contributed to panel A. A_selected_link and within_A_tissue_scope retain this provenance in the cell source table.",
        a_tissue_note,
        "Additional tissue LINK does not establish molecular-GWAS integration in those tissues and does not increase the dog count for the original A tissue-specific record. The complete A-scope pair matrix is provided separately as a supplementary figure.",
        common_caption
      ), file.path(out, "Figure5B_molecular_partner_support_caption.txt"))
    }
    capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
    writeLines(c("execution_status=success", paste0("output_dir=", out), paste0("master_log=", log_path)),
               file.path(out, "execution_completion.marker"))
    log_stage("task completed / done; execution_status=success")
    log_stage("please check the figure layout separately, source table, status count and biological interpretation; this script does not perform post-run quality auditing.")
  }, error = function(e) {
    log_stage(paste("execution_status=failed;", conditionMessage(e)))
    stop(e)
  })
}

main()
