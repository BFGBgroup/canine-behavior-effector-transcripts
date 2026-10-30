#!/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/Rscript


Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
ROOT <- "/path/to/project/analysis_20260817"
OUT_ROOT <- file.path(ROOT, "article/figure/figure4/B_candidate_pair_support")
LINK_ROOT <- file.path(ROOT, "5-dogs/7.Allelome.LINK")
ENDPOINT_RUN <- file.path(LINK_ROOT, "1.2_build_dog5_allelome_locus_tables/20260907_202637")
ANNOTATION_RUN <- file.path(LINK_ROOT, "1.4_annotate_and_merge_dog5_allelome_links/20260907_213429")
MANIFEST <- file.path(ENDPOINT_RUN, "dog5_allelome_locus_table_manifest.tsv")
LINKS <- file.path(ANNOTATION_RUN, "dog5_allelome_links.annotated.tsv")
ANNOTATION_TASKS <- file.path(ANNOTATION_RUN, "dog5_allelome_annotation.task_manifest.tsv")
SUCCESS <- file.path(LINK_ROOT, "scripts/1.3_run_dog5_allelome_link/log/1.3_run_dog5_allelome_link.20260907_205139.success.tsv")
MIN_SUPPORT <- 10
WINDOW_BP <- 500000
BIAS_UPPER <- 0.6
BIAS_LOWER <- 1 - BIAS_UPPER
MAX_DISPLAY_PAIRS <- 20L

WIDTH_MM <- 94
HEIGHT_MM <- 98
FONT <- "Arial"
DOGS <- c("Ace_396226", "CodRou_346627", "Crak_358336", "Dozer_386519", "SanMcg_349891")
TISSUES <- c("cerebellum", "cortex_frontal", "cortex_occipital", "pituitary_gland")
TISSUE_LABELS <- c("CB", "FC", "OC", "Pit")

main <- function() {
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out <- file.path(OUT_ROOT, stamp)
  if (!dir.create(out, recursive = TRUE, showWarnings = FALSE)) stop("output directory already exists or cannot create: ", out)
  log_dir <- file.path(OUT_ROOT, "log")
  dir.create(log_dir, showWarnings = FALSE)
  log_path <- file.path(log_dir, paste0("plot_figure4B.", stamp, ".", Sys.getpid(), ".log"))
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
    read_tsv <- function(p) fread(file = p, sep = "\t", na.strings = "NA", showProgress = FALSE)
    write_tsv <- function(x, name) fwrite(x, file.path(out, name), sep = "\t", na = "NA")
    collapse_values <- function(x) paste(sort(unique(x[!is.na(x) & nzchar(x)])), collapse = ";")
    bias <- function(x) x >= BIAS_UPPER | x <= BIAS_LOWER
    log_stage("[1/5] read endpoint list and actual link; concurrent_tasks=1 threads_per_task=1")
    cat("MANIFEST=", MANIFEST, "\nLINKS=", LINKS, "\nSUCCESS=", SUCCESS,
        "\nANNOTATION_TASKS=", ANNOTATION_TASKS, "\nOUTPUT=", out, "\n", sep = "")
    manifest <- read_tsv(MANIFEST)
    manifest <- manifest[chromosome %in% paste0("chr", 1:38)]
    success <- read_tsv(SUCCESS)
    annotated_tasks <- read_tsv(ANNOTATION_TASKS)

    si <- match(manifest$task_id, success$task_id)
    ai <- match(manifest$task_id, annotated_tasks$task_id)
    if (anyNA(si) || anyNA(ai) || any(success$execution_status[si] != "success") ||
        any(success$input_file[si] != manifest$locus_table) ||
        any(annotated_tasks$source_metadata[ai] != manifest$locus_metadata)) {
      stop("upstream success record, annotation tasks do not correspond to the endpoint list; missing tasks must not be treated as zero links.")
    }
    links <- read_tsv(LINKS)
    links <- links[pair_type == "candidate_new_gene--known_protein_coding"]
    if (!nrow(links)) stop("The current input contains no candidate–protein-coding gene pairs; cannot define the set of previously detected pairs.")
    if (any(!links$task_id %in% manifest$task_id)) stop("link task exceeds current endpoint list.")

    candidate_is_base <- links$base_endpoint_type == "candidate_new_gene"
    links[, `:=`(
      candidate_gene_id = fifelse(candidate_is_base, base_gene_id, target_gene_id),
      partner_gene_id = fifelse(candidate_is_base, target_gene_id, base_gene_id),
      candidate_endpoint_id = fifelse(candidate_is_base, name_base, name_target),
      partner_endpoint_id = fifelse(candidate_is_base, name_target, name_base),
      chromosome = fifelse(candidate_is_base, chr_base, chr_target),
      partner_symbol = fifelse(candidate_is_base, target_gene_symbol, base_gene_symbol),
      candidate_class_code = fifelse(candidate_is_base, base_class_code, target_class_code),
      candidate_novelty = fifelse(candidate_is_base, base_novelty_category, target_novelty_category))]
    links[, pair_id := paste(candidate_gene_id, partner_gene_id, sep = "::")]
    catalog <- links[, .(chromosome = collapse_values(chromosome),
      partner_symbol = collapse_values(partner_symbol),
      candidate_class_code = collapse_values(candidate_class_code),
      candidate_novelty = collapse_values(candidate_novelty)),
      by = .(pair_id, candidate_gene_id, partner_gene_id)]
    if (any(grepl(";", catalog$chromosome, fixed = TRUE))) stop("the same pair corresponding multiple chromosome, requires first specify endpoint identity.")
    observed <- links[, .(observed_link_rows = .N,
      noflag_link_rows = sum(evidence_stratum == "nonoverlap_no_multigene_flag"),
      evidence_strata = collapse_values(evidence_stratum),
      reported_mechanisms = collapse_values(mechanism),
      observed_candidate_endpoint = collapse_values(candidate_endpoint_id),
      observed_partner_endpoint = collapse_values(partner_endpoint_id),
      source_link_files = collapse_values(source_link_file),
      source_link_lines = collapse_values(as.character(source_link_line))),
      by = .(animal_id, tissue, pair_id)]
    log_stage(paste("previously detect pair count: ", nrow(catalog), "; read all pair, main figure at most display", MAX_DISPLAY_PAIRS, "row"))

    log_stage("[2/5] by two specific endpoint reconstruct 20 dogs—tissue combination of pair eligibility")
    samples <- CJ(animal_id = DOGS, tissue = TISSUES)
    evidence <- vector("list", nrow(samples))
    for (i in seq_len(nrow(samples))) {
      animal <- samples$animal_id[i]
      tissue_name <- samples$tissue[i]
      tasks <- manifest[animal_id == animal & tissue == tissue_name]
      d <- copy(catalog)
      d[, `:=`(animal_id = animal, tissue = tissue_name, sample_present = nrow(tasks) > 0L)]
      endpoint_parts <- vector("list", nrow(tasks))
      for (j in seq_len(nrow(tasks))) {
        task <- tasks[j]
        ep <- read_tsv(task$locus_table)
        meta <- read_tsv(task$locus_metadata)
        mi <- match(ep$name, meta$group_id)
        if (anyNA(mi) || anyDuplicated(meta$group_id)) stop("metadata endpoint matching is not unique: ", task$task_id)
        ep[, `:=`(gene_id = meta$gene_id[mi], record_type = meta$record_type[mi],
          gene_biotype = meta$gene_biotype[mi], multi_gene_variant_count = meta$multi_gene_variant_count[mi],
          variant_count = meta$variant_count[mi], gene_strand = meta$gene_strand[mi],
          source_gene_file = meta$input_file[mi], source_gene_line = meta$input_line[mi],
          count_unit = meta$count_unit[mi], source_locus_table = task$locus_table,
          source_metadata = task$locus_metadata, task_id = task$task_id)]
        endpoint_parts[[j]] <- ep[gene_id %in% c(catalog$candidate_gene_id, catalog$partner_gene_id)]
      }

      if (nrow(tasks) && !setequal(tasks$chromosome, paste0("chr", 1:38))) stop("sample-chromosome tasks are incomplete: ", animal, " ", tissue_name)
      fields <- c("name", "start", "end", "A1_reads", "A2_reads", "total_reads", "allelic_ratio",
        "record_type", "gene_biotype", "variant_count", "multi_gene_variant_count", "gene_strand",
        "count_unit", "source_gene_file", "source_gene_line", "source_locus_table", "source_metadata", "task_id")
      if (nrow(tasks)) {
        ep <- rbindlist(endpoint_parts, use.names = TRUE)
        ep[, endpoint_key := paste(chr, gene_id, sep = "::")]
        if (anyDuplicated(ep$endpoint_key)) stop("Multiple endpoints share the same chromosome and gene ID within a sample; eligibility and LINK evidence cannot be combined across endpoints.")
        ci <- match(paste(d$chromosome, d$candidate_gene_id, sep = "::"), ep$endpoint_key)
        pi <- match(paste(d$chromosome, d$partner_gene_id, sep = "::"), ep$endpoint_key)
        for (nm in fields) {
          set(d, j = paste0("candidate_", nm), value = ep[[nm]][ci])
          set(d, j = paste0("partner_", nm), value = ep[[nm]][pi])
        }
      } else {
        for (nm in fields) {
          set(d, j = paste0("candidate_", nm), value = NA)
          set(d, j = paste0("partner_", nm), value = NA)
        }
      }
      d[, `:=`(candidate_present = !is.na(candidate_name), partner_present = !is.na(partner_name))]

      d[, window_pass := candidate_present & partner_present &
        partner_end >= candidate_start - WINDOW_BP & partner_start <= candidate_end + WINDOW_BP]
      d[, endpoint_types_pass := candidate_record_type == "manual_new_gene" &
        !partner_record_type %in% c("manual_new_gene", "unannotated_variant") & partner_gene_biotype == "protein_coding"]
      d[, pair_eligible := sample_present & candidate_present & partner_present &
        candidate_total_reads >= MIN_SUPPORT & partner_total_reads >= MIN_SUPPORT & window_pass & endpoint_types_pass]
      d[is.na(pair_eligible), pair_eligible := FALSE]
      d[, ineligible_reason := fcase(!sample_present, "sample_missing",
        !candidate_present & !partner_present, "both_endpoints_missing",
        !candidate_present, "candidate_endpoint_missing", !partner_present, "partner_endpoint_missing",
        !is.na(endpoint_types_pass) & !endpoint_types_pass, "endpoint_type_mismatch",
        candidate_total_reads < MIN_SUPPORT & partner_total_reads < MIN_SUPPORT, "both_support_below_10",
        candidate_total_reads < MIN_SUPPORT, "candidate_support_below_10",
        partner_total_reads < MIN_SUPPORT, "partner_support_below_10",
        !is.na(window_pass) & !window_pass, "outside_window", pair_eligible, "eligible", default = "missing_endpoint_values")]
      d[, `:=`(candidate_bias_pass = bias(candidate_allelic_ratio), partner_bias_pass = bias(partner_allelic_ratio))]
      d <- merge(d, observed[animal_id == animal & tissue == tissue_name], by = c("animal_id", "tissue", "pair_id"), all.x = TRUE, sort = FALSE)
      d[is.na(observed_link_rows), observed_link_rows := 0L]
      d[is.na(noflag_link_rows), noflag_link_rows := 0L]

      pos <- which(d$observed_link_rows > 0L)
      if (length(pos) && (any(!d$pair_eligible[pos]) || anyNA(d$candidate_bias_pass[pos]) ||
          anyNA(d$partner_bias_pass[pos]) || any(!d$candidate_bias_pass[pos] | !d$partner_bias_pass[pos]) ||
          any(d$observed_candidate_endpoint[pos] != d$candidate_name[pos]) ||
          any(d$observed_partner_endpoint[pos] != d$partner_name[pos]))) {
        stop("Detected positive pair is inconsistent with the eligibility, bias or identity of its two endpoints: ", animal, " ", tissue_name)
      }
      d[, `:=`(link_detected = fifelse(pair_eligible, observed_link_rows > 0L, NA),
        link_noflag = fifelse(pair_eligible, noflag_link_rows > 0L, NA),
        phase_comparability = "not_assessed_in_this_step")]
      d[, pair_state := fifelse(!pair_eligible, "not_evaluable", fifelse(link_detected, "detected", "not_detected"))]
      evidence[[i]] <- d
      log_stage(paste(i, "/", nrow(samples), animal, tissue_name, "; sample exists: ", nrow(tasks) > 0L))
    }
    evidence <- rbindlist(evidence, use.names = TRUE)

    log_stage("[3/5] summary each pair of per tissue dog count, save complete underlying table and display sort")
    summary <- evidence[, .(n_sample_present = sum(sample_present), n_evaluable = sum(pair_eligible),
      n_linked = sum(link_detected, na.rm = TRUE), n_noflag = sum(link_noflag, na.rm = TRUE),
      evaluable_dogs = collapse_values(animal_id[pair_eligible]),
      linked_dogs = collapse_values(animal_id[which(link_detected)]),
      noflag_dogs = collapse_values(animal_id[which(link_noflag)])), by = .(pair_id, tissue)]
    summary[, linked_pct := fifelse(n_evaluable > 0L, 100 * n_linked / n_evaluable, NA_real_)]
    summary[, tissue_index := match(tissue, TISSUES)]
    ranking <- summary[, .(max_noflag_in_one_tissue = max(n_noflag),
      max_linked_in_one_tissue = max(n_linked), n_tissues_detected = sum(n_linked > 0L),
      evaluable_dog_tissue_count = sum(n_evaluable)), by = pair_id]
    dog_support <- evidence[, .(n_distinct_linked_dogs = uniqueN(animal_id[which(link_detected)])), by = pair_id]
    ranking <- merge(merge(catalog, ranking, by = "pair_id"), dog_support, by = "pair_id")

    setorderv(ranking, c("max_noflag_in_one_tissue", "max_linked_in_one_tissue", "n_distinct_linked_dogs",
      "n_tissues_detected", "evaluable_dog_tissue_count", "candidate_gene_id", "partner_gene_id"),
      c(-1L, -1L, -1L, -1L, -1L, 1L, 1L))
    ranking[, `:=`(display_rank = seq_len(.N), selected_for_plot = seq_len(.N) <= MAX_DISPLAY_PAIRS)]

    selected <- ranking[selected_for_plot == TRUE]
    write_tsv(evidence, "Figure4B_pair_dog_evidence.tsv")
    write_tsv(summary, "Figure4B_pair_tissue_summary.tsv")
    write_tsv(ranking, "Figure4B_pair_ranking.tsv")
    write_tsv(links, "Figure4B_observed_links_source.tsv")
    write_tsv(manifest, "Figure4B_input_manifest.tsv")

    log_stage("[4/5] plot pair×tissue support matrix; numbers indicate the count of dogs with detected links/evaluable dog count")
    plot_data <- merge(summary, selected, by = "pair_id")
    plot_data[, `:=`(plot_y = nrow(selected) + 1L - display_rank,
      fill_level = fifelse(n_evaluable > 0L, as.character(n_linked), "NA"),
      cell_label = fifelse(n_evaluable > 0L, paste0(n_linked, "/", n_evaluable,
        fifelse(n_linked > n_noflag, "*", "")), "NA"),
      text_colour = fifelse(n_evaluable > 0L & n_linked >= 3L, "white", "#222222"))]
    selected[, partner_label := fifelse(nzchar(partner_symbol), partner_symbol, partner_gene_id)]
    selected[, row_label := paste0(candidate_gene_id, " \u2013 ", partner_label,
      fifelse(nzchar(candidate_class_code), paste0(" [", candidate_class_code, "]"), ""))]


    selected[, candidate_display := paste0(sub("^DOG5_NEWGENE_", "", candidate_gene_id),
      fifelse(nzchar(candidate_class_code), paste0(" [", candidate_class_code, "]"), ""))]
    headers <- data.frame(x = c(-3.3, -1.6, 1:4), y = nrow(selected) + 1,
      label = c("Candidate", "Coding partner", TISSUE_LABELS),
      hjust = c(0, 0, rep(0.5, 4)))
    colours <- c("0" = "#FFFFFF", "1" = "#DEEBF7", "2" = "#9ECAE1", "3" = "#4292C6",
                 "4" = "#2171B5", "5" = "#084594", "NA" = "#D9D9D9")

    legend_levels <- names(colours)[names(colours) %in% plot_data$fill_level]
    legend_labels <- ifelse(legend_levels == "NA", "Not evaluable", legend_levels)


    panel_caption <- "Cells: detected/evaluable dogs."
    if (any(plot_data$n_evaluable > 0L & plot_data$n_linked > plot_data$n_noflag)) {
      panel_caption <- paste0(panel_caption,
        "\n* At least one detected dog has flagged-only links.")
    }
    p <- ggplot(plot_data, aes(x = tissue_index, y = plot_y)) +
      geom_tile(aes(fill = fill_level), colour = "#E5E5E5", linewidth = 0.25, width = 0.96, height = 0.94) +
      geom_text(aes(label = cell_label, colour = text_colour), size = 7, size.unit = "pt", family = FONT) +
      geom_text(data = selected, aes(x = -3.3, y = nrow(selected) + 1L - display_rank,
        label = candidate_display), inherit.aes = FALSE, hjust = 0,
        size = 6.5, size.unit = "pt", family = FONT, colour = "black") +
      geom_text(data = selected, aes(x = -1.6, y = nrow(selected) + 1L - display_rank,
        label = partner_label), inherit.aes = FALSE, hjust = 0,
        size = 6.5, size.unit = "pt", family = FONT, colour = "black") +
      geom_text(data = headers, aes(x = x, y = y, label = label, hjust = hjust),
        inherit.aes = FALSE, size = 7, size.unit = "pt", family = FONT,
        fontface = "bold", colour = "black") +

      annotate("text", x = -3.3, y = nrow(selected) + 2.05,
        label = "Candidate prefix: DOG5_NEWGENE_", hjust = 0,
        size = 6.5, size.unit = "pt", family = FONT, colour = "#444444") +
      scale_colour_identity() +
      scale_fill_manual(values = colours, breaks = legend_levels, limits = legend_levels, drop = FALSE,
        name = "Dogs with a detected link", labels = legend_labels,
        guide = guide_legend(nrow = 1, direction = "horizontal",
          theme = theme(legend.title.position = "top"))) +
      scale_x_continuous(breaks = NULL, limits = c(-3.4, 4.5), expand = c(0, 0)) +
      scale_y_continuous(breaks = NULL,
        limits = c(0.5, nrow(selected) + 2.65), expand = c(0, 0)) +
      labs(x = NULL, y = NULL, tag = "B",
        caption = panel_caption) +
      theme_minimal(base_size = 8, base_family = FONT) +
      theme(panel.grid = element_blank(), axis.text.x = element_text(colour = "black", size = 8),
        axis.text.y = element_text(colour = "black", size = 7),
        legend.position = "bottom", legend.justification = "left", legend.title = element_text(size = 7),
        legend.text = element_text(size = 7), legend.key.width = grid::unit(4, "mm"),
        legend.key.height = grid::unit(3, "mm"), plot.caption = element_text(size = 6.5, hjust = 0),
        plot.tag = element_text(face = "bold", size = 10), plot.tag.position = "topleft",
        plot.margin = margin(5, 5, 5, 5))
    height_mm <- HEIGHT_MM
    prefix <- file.path(out, "Figure4B_candidate_pair_support")
    ggsave(paste0(prefix, ".pdf"), p, device = grDevices::cairo_pdf, width = WIDTH_MM, height = height_mm, units = "mm", bg = "white")
    ggsave(paste0(prefix, ".svg"), p, device = svglite::svglite, width = WIDTH_MM, height = height_mm, units = "mm", bg = "white")
    ggsave(paste0(prefix, ".png"), p, device = "png", type = "cairo", width = WIDTH_MM, height = height_mm, units = "mm", dpi = 600, bg = "white")
    ggsave(paste0(prefix, ".tiff"), p, device = "tiff", type = "cairo", compression = "lzw", width = WIDTH_MM, height = height_mm, units = "mm", dpi = 600, bg = "white")
    writeLines(c("width_mm=94; height_mm=98; font=Arial",
      "font_pt: panel=10; tissue_headers/cells/legend=7; gene_labels/prefix/caption=6.5",
      "placement=Figure4 right column; intended_scale=1",
      "assembly_and_caption_page_fit=NOT_CHECKED",
      paste("font_paths", paste(font_paths, collapse=";"))),
      file.path(out, "Figure4B_layout.txt"))
    write_tsv(plot_data, "Figure4B_plotted_cells.tsv")
    write_tsv(selected, "Figure4B_selected_pairs.tsv")
    write_tsv(data.table(parameter = c("manifest", "annotated_links", "annotation_tasks", "success_list",
      "pair_universe", "min_support", "window_bp", "bias_lower", "bias_upper", "max_display_pairs",
      "width_mm", "height_mm", "concurrent_tasks", "threads_per_task", "phase_comparability"),
      value = as.character(c(MANIFEST, LINKS, ANNOTATION_TASKS, SUCCESS, "ever_detected_candidate_pc_pairs",
        MIN_SUPPORT, WINDOW_BP, BIAS_LOWER, BIAS_UPPER, MAX_DISPLAY_PAIRS, WIDTH_MM, height_mm, 1, 1, "not_assessed"))),
      "Figure4B_parameters.tsv")
    writeLines(c(
      "Figure 4B. Cross-dog support for candidate–protein-coding gene links.",
      "Rows represent specific candidate–partner pairs previously detected in at least one dog-tissue sample. Columns are CB (cerebellum), FC (frontal cortex), OC (occipital cortex), and Pit (pituitary gland). This selected-pair display is not an estimate of discovery rates among all possible pairs.",
      "For each pair and sample, eligibility requires both specific endpoints to be present, each with site-level allele support >=10, and the protein-coding endpoint to intersect the candidate interval expanded by 500 kb on each side on the same autosome (original inclusive LINK comparisons). Bias is not required for this denominator. Stored rounded allelic ratios >=0.6 or <=0.4 define upstream bias flags; actual link calls are read from the annotated output, not inferred anew.",
      "Cell text gives detected/evaluable dogs. Fill encodes the number of detected dogs (0–5), not probability or significance. NA denotes no evaluable dogs; 0/n denotes evaluable but undetected. Missing samples are not counted as negative. Upstream task failures are not converted to biological non-detection.",
      "An asterisk is a review flag, not a significance symbol: at least one detected dog lacks a nonoverlap_no_multigene_flag link for that pair in that tissue. A dog with both flagged and unflagged records counts as unflagged-supported. Full stratum labels, unflagged dog counts and dog identities are available in the source tables. Lack of a flag is not functional validation.",
      "Candidate IDs identify the candidate catalogue, not necessarily independent novel genes. Bracketed annotation class codes are retained from upstream candidate annotations: [u], intergenic novel transcript; [p], possible polymerase run-on or readthrough; [x], antisense novel transcript. Multiple codes are separated by semicolons. Partner symbols are display labels only; gene IDs define pairs.",
      "The shared prefix DOG5_NEWGENE_ appears on a separate line above the column headers; for example, 03157 denotes DOG5_NEWGENE_03157. Leading zeros are retained. Bracketed classes follow candidate numbers, not partner names. Candidate and protein-coding partner labels occupy separate columns; full identifiers remain in the source tables.",
      paste0("Displayed ", nrow(selected), " of ", nrow(catalog), " previously detected pairs. Ranking uses descending maximum unflagged-support dogs within a tissue, maximum detected dogs within a tissue, distinct detected dogs across tissues, detected tissue count, and evaluable dog-tissue count, followed by gene IDs. This is a descriptive display ordering, not a high-confidence threshold. No GWAS/PPA filter is used; all pairs are retained in source tables."),
      "Support counts are summed site-level allele counts, not independent RNA molecules. Raw enhancing/repressing labels are retained for provenance only; phase comparability is not independently established here, and neither pair direction nor a causal regulatory mechanism is claimed. Dogs are counted once per pair within each tissue; observations across tissues from one dog are not independent animals.",
      "Visual reference: https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#cell_fun ; source: https://github.com/jokergoo/ComplexHeatmap-reference/blob/master/02-single_heatmap.Rmd . Independent ggplot2 implementation.",
      "Method reference: https://doi.org/10.1093/nar/gkaf912 (cross-individual support; thresholds are not copied from the human study)."
    ), file.path(out, "Figure4B_caption.txt"))
    capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
    log_stage("[5/5] main task completed normally, record execution status")
    writeLines(c("execution_status=success", paste0("output_dir=", out), paste0("master_log=", log_path)),
      file.path(out, "execution_completion.marker"))
    log_stage("task completed / done; execution_status=success; please check separately output completeness, source table and final dimensions layout.")
  }, error = function(e) {log_stage(paste("execution_status=failed;", conditionMessage(e))); stop(e)})
}

main()
