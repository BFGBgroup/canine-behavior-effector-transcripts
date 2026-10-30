#!/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/Rscript


Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
ROOT <- "/path/to/project/analysis_20260817"
OUT_ROOT <- file.path(ROOT, "article/figure/figure4/A_candidate_link_coverage")
LINK_ROOT <- file.path(ROOT, "5-dogs/7.Allelome.LINK")
ENDPOINT_RUN <- file.path(LINK_ROOT, "1.2_build_dog5_allelome_locus_tables/20260907_202637")
MANIFEST <- file.path(ENDPOINT_RUN, "dog5_allelome_locus_table_manifest.tsv")
ANNOTATED <- file.path(LINK_ROOT, "1.4_annotate_and_merge_dog5_allelome_links/20260907_213429/dog5_allelome_links.annotated.tsv")
SUCCESS <- file.path(LINK_ROOT, "scripts/1.3_run_dog5_allelome_link/log/1.3_run_dog5_allelome_link.20260907_205139.success.tsv")
MIN_SUPPORT <- 10
WINDOW_BP <- 500000
BIAS_UPPER <- 0.6
BIAS_LOWER <- 1 - BIAS_UPPER
WIDTH_MM <- 80

HEIGHT_MM <- 40
FONT <- "Arial"

main <- function() {

  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out <- file.path(OUT_ROOT, stamp)
  if (!dir.create(out, recursive = TRUE, showWarnings = FALSE)) stop("output directory already exists or cannot create: ", out)
  log_dir <- file.path(OUT_ROOT, "log")
  dir.create(log_dir, showWarnings = FALSE)
  log_path <- file.path(log_dir, paste0("plot_figure4A.", stamp, ".", Sys.getpid(), ".log"))
  con <- file(log_path, open = "wt")
  sink(con, type = "output")
  sink(con, type = "message")
  on.exit({sink(type = "message"); sink(type = "output"); close(con)}, add = TRUE)
  log_stage <- function(x) {cat(format(Sys.time(), "%H:%M:%S"), x, "\n"); flush(con)}
  tryCatch({
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(ggplot2))
    setDTthreads(1)

    font_paths <- systemfonts::match_fonts(c(FONT,FONT),weight=c("normal","bold"))$path
    expected_fonts <- c("/path/to/fonts/Arial/Arial.TTF",
                        "/path/to/fonts/Arial/Arialbd.TTF")
    if (!all(file.exists(expected_fonts)) ||
        !identical(normalizePath(font_paths),normalizePath(expected_fonts)))
      stop("Arial regular or bold font was not resolved correctly, please check font installation")
    log_stage("[1/5] read actual input; serial 1 task, per-task 1 thread")
    cat("MANIFEST=", MANIFEST, "\nANNOTATED=", ANNOTATED, "\nSUCCESS=", SUCCESS, "\nOUTPUT=", out, "\n", sep = "")
    read_tsv <- function(path) fread(file = path, sep = "\t", na.strings = "NA", showProgress = FALSE)
    write_tsv <- function(x, name) fwrite(x, file.path(out, name), sep = "\t", na = "NA")
    manifest <- read_tsv(MANIFEST)
    manifest <- manifest[chromosome %in% paste0("chr", 1:38)]
    status <- read_tsv(SUCCESS)

    si <- match(manifest$task_id, status$task_id)
    if (anyNA(si) || any(status$execution_status[si] != "success") ||
        any(status$input_file[si] != manifest$locus_table)) {
      stop("Endpoint manifest does not match the completed LINK tasks; coverage cannot be calculated.")
    }
    links <- read_tsv(ANNOTATED)
    links <- links[pair_type == "candidate_new_gene--known_protein_coding"]
    links[, candidate_gene_id := fifelse(base_endpoint_type == "candidate_new_gene", base_gene_id, target_gene_id)]
    links[, partner_gene_id := fifelse(base_endpoint_type == "candidate_new_gene", target_gene_id, base_gene_id)]
    link_counts <- links[, .(
      observed_link_rows = .N,
      linked_pc_partners = uniqueN(partner_gene_id),
      noflag_link_rows = sum(evidence_stratum == "nonoverlap_no_multigene_flag"),
      noflag_pc_partners = uniqueN(partner_gene_id[evidence_stratum == "nonoverlap_no_multigene_flag"])
    ), by = .(animal_id, tissue, gene_id = candidate_gene_id)]

    log_stage("[2/5] from endpoint reconstruct evaluable candidate; do not filter the denominator by positive allelic bias")
    loci <- vector("list", nrow(manifest))
    for (i in seq_len(nrow(manifest))) {
      task <- manifest[i]
      endpoint <- read_tsv(task$locus_table)
      meta <- read_tsv(task$locus_metadata)

      mi <- match(endpoint$name, meta$group_id)
      if (anyNA(mi) || anyDuplicated(meta$group_id) || anyDuplicated(endpoint$name)) {
        stop("endpoint and metadata cannot one-to-one matching: ", task$task_id)
      }
      endpoint[, `:=`(gene_id = meta$gene_id[mi], record_type = meta$record_type[mi],
                      gene_biotype = meta$gene_biotype[mi])]
      candidate <- endpoint[record_type == "manual_new_gene"]

      pc <- endpoint[!record_type %in% c("manual_new_gene", "unannotated_variant") &
                       gene_biotype == "protein_coding" & total_reads >= MIN_SUPPORT]


      partner_n <- vapply(seq_len(nrow(candidate)), function(j) {
        uniqueN(pc$gene_id[pc$chr == candidate$chr[j] & pc$name != candidate$name[j] &
                          pc$end >= candidate$start[j] - WINDOW_BP &
                          pc$start <= candidate$end[j] + WINDOW_BP])
      }, integer(1))
      ci <- match(candidate$name, meta$group_id)
      loci[[i]] <- candidate[, .(gene_id, endpoint_name = name, chromosome = chr,
        start0 = start, end0 = end, hap1_support = A1_reads, hap2_support = A2_reads,
        total_support = total_reads, allelic_ratio_used = allelic_ratio)]
      loci[[i]][, `:=`(animal_id = task$animal_id, tissue = task$tissue, task_id = task$task_id,
        eligible_pc_partners = partner_n, support_pass = total_support >= MIN_SUPPORT,
        bias_pass = allelic_ratio_used >= BIAS_UPPER | allelic_ratio_used <= BIAS_LOWER,
        allelic_ratio_raw = meta$allelic_ratio_raw[ci], variant_count = meta$variant_count[ci],
        multi_gene_variant_count = meta$multi_gene_variant_count[ci],
        count_unit = meta$count_unit[ci], source_locus_table = task$locus_table,
        source_metadata = task$locus_metadata, source_gene_table = meta$input_file[ci],
        source_gene_line = meta$input_line[ci])]
      if (i == 1L || i == nrow(manifest) || task$sample_key != manifest$sample_key[i + 1L]) {
        log_stage(paste("already process endpoint task", i, "/", nrow(manifest), task$sample_key))
      }
    }
    loci <- rbindlist(loci)
    loci[, eligible := support_pass & eligible_pc_partners > 0L]
    loci <- merge(loci, link_counts, by = c("animal_id", "tissue", "gene_id"), all.x = TRUE, sort = FALSE)
    count_columns <- c("observed_link_rows", "linked_pc_partners", "noflag_link_rows", "noflag_pc_partners")
    for (nm in count_columns) set(loci, which(is.na(loci[[nm]])), nm, 0L)
    loci[, `:=`(linked = observed_link_rows > 0L, linked_noflag = noflag_link_rows > 0L)]

    log_stage("[3/5] by candidate gene_id deduplicate summary; save denominator, bias and link level")
    summary <- loci[, .(
      N_candidate_observed = uniqueN(gene_id),
      N_eligible = uniqueN(gene_id[eligible]),
      N_biased = uniqueN(gene_id[eligible & bias_pass]),
      N_linked = uniqueN(gene_id[eligible & linked]),
      N_linked_noflag = uniqueN(gene_id[eligible & linked_noflag])
    ), by = .(animal_id, tissue)]
    samples <- unique(manifest[, .(animal_id, tissue)])
    summary <- merge(samples, summary, by = c("animal_id", "tissue"), all.x = TRUE)
    for (nm in grep("^N_", names(summary), value = TRUE)) set(summary, which(is.na(summary[[nm]])), nm, 0L)
    summary[, `:=`(
      linked_pct = fifelse(N_eligible > 0, 100 * N_linked / N_eligible, NA_real_),
      linked_noflag_pct = fifelse(N_eligible > 0, 100 * N_linked_noflag / N_eligible, NA_real_),
      biased_pct = fifelse(N_eligible > 0, 100 * N_biased / N_eligible, NA_real_),
      linked_among_biased_pct = fifelse(N_biased > 0, 100 * N_linked / N_biased, NA_real_))]
    dogs <- c("Ace_396226", "CodRou_346627", "Crak_358336", "Dozer_386519", "SanMcg_349891")
    dog_labels <- c("Ace", "CodRou", "Crak", "Dozer", "SanMcg")
    tissues <- c("cerebellum", "cortex_frontal", "cortex_occipital", "pituitary_gland")
    summary[, `:=`(tissue_index = match(tissue, tissues), dog_index = match(animal_id, dogs))]
    summary[, plot_x := tissue_index + c(-0.22, -0.11, 0, 0.11, 0.22)[dog_index]]
    setorder(summary, tissue_index, dog_index)
    medians <- summary[, .(n_dogs = sum(!is.na(linked_pct)),
      median_linked_pct = if (all(is.na(linked_pct))) NA_real_ else median(linked_pct, na.rm = TRUE)),
      by = .(tissue, tissue_index)]
    write_tsv(summary, "Figure4A_sample_summary.tsv")
    write_tsv(medians, "Figure4A_tissue_medians.tsv")
    write_tsv(loci, "Figure4A_candidate_locus_source.tsv")
    write_tsv(links, "Figure4A_observed_candidate_pc_links.tsv")
    write_tsv(manifest, "Figure4A_input_manifest.tsv")
    write_tsv(data.table(parameter = c("endpoint_manifest", "annotated_links", "success_list", "min_support",
      "window_bp", "bias_lower", "bias_upper", "ratio_source", "denominator", "threads", "width_mm", "height_mm"),
      value = as.character(c(MANIFEST, ANNOTATED, SUCCESS, MIN_SUPPORT, WINDOW_BP, BIAS_LOWER, BIAS_UPPER,
        "stored_rounded_allelic_ratio", "supported_candidate_with_supported_local_pc_endpoint", 1, WIDTH_MM, HEIGHT_MM))),
      "Figure4A_parameters.tsv")

    log_stage("[4/5] plot per-dog points and within-tissue medians, export vector and 600 dpi figure")
    palette <- setNames(c("#D55E00", "#0072B2", "#009E73", "#E69F00", "#CC79A7"), dogs)
    shapes <- setNames(c(21, 22, 24, 23, 25), dogs)

    tissue_n <- medians$n_dogs[match(seq_along(tissues), medians$tissue_index)]
    tissue_labels <- paste0(c("CB", "FC", "OC", "Pit"), "\nn = ", tissue_n)
    p <- ggplot(summary, aes(x = plot_x, y = linked_pct)) +
      geom_segment(data = medians[!is.na(median_linked_pct)],
        aes(x = tissue_index - 0.28, xend = tissue_index + 0.28,
            y = median_linked_pct, yend = median_linked_pct),
        inherit.aes = FALSE, linewidth = 0.32, colour = "black") +
      geom_point(aes(colour = animal_id, fill = animal_id, shape = animal_id),
                 size = 1.6, stroke = 0.25, na.rm = TRUE) +
      scale_colour_manual(values = palette, breaks = dogs, labels = dog_labels, name = NULL) +
      scale_fill_manual(values = palette, breaks = dogs, labels = dog_labels, name = NULL) +
      scale_shape_manual(values = shapes, breaks = dogs, labels = dog_labels, name = NULL) +
      scale_x_continuous(breaks = 1:4, labels = tissue_labels, limits = c(0.5, 4.5), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.12))) +

      labs(x = NULL, y = "Eligible candidates\nlinked (%)", tag = "A") +
      guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1),
             shape = guide_legend(nrow = 1)) +
      theme_classic(base_size = 7, base_family = FONT) +
      theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        axis.text = element_text(colour = "black", size = 7), axis.title.y = element_text(size = 7.5, margin = margin(r = 4)),
        axis.text.x = element_text(lineheight = 0.9),
        legend.position = "top", legend.justification = "left", legend.text = element_text(size = 7),
        legend.key.width = grid::unit(3, "mm"), legend.key.height = grid::unit(3, "mm"),
        legend.spacing.x = grid::unit(1, "mm"), legend.margin = margin(0, 0, 2, 0),
        plot.tag = element_text(size = 10, face = "bold"), plot.tag.position = "topleft",
        plot.margin = margin(5, 6, 5, 5))
    prefix <- file.path(out, "Figure4A_candidate_link_coverage")
    ggsave(paste0(prefix, ".pdf"), p, device = grDevices::cairo_pdf, width = WIDTH_MM, height = HEIGHT_MM, units = "mm", bg = "white")
    ggsave(paste0(prefix, ".svg"), p, device = svglite::svglite, width = WIDTH_MM, height = HEIGHT_MM, units = "mm", bg = "white")
    ggsave(paste0(prefix, ".png"), p, device = "png", type = "cairo", width = WIDTH_MM, height = HEIGHT_MM, units = "mm", dpi = 600, bg = "white")
    ggsave(paste0(prefix, ".tiff"), p, device = "tiff", type = "cairo", compression = "lzw", width = WIDTH_MM, height = HEIGHT_MM, units = "mm", dpi = 600, bg = "white")
    writeLines(c("width_mm=80; height_mm=40; font=Arial",
      "font_pt: panel=10; axis_title=7.5; ticks/legend=7",
      "legend=one_row; y_title=two_lines; point_size_mm=1.6",
      "x_labels=tissue_and_observed_dog_n; n_source=nonmissing_linked_pct",
      "Figure4 planned left column=80 mm; right column B=94 mm; gutter=4 mm",
      "PDF/SVG vector; PNG/TIFF 600 dpi; TIFF LZW; white background",
      "final_size_collision_and_assembly_QA=NOT_CHECKED",
      paste("font_paths",paste(font_paths,collapse=";"))),
      file.path(out,"Figure4A.layout_parameters.txt"))
    writeLines(c(
      "Figure 4A. Detection of protein-coding-gene links among eligible candidate loci.",
      "Each point represents one dog-tissue sample; colours, shapes and fixed horizontal offsets identify dogs. Black lines show tissue-wise medians of individual percentages, not pooled locus proportions. The n below each tissue is the number of dogs with a nonmissing plotted percentage, not the number of candidate loci. CB, cerebellum; FC, frontal cortex; OC, occipital cortex; Pit, pituitary gland.",
      "Eligible loci are candidate_new_gene endpoints with total support >=10 and at least one known protein-coding endpoint with support >=10 on the same autosome within the original LINK interval window (500 kb on either side of the candidate interval, including boundary contact). Bias is not required for eligibility.",
      "The plotted percentage is 100 x N_linked / N_eligible. Each candidate is counted once per dog-tissue sample regardless of the number of partners. Eligibility is determined separately in each sample; similar percentages do not imply identical candidate sets or partner identities. This is a conditional detection proportion, not the proportion of all discovered candidates with validated functions.",
      "The numerator counts distinct eligible candidate gene IDs with at least one observed candidate--protein-coding link in the annotated LINK output, considering both endpoint roles. Overlap and multigene-review strata remain included in the main figure. A noflag sensitivity percentage uses the same denominator and requires at least one nonoverlap_no_multigene_flag link; other flagged links do not disqualify that candidate.",
      "Bias categories use the stored two-decimal allelic ratio (>=0.6 or <=0.4), matching the upstream procedure. Support is the sum of site-level allele-support counts, not a deduplicated count of independent RNA molecules. Links are descriptive algorithmic associations, not proof of regulation or an independent significance test.",
      paste0("Numbers of dogs with a plotted percentage by tissue: ",
        paste(c("CB", "FC", "OC", "Pit")[medians$tissue_index], medians$n_dogs, sep = "=", collapse = "; "), "."),
      "Missing samples are not imputed. A zero denominator gives NA. Observations across tissues from the same dog are repeated measurements, not independent animals. No GWAS or molecular PPA filter is applied, and no tissue-difference significance test is shown.",
      "The main percentage includes overlap and multigene-review links. The noflag percentage in the sample source table is a sensitivity summary, not a validated-regulation rate. Contact with a protein-coding partner under the LINK rules does not establish a direct regulatory mechanism.",
      "Layout: 80 x 40 mm; Arial; vector PDF/SVG and 600 dpi PNG/LZW TIFF. Final assembly and text legibility require separate rendered inspection.",
      "References: https://doi.org/10.1093/nar/gkaf912 ; https://doi.org/10.1371/journal.pbio.3000202 ; https://doi.org/10.1371/journal.pbio.1002128",
      "Style-only reference: https://github.com/JoachimGoedhart/PlotsOfData/blob/master/app.R ; the figure is independently implemented and does not reproduce the original paper's denominator."
    ), file.path(out, "Figure4A_caption.txt"))
    capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
    log_stage("[5/5] write execution completion marker; automated auditing of output content was not performed")
    writeLines(c("execution_status=success", paste0("output_dir=", out), paste0("master_log=", log_path)),
               file.path(out, "execution_completion.marker"))
    log_stage("task completed / done; execution_status=success; please check separately output completeness, font layout and biological interpretation.")
  }, error = function(e) {
    log_stage(paste("execution_status=failed;", conditionMessage(e)))
    stop(e)
  })
}

main()
