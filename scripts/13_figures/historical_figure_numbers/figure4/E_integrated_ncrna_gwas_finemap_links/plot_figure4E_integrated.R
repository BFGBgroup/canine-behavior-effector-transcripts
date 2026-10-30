#!/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/Rscript


Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
ROOT <- "/path/to/repository/scripts/13_figures/historical_figure_numbers/figure4/E_integrated_ncrna_gwas_finemap_links"
RAW <- "/path/to/project/analysis_20260817/results/5-dogs/4.allelome-compare-raw-gwas/ncrna-mrna/20260910_174940/dog5_annotated_ncrna_pcgene_raw_gwas_allelome_focus.tsv"
CS <- "/path/to/project/analysis_20260817/results/5-dogs/5.allelome-compare-gwas-finemap/ncrna-mrna/20260910_183054/dog5_annotated_ncrna_pcgene_gwas_finemap_focus.tsv"

WIDTH_MM <- 178
HEIGHT_MM <- 50
FONT <- "Arial"
OUT <- file.path(ROOT, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
if (!dir.create(OUT, recursive=TRUE, showWarnings=FALSE)) stop("cannot create output directory: ", OUT)
con <- file(file.path(OUT, "plot_figure4E_integrated.log"), "wt")
sink(con, split=TRUE)
sink(con, type="message")
log_stage <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep="")
main <- function() {
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(ggplot2))
  setDTthreads(1)
  font_paths <- systemfonts::match_fonts(c(FONT, FONT), weight=c("normal", "bold"))$path
  expected_fonts <- c("/path/to/fonts/Arial/Arial.TTF",
                      "/path/to/fonts/Arial/Arialbd.TTF")
  if (!all(file.exists(expected_fonts)) ||
      !identical(normalizePath(font_paths), normalizePath(expected_fonts)))
    stop("Arial regular or bold font is unavailable; please install the specified font, do not silently substitute.")
  log_stage("[1/4] read two records path; concurrent_tasks=1; threads_per_task=1")
  log_stage("RAW=", RAW, "\nCS=", CS, "\nOUT=", OUT)
  common <- c("ncrna_gene_id", "ncrna_gene_symbol", "pcgene_gene_id", "pcgene_gene_symbol",
              "trait", "gwas_chrom", "gwas_pos", "gwas_unordered_alleles", "gwas_hit_on_ncrna")
  read_input <- function(path, extra) {
    d <- fread(path, sep="\t", na.strings=c("", "NA"), nThread=1)
    fields <- c(common, extra)
    if (length(setdiff(fields, names(d)))) stop("missing fields: ", path)
    if (!nrow(d) || anyNA(d[, ..fields])) stop("input is empty or key fields are missing: ", path)
    d
  }
  a <- read_input(RAW, "gwas_p_lrt")
  b <- read_input(CS, c("gwas_branch", "gwas_analysis_scope", "candidate_type",
                       "gwas_in_credible_set", "gwas_locus_lead_p_lrt"))

  raw_input_rows <- nrow(a)
  cs_input_rows <- nrow(b)
  a <- a[gwas_hit_on_ncrna == "yes" & is.finite(gwas_p_lrt) &
           gwas_p_lrt >= 0 & gwas_p_lrt < 1e-5]
  b <- b[gwas_branch == "main" & gwas_analysis_scope == "primary" &
           candidate_type == "susie_credible_set" & gwas_in_credible_set == "yes" &
           gwas_hit_on_ncrna == "yes" & is.finite(gwas_locus_lead_p_lrt) &
           gwas_locus_lead_p_lrt >= 0 & gwas_locus_lead_p_lrt < 1e-5]
  log_stage("original GWAS filter: ", raw_input_rows, " -> ", nrow(a),
            "; main/primary credible set filter: ", cs_input_rows, " -> ", nrow(b))

  normalized <- function(d, route) {
    z <- copy(d[, ..common])
    z[, allele_set := vapply(strsplit(gwas_unordered_alleles, "/", fixed=TRUE),
                             function(x) paste(sort(unique(x)), collapse="/"), character(1))]
    z[, variant_id := paste(gwas_chrom, gwas_pos, allele_set, sep=":")]
    z[, route_name := route]
    z
  }
  d <- rbindlist(list(normalized(a, "GWAS"), normalized(b, "CS")))
  labels <- unique(d[, .(ncrna_gene_id, ncrna_gene_symbol)])
  if (anyDuplicated(labels$ncrna_gene_id)) stop("the same stable ID exists mismatch display name")
  classify <- function(raw, cs) fifelse(raw & cs, "Shared", fifelse(raw, "GWAS only", "CS only"))
  log_stage("[2/4] deduplicate matching; expanded rows RAW=", nrow(a), "; CS=", nrow(b))
  variants <- d[, .(in_gwas=any(route_name=="GWAS"), in_cs=any(route_name=="CS")),
                by=.(ncrna_gene_id, trait, variant_id, gwas_chrom, gwas_pos, allele_set)]
  variants[, category := classify(in_gwas, in_cs)]
  pairs <- d[, .(in_gwas=any(route_name=="GWAS"), in_cs=any(route_name=="CS"),
                  pcgene=paste(sort(unique(pcgene_gene_symbol)), collapse=";")),
             by=.(ncrna_gene_id, trait, pcgene_gene_id)]
  pairs[, category := classify(in_gwas, in_cs)]
  summary <- variants[, .(gwas_only=sum(category=="GWAS only"), shared=sum(category=="Shared"),
                          cs_only=sum(category=="CS only"), n_variants=.N), by=.(ncrna_gene_id, trait)]
  summary <- merge(summary, pairs[, .(n_partners=.N), by=.(ncrna_gene_id, trait)],
                   by=c("ncrna_gene_id", "trait"))
  summary <- merge(summary, labels, by="ncrna_gene_id")
  summary[, trait_label := gsub("_", " ", sub("^score_", "", trait))]
  summary[, trait_label := paste0(toupper(substr(trait_label, 1, 1)), substring(trait_label, 2))]
  summary[trait == "score_dog_directed_aggression", trait_label := "Dog-directed aggression"]
  summary[trait == "score_dog_directed_fear", trait_label := "Dog-directed fear"]
  summary[trait == "score_separation_related_problems", trait_label := "Separation-related problems"]

  setorder(summary, trait, ncrna_gene_id)
  width_mm <- WIDTH_MM
  height_mm <- HEIGHT_MM

  summary[, row_y := seq(height_mm-19.5, 4.25, length.out=.N)]
  categories <- c("GWAS only", "Shared", "CS only")

  partner_units <- merge(copy(pairs), summary[, .(ncrna_gene_id, trait, row_y)],
                         by=c("ncrna_gene_id", "trait"))
  setorder(partner_units, ncrna_gene_id, trait, pcgene_gene_id)
  partner_units[, unit_index := seq_len(.N), by=.(ncrna_gene_id, trait)]
  partner_units[, x := 95 + (unit_index-1)*3.5]
  variant_units <- merge(copy(variants), summary[, .(ncrna_gene_id, trait, row_y)],
                         by=c("ncrna_gene_id", "trait"))
  variant_units[, category_order := match(category, categories)]
  variant_units[, chromosome_order := suppressWarnings(as.integer(sub("^chr", "", gwas_chrom))) ]
  setorder(variant_units, ncrna_gene_id, trait, category_order, chromosome_order,
           gwas_chrom, gwas_pos, allele_set)
  variant_units[, unit_index := seq_len(.N), by=.(ncrna_gene_id, trait)]
  variant_units[, x := 126 + (unit_index-1)*4.2]
  write_table <- function(x, name) fwrite(x, file.path(OUT, name), sep="\t", na="NA", nThread=1)
  write_table(summary, "Figure4E_integrated_summary.tsv")
  write_table(variants, "Figure4E_variant_membership.tsv")
  write_table(pairs, "Figure4E_pair_membership.tsv")
  write_table(partner_units, "Figure4E_partner_unit_key.tsv")
  write_table(variant_units, "Figure4E_variant_unit_key.tsv")
  write_table(a, "Figure4E_raw_source_rows.tsv")
  write_table(b, "Figure4E_cs_source_rows.tsv")
  write_table(data.table(route=c("GWAS", "CS"), input=c(RAW, CS),
                        input_rows=c(raw_input_rows, cs_input_rows),
                        selected_rows=c(nrow(a), nrow(b))), "input_manifest.tsv")
  log_stage("lncRNA—trait row count=", nrow(summary), "; pair=", nrow(pairs),
            "; trait—variant union=", uniqueN(variants[, .(trait, variant_id)]))
  log_stage("[3/4] plot a single row of unit symbols: one circle per partner, one square per variant")


  cols <- c("GWAS only"="#D55E00", "Shared"="#0072B2", "CS only"="#009E73")
  canvas_width <- WIDTH_MM-6
  p <- ggplot() +
    geom_point(data=partner_units, aes(x=x, y=row_y), shape=16, size=1.65,
               colour="#24485C") +
    geom_point(data=variant_units, aes(x=x, y=row_y, colour=category),
               shape=15, size=2.5) +
    geom_text(data=summary, aes(x=0, y=row_y, label=ncrna_gene_symbol),
              hjust=0, size=7, size.unit="pt", family=FONT) +
    geom_text(data=summary, aes(x=30, y=row_y, label=trait_label),
              hjust=0, size=7, size.unit="pt", family=FONT) +
    annotate("text", x=0, y=height_mm-2.5, label="E", hjust=0,
             size=10, size.unit="pt", family=FONT, fontface="bold") +
    annotate("text", x=6, y=height_mm-2.5, label="GWAS-linked annotated lncRNAs",
             hjust=0, size=8, size.unit="pt", family=FONT, fontface="bold") +
    annotate("text", x=70, y=height_mm-8.5, label="Variant source:",
             hjust=0, size=7, size.unit="pt", family=FONT, colour="#555555") +
    annotate("point", x=c(94,122,148), y=height_mm-8.5,
             shape=15, size=1.9, colour=unname(cols)) +
    annotate("text", x=c(96.5,124.5,150.5), y=height_mm-8.5,
             label=categories, hjust=0, size=7, size.unit="pt", family=FONT) +
    annotate("text", x=c(0,30), y=height_mm-13,
             label=c("lncRNA", "Behavioral trait"), hjust=0,
             size=8, size.unit="pt", family=FONT, fontface="bold") +
    annotate("text", x=c(104,146), y=height_mm-13,
             label=c("Coding partners", "Overlapping variants"),
             size=8, size.unit="pt", family=FONT, fontface="bold") +
    annotate("text", x=c(104,146), y=height_mm-16,
             label=c("1 dot = 1 partner", "1 tile = 1 variant"),
             size=6.5, size.unit="pt", family=FONT, colour="#666666") +
    annotate("segment", x=0, xend=canvas_width, y=height_mm-18.2,
             yend=height_mm-18.2, colour="#AAB4BA", linewidth=.25) +
    scale_colour_manual(values=cols, limits=categories, guide="none") +
    scale_x_continuous(limits=c(0, canvas_width), expand=c(0,0)) +
    scale_y_continuous(limits=c(0, height_mm), expand=c(0,0)) +
    theme_void(base_size=8, base_family=FONT) +
    theme(legend.position="none", plot.margin=margin(0,3,0,3,unit="mm"))
  stem <- file.path(OUT, "Figure4E_integrated")
  grDevices::cairo_pdf(paste0(stem,".pdf"), width=width_mm/25.4, height=height_mm/25.4, family=FONT)
  print(p); grDevices::dev.off()
  grDevices::png(paste0(stem,".png"), width=width_mm, height=height_mm, units="mm", res=600, type="cairo", bg="white")
  print(p); grDevices::dev.off()
  grDevices::tiff(paste0(stem,".tiff"), width=width_mm, height=height_mm, units="mm", res=600, type="cairo", compression="lzw", bg="white")
  print(p); grDevices::dev.off()
  if (requireNamespace("svglite", quietly=TRUE)) {
    svglite::svglite(paste0(stem,".svg"), width=width_mm/25.4, height=height_mm/25.4)
    print(p); grDevices::dev.off()
  } else log_stage("SVG skipped: optional package svglite unavailable")
  writeLines(c("width_mm=178; height_mm=50; font=Arial",
    "font_pt: panel=10; panel_title/column_headers=8; row_labels/legend=7; unit_explanations=6.5",
    "placement=Figure4 bottom full-width panel; intended_scale=1",
    "assembly_and_caption_page_fit=NOT_CHECKED",
    paste("font_paths",paste(font_paths,collapse=";"))),
    file.path(OUT,"Figure4E_layout.txt"))
  writeLines(c(
    "Annotated lncRNAs linked to behavioral GWAS through two overlapping selection routes.",
    sprintf("The panel includes %d annotated lncRNAs across %d behavioral traits, %d partner links and %d distinct trait-variant combinations (%d GWAS only, %d Shared and %d CS only row-level variant memberships).",
      uniqueN(summary$ncrna_gene_id), uniqueN(summary$trait), nrow(pairs),
      uniqueN(variants[, .(trait, variant_id)]), sum(summary$gwas_only),
      sum(summary$shared), sum(summary$cs_only)),
    "Each row is an lncRNA-trait combination. Each equally sized dark circle represents one distinct coding partner in the union of routes.",
    "Each equally sized colored square represents one distinct overlapping variant, classified as GWAS only, Shared or CS only. Shared variants are drawn once.",
    "Rows are ordered by trait and lncRNA ID, not evidence strength. Variants are arranged by route category and then chromosome, position and allele set; partners by stable gene ID.",
    "Symbol positions are counting slots, not genomic coordinates or shared identities across rows. Partner circles and variant squares are not paired one-to-one. Unit-key tables identify every plotted member.",
    "Colors follow the existing Figure 4 project palette: GWAS only #D55E00, Shared #0072B2, CS only #009E73. Colors encode route membership, not evidence strength.",
    "GWAS route: variant P_LRT < 1e-5 and lncRNA gene-span overlap. CS route: main branch, primary region lead P_LRT < 1e-5, credible-set membership and lncRNA gene-span overlap; no member P cutoff.",
    "Variants are matched by lncRNA ID, trait, chromosome, position and unordered allele set. Coding partners are deduplicated by stable gene ID within each row.",
    "Shared denotes membership in both selected results, not independent replication. Route-specific absence does not prove lack of association or exclusion by the fine-mapping model.",
    "Variant counts do not indicate independent signals, full credible-set sizes or causal probabilities. LINK is not demonstrated causal regulation.",
    sprintf("Figure dimensions: %.1f x %.1f mm.", width_mm, height_mm)), file.path(OUT,"Figure4E_caption.txt"))
  writeLines(capture.output(sessionInfo()), file.path(OUT,"sessionInfo.txt"))
  log_stage("[4/4] task completed / done; execution_status=success")
  writeLines("execution_status=success", file.path(OUT,"execution_completion.marker"))
  log_stage("please check separately source table statistics, PDF layout and biological interpretation; completion marker only indicates execution success.")
}
status <- tryCatch({main(); 0L}, error=function(e) {
  log_stage("execution_status=failed; ", conditionMessage(e)); 1L
})
sink(type="message"); sink(); close(con)
quit(save="no", status=status)
