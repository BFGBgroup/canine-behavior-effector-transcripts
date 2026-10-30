#!/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/Rscript


Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
           NUMEXPR_NUM_THREADS = "1")

if (!nzchar(Sys.setlocale("LC_CTYPE", "C.utf8")))
  stop("Cannot enable C.utf8; UTF-8 input and log output require this locale.")
options(stringsAsFactors = FALSE, warn = 1)
library(grid)

INPUT <- "/path/to/project/analysis_20260817/results/5-dogs/8.new-gene-allelome-eqtl-gwas-finemap-integration/1.1_integrate_dog5_newgene_allelome_eqtl_gwas_finemap/20260910_195332/dog5_newgene_eqtl_gwas_finemap_all_candidates.tsv"
PARTNER_INPUT <- "/path/to/project/analysis_20260817/results/5-dogs/8.new-gene-allelome-eqtl-gwas-finemap-integration/1.1_integrate_dog5_newgene_allelome_eqtl_gwas_finemap/20260910_195332/dog5_candidate_newgene_pcGene_mechanism_pairs.tsv"
OUTPUT_ROOT <- "/path/to/repository/scripts/13_figures/historical_figure_numbers/figure5/A_integrated_candidate_evidence"
PRODUCER <- "/path/to/repository/scripts/12_evidence_integration/8.new-gene-allelome-eqtl-gwas-finemap-integration/scripts/1.1_integrate_dog5_newgene_allelome_eqtl_gwas_finemap/1.1.integrate_dog5_newgene_allelome_eqtl_gwas_finemap.py"
POSITION_INPUT <- "/path/to/project/analysis_20260817/results/5-dogs/5.allelome-compare-gwas-finemap/summary/20260910_182124/5dogs_gwas_finemap_allelome_main_focus.tsv"
RAW_INPUT <- "/path/to/project/analysis_20260817/results/5-dogs/7.new-gene-allelome-eqtl-finemap-raw-gwas-integration/summary/results/20260917_202456/dog5_candidate_blocks.tsv"
LOCUS_INPUT <- "/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.10_run_dog_behavior_susie_rss_finemapping/20260802_205352/locus_run_manifest.tsv"
POSITION_PRODUCER <- "/path/to/repository/scripts/12_evidence_integration/5.allelome-compare-gwas-finemap/summary/scripts/build_dog5_allelome_gwas_finemap_main_summary.py"
RAW_PRODUCER <- "/path/to/repository/scripts/12_evidence_integration/7.new-gene-allelome-eqtl-finemap-raw-gwas-integration/summary/scripts/1.1_build_dog5_newgene_allelome_eqtl_raw_gwas_summary.py"
SCRIPT <- file.path(OUTPUT_ROOT, "plot_figure5A.R")
WIDTH_MM <- 178
LEAD_P_THRESHOLD <- 1e-5
FONT <- "Arial"
FONT_PT <- 7
BLUE <- "#0072B2"
TEAL <- "#00856A"

log_stage <- function(text) {
  cat(sprintf("[%s ] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text))
  flush.console()
}


format_probability <- function(x, scientific_small = FALSE) {
  vapply(x, function(value) {
    if (is.na(value)) return("-")
    if (!scientific_small || value == 0 || value >= 0.001) return(sprintf("%.3f", value))

    sprintf("%.2e", value)
  }, character(1))
}

prepare_records <- function() {
  d <- read.delim(INPUT, sep = "\t", header = TRUE, quote = "", comment.char = "",
                  colClasses = "character", check.names = FALSE, na.strings = c("NA", ""),
                  encoding = "UTF-8", fileEncoding = "UTF-8")
  required <- c("candidate_id", "new_gene_id", "eqtl_tissue", "gwas_trait",
                "gwas_locus_id", "gwas_overlap_snp_id", "gwas_analysis_scope",
                "gwas_overlap_snp_in_credible_set", "allelome_block_match",
                "pcgene_link_count", "supporting_dog_count", "eqtl_ppa",
                "gwas_overlap_snp_pip", "evidence_tier", "eqtl_is_credible_set_member",
                "eqtl_priority_ppa_ge_0.1", "gwas_branch",
                "gwas_locus_lead_p_lrt", "gwas_overlap_snp_p_lrt")

  missing <- setdiff(required, names(d))
  if (length(missing)) stop("input missing plot field: ", paste(missing, collapse = ", "))
  if (!nrow(d) || anyNA(d$candidate_id) || anyDuplicated(d$candidate_id))
    stop("input is empty or candidate_id missing/duplicate; cannot create record level mapping.")
  tiers <- c("A_BOTH_95PCT_CREDIBLE_PLUS_ALLELOME",
             "B_GWAS_NONCREDIBLE_PLUS_ALLELOME",
             "C_BOTH_95PCT_CREDIBLE_NO_ALLELOME",
             "D_GWAS_NONCREDIBLE_NO_ALLELOME")
  tissue <- c(cerebellum = "CB", cortex_frontal = "FC", cortex_occipital = "OC",
              pituitary_gland = "Pit")
  scope <- c(primary = "PRI", exploratory = "EXP")
  if (anyNA(d$evidence_tier) || any(!d$evidence_tier %in% tiers) ||
      anyNA(d$eqtl_tissue) || any(!d$eqtl_tissue %in% names(tissue)) ||
      anyNA(d$gwas_analysis_scope) || any(!d$gwas_analysis_scope %in% names(scope)))
    stop("undefined values encountered tier, tissue or GWAS region span.")
  if (anyNA(d$gwas_branch) || any(d$gwas_branch != "main") ||
      anyNA(d$eqtl_is_credible_set_member) || any(d$eqtl_is_credible_set_member != "1") ||
      anyNA(d$eqtl_priority_ppa_ge_0.1) || any(d$eqtl_priority_ppa_ge_0.1 != "1"))
    stop("current script figure legend requirement main branch, molecular CS member and PPA>=0.1; input conditions have changed.")

  input_count <- nrow(d)
  d$source_row <- seq_len(nrow(d))
  lead_p <- suppressWarnings(as.numeric(d$gwas_locus_lead_p_lrt))
  ppa_value <- suppressWarnings(as.numeric(d$eqtl_ppa))
  valid <- !is.na(lead_p) & is.finite(lead_p) & lead_p > 0 & lead_p <= 1
  rules <- data.frame(
    new_gene = !is.na(d$new_gene_id) & grepl("^DOG5_NEWGENE_[0-9]+$", d$new_gene_id),
    regional_lead_p_lt_1e5 = valid & lead_p < LEAD_P_THRESHOLD,
    molecular_cs_ppa_ge_0.1 = d$eqtl_is_credible_set_member %in% "1" &
      !is.na(ppa_value) & is.finite(ppa_value) & ppa_value >= 0.1 & ppa_value <= 1,
    gwas_cs = d$gwas_overlap_snp_in_credible_set %in% "yes",
    matched_link = d$allelome_block_match %in% "yes")
  keep <- rep(TRUE, nrow(d))
  d$exclusion_reason <- ""
  for (rule in names(rules)) {
    passed <- rules[[rule]]
    d[[paste0("passes_", rule)]] <- passed
    keep <- keep & passed
    failed <- which(!passed)
    d$exclusion_reason[failed] <- ifelse(nzchar(d$exclusion_reason[failed]),
      paste(d$exclusion_reason[failed], paste0("fails_", rule), sep = ";"), paste0("fails_", rule))
  }
  excluded <- d[!keep, , drop = FALSE]
  d <- d[keep, , drop = FALSE]
  if (!nrow(d)) stop("No records remain after the combined dual-CS and LINK filter for new genes; cannot plot the evidence figure.")
  d <- d[order(match(d$evidence_tier, tiers), d$new_gene_id,
               match(d$eqtl_tissue, names(tissue)), d$gwas_trait,
               d$gwas_locus_id, d$gwas_overlap_snp_id, d$candidate_id), , drop = FALSE]

  d$plot_column <- seq_len(nrow(d))
  d$record_id <- sprintf("R%02d", d$plot_column)
  d$display_candidate <- sub("^DOG5_NEWGENE_", "", d$new_gene_id)
  d$display_tissue <- unname(tissue[d$eqtl_tissue])
  d$display_scope <- unname(scope[d$gwas_analysis_scope])


  d$display_record <- paste(d$display_candidate, d$display_tissue, sep = " / ")
  d$display_trait <- gsub("_", " ", sub("^score_", "", d$gwas_trait))
  substr(d$display_trait, 1, 1) <- toupper(substr(d$display_trait, 1, 1))
  short_trait <- c("Dog directed aggression" = "Dog-dir. aggression",
                   "Stranger directed fear" = "Stranger-dir. fear",
                   "Dog directed fear" = "Dog-dir. fear",
                   "Separation related problems" = "Separation-related problems")
  matched_trait <- match(d$display_trait, names(short_trait))
  d$display_trait[!is.na(matched_trait)] <- unname(short_trait[matched_trait[!is.na(matched_trait)]])
  d$display_variant <- paste0(sub("^Affx-", "", d$gwas_overlap_snp_id), " / L",
                              as.integer(sub(".*_", "", d$gwas_locus_id)))
  d$display_group <- c("GWAS CS + LINK", "Non-CS + LINK", "GWAS CS, no matched LINK",
                        "Non-CS, no matched LINK")[match(d$evidence_tier, tiers)]
  for (field in c("gwas_overlap_snp_in_credible_set", "allelome_block_match")) {
    if (any(!is.na(d[[field]]) & !d[[field]] %in% c("yes", "no")))
      stop("unknown status encoding: ", field)
  }
  d$plot_gwas_cs <- match(d$gwas_overlap_snp_in_credible_set, c("no", "yes")) - 1L
  d$plot_any_link <- match(d$allelome_block_match, c("no", "yes")) - 1L
  pc <- as.numeric(d$pcgene_link_count)
  ppa <- as.numeric(d$eqtl_ppa)
  pip <- as.numeric(d$gwas_overlap_snp_pip)
  if (any(!is.na(ppa) & (!is.finite(ppa) | ppa < 0 | ppa > 1)) ||
      any(!is.na(pip) & (!is.finite(pip) | pip < 0 | pip > 1)))
    stop("posterior numeric value cannot interpretation is 0 to 1 between of probability.")


  partners <- read.delim(PARTNER_INPUT, sep = "\t", quote = "", comment.char = "",
                         colClasses = "character", check.names = FALSE,
                         na.strings = c("NA", ""), fileEncoding = "UTF-8")
  pairs <- unique(partners[, c("candidate_id", "partner_gene_id"), drop = FALSE])
  pairs <- pairs[!is.na(pairs$candidate_id) & !is.na(pairs$partner_gene_id), , drop = FALSE]
  counts <- table(pairs$candidate_id)
  d$coding_partner_count <- as.integer(counts[match(d$candidate_id, names(counts))])

  d$coding_partner_count[is.na(d$coding_partner_count) & !is.na(pc) & pc == 0] <- 0L
  d$display_coding_partners <- ifelse(is.na(d$coding_partner_count), "-", as.character(d$coding_partner_count))

  d$pc_link_review <- vapply(d$candidate_id, function(id) {
    z <- partners[which(partners$candidate_id == id & !is.na(partners$partner_gene_id)), , drop = FALSE]
    any(grepl("multigene_review", z$allelome_evidence_strata), na.rm = TRUE)
  }, logical(1))
  d$plot_pc_link <- as.integer(pc > 0)
  d$display_ppa <- format_probability(ppa)
  d$display_pip <- format_probability(pip, TRUE)

  d$display_lead_p <- sprintf("%.2e", as.numeric(d$gwas_locus_lead_p_lrt))
  d$display_variant_p <- sprintf("%.2e", as.numeric(d$gwas_overlap_snp_p_lrt))
  d$display_lead_p[is.na(d$gwas_locus_lead_p_lrt)] <- "-"
  d$display_variant_p[is.na(d$gwas_overlap_snp_p_lrt)] <- "-"


  d$coding_link_dogs <- vapply(d$candidate_id, function(id) {
    z <- partners[which(partners$candidate_id == id & !is.na(partners$partner_gene_id)), , drop = FALSE]
    values <- z$supporting_dogs[!is.na(z$supporting_dogs)]
    dogs <- unique(unlist(strsplit(values, "|", fixed = TRUE), use.names = FALSE))
    dogs <- dogs[!is.na(dogs) & nzchar(dogs) & dogs != "NA"]
    if (!length(dogs)) return(NA_character_)
    paste(sort(dogs), collapse = "|")
  }, character(1))
  d$coding_link_dog_count <- vapply(d$coding_link_dogs, function(x) {
    if (is.na(x)) return(NA_integer_)
    length(strsplit(x, "|", fixed = TRUE)[[1]])
  }, integer(1))
  d$coding_link_dog_count[!is.na(d$coding_partner_count) & d$coding_partner_count == 0L] <- 0L
  d$display_dogs <- ifelse(is.na(d$coding_link_dog_count), "-", as.character(d$coding_link_dog_count))

  d$plot_link_class <- "unknown"
  d$plot_link_class[which(d$plot_any_link == 0L & d$plot_pc_link == 0L)] <- "none"
  d$plot_link_class[which(d$plot_any_link == 1L & d$plot_pc_link == 0L)] <- "other"
  d$plot_link_class[which(d$plot_any_link == 1L & d$plot_pc_link == 1L)] <- "pc"
  d$plot_ppa <- ppa
  d$plot_pip <- pip
  list(records = d, excluded = excluded, input_count = input_count)
}


read_source <- function(path) {
  read.delim(path, sep = "\t", quote = "", comment.char = "", colClasses = "character",
             check.names = FALSE, na.strings = c("NA", ""), fileEncoding = "UTF-8")
}
collapse_values <- function(x) paste(sort(unique(x[!is.na(x) & nzchar(x)])), collapse = ";")

variant_identity <- function(x) {
  vapply(strsplit(gsub(" ", ":", trimws(x), fixed = TRUE), ":", fixed = TRUE), function(z) {
    if (length(z) != 3L) stop("cannot interpretation variant identity: ", paste(z, collapse = ":"))
    a <- strsplit(z[3], "/", fixed = TRUE)[[1]]
    if (length(a) != 2L) stop("requires biallelic variant identities")
    paste(z[1], z[2], paste(sort(a), collapse = "/"), sep = ":")
  }, character(1))
}
combine_routes <- function(prepared) {
  d <- prepared$records
  d$source_module <- "8"
  d$source_file <- INPUT
  d$source_record_ids <- d$candidate_id
  d$inclusion_status <- "Primary"
  d$evidence_route <- "Molecular"
  d$display_group <- "Primary"
  d$molecular_status <- "original_output_not_finalized"
  d$link_sample_tissues <- NA_character_
  d$link_support_scope <- "protein_coding_partner_links_within_molecular_block"
  d$raw_source_candidate_id <- NA_character_
  d$raw_matched_p <- NA_real_
  d$raw_route_status <- "not_matched"
  d$selection_reason <- "regional_lead_P_lt_1e-5_and_both_CS_and_LINK"
  d$predicted_biotype <- NA_character_

  raw <- read_source(RAW_INPUT)
  raw$source_row <- seq_len(nrow(raw))
  raw_key <- paste(raw$new_gene_id, raw$tissue, raw$gwas_trait,
                   variant_identity(raw$eqtl_raw_gwas_variant), sep = "|")
  if (anyDuplicated(raw_key)) stop("No. 7 module identity duplicate key, cannot automatic select one records source")
  key <- paste(d$new_gene_id, d$eqtl_tissue, d$gwas_trait,
               variant_identity(d$normalized_variant_key), sep = "|")
  m <- match(key, raw_key)
  raw$plot_column <- match(raw_key, key)
  raw$figure_action <- ifelse(is.na(raw$plot_column), "not_selected_for_this_figure", "annotation_only_no_duplicate_column")
  for (i in which(!is.na(m))) {
    r <- raw[m[i], ]
    d$source_module[i] <- "8;7"
    d$raw_source_candidate_id[i] <- r$candidate_id
    d$predicted_biotype[i] <- r$new_gene_feelnc_biotypes
    d$raw_matched_p[i] <- as.numeric(r$gwas_p_lrt)
    if (!is.finite(d$raw_matched_p[i])) stop("No. 7 module matching P value missing")
    d$raw_route_status[i] <- if (d$raw_matched_p[i] < LEAD_P_THRESHOLD) "Pass" else "Exploratory"
  }
  pos <- read_source(POSITION_INPUT)
  pos$source_row <- seq_len(nrow(pos))

  snp <- sub(" .*", "", pos$gwas_variant_display)
  pos$key <- paste(pos$manual_new_gene, pos$gwas_trait, snp, pos$gwas_locus_id, sep = "|")
  selected <- pos$gwas_branch %in% "main" &
    !is.na(pos$manual_new_gene) & grepl("^DOG5_NEWGENE_[0-9]+$", pos$manual_new_gene) &
    !is.na(pos$gwas_credible_sets) & nzchar(pos$gwas_credible_sets)
  pos$figure_action <- ifelse(selected, "retained_from_selected_focus_file", "not_main_newgene_CS_relation")
  selected_keys <- sort(unique(pos$key[selected]))
  locus <- read_source(LOCUS_INPUT)
  locus <- locus[locus$branch %in% "main", , drop = FALSE]
  if (anyDuplicated(locus$locus_id)) stop("Primary-region manifest contains a duplicate region key")

  for (k in selected_keys) {
    rows <- which(selected & pos$key == k)
    if (!length(rows)) stop("specified exploratory relationship not in input in: ", k)
    z <- pos[rows, , drop = FALSE]
    first <- z[1, ]
    j <- match(first$gwas_locus_id, locus$locus_id)
    if (is.na(j)) stop("cannot find corresponding region lead: ", first$gwas_locus_id)
    if (any(z$gwas_branch != "main") || anyNA(z$gwas_credible_sets)) stop("position path missing main/CS identity")

    e <- d[1, , drop = FALSE]
    e[] <- lapply(e, function(x) x[NA_integer_])
    e$candidate_id <- paste0("module5::", k)
    e$new_gene_id <- first$manual_new_gene
    e$gwas_trait <- first$gwas_trait
    e$gwas_locus_id <- first$gwas_locus_id
    e$gwas_branch <- first$gwas_branch
    e$gwas_analysis_scope <- first$gwas_analysis_scope
    e$gwas_overlap_snp_id <- sub(" .*", "", first$gwas_variant_display)
    e$normalized_variant_key <- variant_identity(sub("^.* \\| ", "", first$gwas_variant_display))
    e$gwas_overlap_snp_p_lrt <- first$gwas_variant_p_lrt
    e$gwas_overlap_snp_pip <- first$gwas_pip
    e$gwas_locus_lead_p_lrt <- locus$gwas_lead_p_lrt[j]
    e$gwas_locus_lead_snp_id <- locus$gwas_lead_variant_id[j]
    e$gwas_overlap_snp_in_credible_set <- "yes"
    e$gwas_credible_sets <- first$gwas_credible_sets
    e$allelome_block_match <- NA_character_
    e$multilayer_integration_scope <- "GWAS_CS_variant_in_newgene_span_with_LINK"
    e$allelome_evidence_strata <- collapse_values(z$evidence_stratum)
    e$source_module <- "5"
    e$source_file <- POSITION_INPUT
    e$source_record_ids <- collapse_values(z$candidate_key)
    e$source_row <- paste(z$source_row, collapse = ";")
    e$inclusion_status <- if (as.numeric(locus$gwas_lead_p_lrt[j]) < LEAD_P_THRESHOLD) "Primary" else "Exploratory"
    e$evidence_route <- "Positional"
    e$display_group <- paste("Positional", e$inclusion_status)
    e$molecular_status <- "not_used_in_this_route"
    e$link_sample_tissues <- collapse_values(z$sample_tissue)
    e$link_support_scope <- "selected_position_overlap_pc_pairs"
    e$selection_reason <- if (e$inclusion_status == "Primary") "position_route_regional_lead_P_lt_1e-5" else "user_retained_original_focus_exploratory_relation"
    e$predicted_biotype <- first$new_gene_feelnc_biotype
    e$raw_route_status <- "-"
    e$pc_link_review <- any(grepl("multigene_review", z$evidence_stratum[!is.na(z$protein_coding_gene_id)]), na.rm = TRUE)
    e$coding_partner_count <- length(unique(z$protein_coding_gene_id[!is.na(z$protein_coding_gene_id)]))

    pc_rows <- !is.na(z$protein_coding_gene_id) & !is.na(z$sample_tissue)
    pc_dogs <- sort(unique(sub("\\|.*", "", z$sample_tissue[pc_rows])))
    e$coding_link_dogs <- paste(pc_dogs, collapse = "|")
    e$coding_link_dog_count <- length(pc_dogs)
    e$display_candidate <- sub("^DOG5_NEWGENE_", "", e$new_gene_id)
    e$display_tissue <- "-"
    e$display_trait <- gsub("_", " ", sub("^score_", "", first$gwas_trait))
    substr(e$display_trait, 1, 1) <- toupper(substr(e$display_trait, 1, 1))
    if (e$display_trait == "Separation related problems") e$display_trait <- "Separation-related problems"
    e$display_record <- e$display_candidate
    e$display_coding_partners <- as.character(e$coding_partner_count)
    e$display_dogs <- as.character(e$coding_link_dog_count)
    e$display_lead_p <- sprintf("%.2e", as.numeric(e$gwas_locus_lead_p_lrt))
    e$display_variant_p <- sprintf("%.2e", as.numeric(e$gwas_overlap_snp_p_lrt))
    e$plot_gwas_cs <- 1L
    e$plot_pc_link <- 1L
    e$plot_any_link <- 1L
    e$plot_link_class <- "pc"
    e$plot_pip <- as.numeric(e$gwas_overlap_snp_pip)
    e$display_pip <- format_probability(e$plot_pip, TRUE)
    e$passes_regional_lead_p_lt_1e5 <- as.numeric(e$gwas_locus_lead_p_lrt) < LEAD_P_THRESHOLD
    e$passes_new_gene <- TRUE
    e$passes_gwas_cs <- TRUE
    e$passes_matched_link <- TRUE
    d <- rbind(d, e)
  }
  d$plot_column <- seq_len(nrow(d))
  d$plot_row <- seq_len(nrow(d))
  d$record_id <- sprintf("R%02d", d$plot_row)

  region_key <- paste(d$gwas_trait, d$gwas_locus_id, sep = "|")
  d$shared_gwas_region <- duplicated(region_key) | duplicated(region_key, fromLast = TRUE)
  prepared$records <- d
  prepared$position_rows <- pos
  prepared$raw_rows <- raw
  prepared
}


ROW_MM <- 6
GROUP_GAP_MM <- 4.5

draw_evidence <- function(d, height_mm) {
  grid.newpage()
  text_at <- function(label, x, y, size = FONT_PT, just = "left",
                      colour = "#222222", face = "plain") {
    grid.text(label, unit(x, "mm"), unit(height_mm-y, "mm"), just = just,
              gp = gpar(fontfamily = FONT, fontsize = size, col = colour, fontface = face))
  }

  p_at <- function(value, x, y) {
    if (is.na(value) || value == "-") return(text_at("-", x, y, just="centre"))
    parts <- strsplit(value, "e", fixed=TRUE)[[1]]
    mantissa <- paste0(parts[1], " \u00d7 10")
    exponent <- as.character(as.integer(parts[2]))
    w1 <- convertWidth(grobWidth(textGrob(mantissa, gp=gpar(fontfamily=FONT, fontsize=7))), "mm", valueOnly=TRUE)
    w2 <- convertWidth(grobWidth(textGrob(exponent, gp=gpar(fontfamily=FONT, fontsize=6.5))), "mm", valueOnly=TRUE)
    left <- x-(w1+w2)/2
    text_at(mantissa, left, y, size=7, colour="#444444")
    text_at(exponent, left+w1, y-1.05, size=6.5, colour="#444444")
  }
  line_at <- function(x0, x1, y0, y1 = y0, colour = "#E6E6E6", lwd = 0.4) {
    grid.segments(unit(x0, "mm"), unit(height_mm-y0, "mm"),
                  unit(x1, "mm"), unit(height_mm-y1, "mm"), gp = gpar(col = colour, lwd = lwd))
  }
  starts <- c(TRUE, head(d$display_group, -1) != tail(d$display_group, -1))

  ys <- 18.5 + (seq_len(nrow(d))-1)*ROW_MM + (cumsum(starts)-1)*GROUP_GAP_MM
  axis_y <- max(ys) + 2.8


  xleft <- 112; xright <- 145
  pip_max <- max(0.20, ceiling(max(d$plot_pip, na.rm = TRUE)/0.05)*0.05)
  pip_x <- function(v) xleft + (xright-xleft)*v/pip_max
  text_at("A", 1.76, 2.8, size = 10, face = "bold")
  text_at("Candidate locus / trait", 6, 6, size = 7, face = "bold")
  text_at("DOG5_NEWGENE_", 6, 9.5, size = 6.5, colour = "#555555")

  text_at("GWAS P", 69, 6, size = 7, just = "centre", face = "bold")
  text_at("Regional lead", 58, 9.5, size = 7, just = "centre")
  text_at("Matched SNP", 80, 9.5, size = 7, just = "centre")
  text_at("Mol. PPA", 99, 7, size = 7, just = "centre", face = "bold")
  text_at("Matched-variant\nGWAS PIP", (xleft+xright)/2, 7, size = 7, just = "centre", face = "bold")
  text_at("PC-link support", 162, 6, size = 7, just = "centre", face = "bold")

  support_header_y <- 9.5
  text_at("Partners", 155, support_header_y, size = 6.5, just = "centre")
  text_at("Dogs", 169, support_header_y, size = 6.5, just = "centre")
  line_at(6, 175, 12, colour = "#BBBBBB")

  for (i in seq_len(nrow(d))) {
    if (!starts[i]) line_at(6, 175, ys[i]-ROW_MM/2)
  }
  for (v in seq(0, pip_max, length.out = 5)) {
    x <- pip_x(v)
    for (g in unique(d$display_group)) {
      yy <- ys[d$display_group == g]
      line_at(x, x, min(yy)-2, max(yy)+2, colour = "#EFEFEF")
    }
    line_at(x, x, axis_y, axis_y+1, colour = "#666666")
    text_at(if (v == 0) "0" else sprintf("%.2f", v), x, axis_y+1.8,
            size = 6.5, just = "centre", colour = "#555555")
  }
  line_at(xleft, xright, axis_y, colour = "#666666")
  for (i in seq_len(nrow(d))) {
    y <- ys[i]
    if (starts[i]) {
      route <- if (d$evidence_route[i] == "Molecular")
        "Molecular-GWAS CS overlap" else "GWAS CS variant within candidate span"
      text_at(paste0(route, "  |  ", d$inclusion_status[i], " regions"),
              6, y-4.5, size = 7, colour = "#555555", face = "plain")
    }

    label <- sub("^DOG5_NEWGENE_", "", d$new_gene_id[i])
    if (d$evidence_route[i] == "Molecular") label <- paste(label, d$display_tissue[i], sep = " / ")

    text_at(label, 6, y-1.3, size = 7, face = "bold")
    text_at(d$display_trait[i], 6, y+1.3, size = 7)
    p_at(d$display_variant_p[i], 80, y)
    p_at(d$display_lead_p[i], 58, y)
    text_at(if (is.na(d$plot_ppa[i])) "-" else d$display_ppa[i], 99, y,
            size = 7, just = "centre", colour = "#444444")
    if (is.na(d$plot_pip[i])) {
      text_at("-", (xleft+xright)/2, y, just = "centre")
    } else {
      x <- pip_x(d$plot_pip[i])
      grid.points(unit(x, "mm"), unit(height_mm-y, "mm"), pch = 16, size = unit(1.7, "mm"),
                  gp = gpar(col = BLUE))
      to_left <- x > xright-10
      text_at(sprintf("%.3g", d$plot_pip[i]), x + if (to_left) -1.8 else 1.8, y,
              size = 6.5, just = if (to_left) "right" else "left")
    }
    text_at(d$display_coding_partners[i], 155, y, size = 7, just = "centre")
    text_at(d$display_dogs[i], 169, y, size = 7, just = "centre")
  }
  text_at("- Molecular evidence not used for nomination", 6, axis_y+1.8, size = 6.5, colour = "#666666")

}

main <- function() {
  run_id <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
  run_dir <- file.path(OUTPUT_ROOT, run_id)
  if (!dir.create(run_dir, recursive = FALSE)) stop("cannot create independent run directory: ", run_dir)
  log_file <- file.path(run_dir, "plot_figure5A.log")
  con <- file(log_file, open = "wt", encoding = "UTF-8")
  sink(con, type = "output"); sink(con, type = "message")
  on.exit({sink(type = "message"); sink(type = "output"); close(con)}, add = TRUE)
  status <- tryCatch({
    log_stage("[1/4] start; concurrent_tasks=1; threads_per_task=1; R/grid independent plot.")
    log_stage(paste("input=", INPUT)); log_stage(paste("run_dir=", run_dir))
    if (!capabilities("cairo")) stop("current R missing Cairo support, cannot export PDF/PNG/TIFF.")
    font_paths <- systemfonts::match_fonts(c(FONT,FONT), weight=c("normal","bold"))$path
    expected_fonts <- c("/path/to/fonts/Arial/Arial.TTF",
                        "/path/to/fonts/Arial/Arialbd.TTF")
    if (!all(file.exists(expected_fonts)) ||
        !identical(normalizePath(font_paths), normalizePath(expected_fonts)))
      stop("Arial regular or bold font is unavailable; do not silently substitute fonts.")
    prepared <- combine_routes(prepare_records())
    log_stage(paste("partner_input=", PARTNER_INPUT))
    d <- prepared$records
    excluded <- prepared$excluded
    log_stage(sprintf("[2/4] No. 8 module%d input records; retain%d records primary candidate, %d record No. 5 module exploratory relationship.No. 7 module only increase source annotation.",
      prepared$input_count, sum(d$inclusion_status == "Primary"), sum(d$inclusion_status == "Exploratory")))
    log_stage(paste("position_input=", POSITION_INPUT))
    log_stage(paste("raw_input=", RAW_INPUT))
    log_stage(paste("locus_lead_input=", LOCUS_INPUT))

    height_mm <- 26.5 + ROW_MM*(nrow(d)-1) + GROUP_GAP_MM*(length(unique(d$display_group))-1)
    prefix <- file.path(run_dir, "Figure5A_integrated_candidate_evidence")
    write.table(d, paste0(prefix, "_source.tsv"), sep = "\t", quote = FALSE,
                 row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
    write.table(excluded, paste0(prefix, "_excluded.tsv"), sep = "\t", quote = FALSE,
                row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
    write.table(prepared$position_rows, paste0(prefix, "_position_details.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    write.table(prepared$raw_rows, paste0(prefix, "_raw_route_annotations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    log_stage(sprintf("[3/4] export%g mm width PDF, 600 dpi PNG and TIFF.", WIDTH_MM))
    for (kind in c("pdf", "svg", "png", "tiff")) {
      output <- paste0(prefix, ".", kind)
      if (kind == "pdf") cairo_pdf(output, width = WIDTH_MM/25.4, height = height_mm/25.4,
                                    family = "sans", pointsize = 6.5, bg = "white")
      if (kind == "svg") svglite::svglite(output, width=WIDTH_MM/25.4,
                                           height=height_mm/25.4, pointsize=7, bg="white")
      if (kind == "png") png(output, width = WIDTH_MM, height = height_mm, units = "mm",
                               res = 600, type = "cairo", bg = "white")
      if (kind == "tiff") tiff(output, width = WIDTH_MM, height = height_mm, units = "mm",
                                 res = 600, type = "cairo", compression = "lzw", bg = "white")
      tryCatch(draw_evidence(d, height_mm), finally = dev.off())
      log_stage(paste("export=", output))
    }
    writeLines(c(paste0("width_mm=",WIDTH_MM,"; height_mm=",height_mm,"; font=Arial"),
      "font_pt: panel=10; headers/values=7; dense_labels/P_exponents=6.5",
      "placement=Figure5 top full-width panel; intended_scale=1",
      "row_pitch_mm=6; group_gap_mm=4.5; first_row_y_from_top_mm=18.5; support_dog_header=Dogs",
      "assembly_and_caption_page_fit=NOT_CHECKED",
      paste("font_paths",paste(font_paths,collapse=";"))),file.path(run_dir,"Figure5A_layout.txt"))
    log_stage("[4/4] save compact figure legend, detailed evidence description, source and run environment.")

    evidence_notes <- c(
      "Figure 5A. Candidate new loci supported through molecular or positional GWAS integration.",
      sprintf("The panel shows %d integration records: %d primary molecular-route records and %d explicitly retained exploratory positional records. These are selected evidence combinations, not calibrated confidence tiers or an exhaustive candidate inventory.", nrow(d), sum(d$inclusion_status == "Primary"), sum(d$inclusion_status == "Exploratory")),
      sprintf("Primary molecular route (module 8): new-gene identity, main GWAS branch, regional lead P < %g, molecular PPA >= 0.1, membership of both molecular and GWAS credible sets, and candidate/block-level LINK. Matched SNP P is displayed in its own aligned numeric column and retained in the source table; it is not used as an additional molecular-route threshold.", LEAD_P_THRESHOLD),
      "Positional (module 5): relations are read and grouped directly from the specified original focus TSV, restricted to main-branch new-gene records with GWAS CS membership. No candidate list or statistic is hard-coded. Regional lead P < 1e-5 is marked Primary; other retained focus relations are explicitly marked Exploratory under the user-selected file scope, without changing the global threshold. Their identities, regions and original expanded rows are retained in the source tables.",
      "Raw-GWAS route (module 7): existing records are matched by gene, molecular tissue, trait and chromosome:position:unordered-alleles, not by module-specific candidate IDs. Pass denotes matched SNP P < 1e-5; Exploratory denotes the retained near-threshold 02502 record. These annotations remain in the source table and do not create duplicate rows. The 01827 records are not selected. Module 4 is not included.",
      "Molecular denotes molecular/GWAS fine-mapping intersection plus candidate/block LINK; Positional denotes a GWAS CS variant within the new-gene span plus LINK. Molecular PPA and tissue are not used in the positional route and are represented by a dash in the PPA column, not as zero or negative evidence. No molecular block is invented for positional records.",
      "Primary and Exploratory describe regional GWAS association scope, not the main/strict fine-mapping branches; all plotted records use the main branch. Primary denotes regional lead P < 1e-5. Retained exploratory records are exceptions under the selected input-file scope, not a relaxation of the primary threshold.",
      "Regional lead P is read from the main-branch locus_run_manifest.tsv for positional records, not substituted with the credible-set minimum P. Original regional and matched-SNP P values remain in the source table; printed values have three significant digits.",
      "Equal-size blue points show matched-variant GWAS PIP on a shared linear axis, with numeric labels. Molecular PPA is shown in a separate aligned numeric column; a dash indicates that it was not used in the positional nomination route. All displayed records have GWAS CS membership and PC LINK support; redundant status squares are omitted. Point positions are not gene confidence scores and rows are not ranked by PIP or partner counts.",
      "Protein-coding partners counts unique partner gene IDs. The Dogs subcolumn under PC-link support counts distinct dog IDs supporting these protein-coding partner links in both routes; counts are not summed across partners or tissues. Molecular records use the partner table for that candidate/block; positional records pool the selected PC-pair overlap rows across tissues. The scope differs in tissue coverage, so counts are descriptive, not a replication-rate comparison. coding_link_dogs and coding_link_dog_count retain the derived IDs and counts; upstream all-LINK fields are not substituted. Neither number is model sample size or an evaluability denominator. Partner counts are not ranking scores.",
      "Each candidate row displays the matched SNP GWAS P_LRT from gwas_overlap_snp_p_lrt, corresponding to the plotted GWAS PIP; full SNP identifiers remain in the source table; positional rows are grouped by candidate, trait, SNP and region. The shared prefix DOG5_NEWGENE_ appears in the candidate-column header; five-digit labels retain leading zeros (01873 denotes DOG5_NEWGENE_01873). Full identifiers remain in the source tables. Regional lead SNP P is displayed in its own aligned numeric column and may refer to a different SNP from the matched-variant GWAS PIP column. Molecular tissue applies only to the molecular model (CB, cerebellum; FC, frontal cortex). Positional support can span several tissues; all dog-tissue combinations and source rows are retained in the positional detail table. The selected 01178 record is predicted mRNA, so the panel is not labelled exclusively lncRNA.",
      sprintf("Row identities map through plot_row and record_id in %s_source.tsv. Original candidate IDs and provenance are retained. The exclusion table covers module 8; the position_details table preserves module 5 expansion and selection status; raw_route_annotations preserves all four module 7 source records with their figure actions.", basename(prefix)),
      "Shared molecular/GWAS CS membership is not formal colocalisation or a jointly calibrated 95% set. LINK joins gene endpoints and does not independently validate the shared SNP or prove regulatory direction. GWAS P, PIP and CS are related summaries, not independent validations.",
      "Multigene attribution review information for displayed PC links is retained in the source table (pc_link_review), derived from the PC-partner table for molecular records and selected positional rows for positional records. It is not marked by a symbol in the figure. Candidate-level flags can include other links; absence of a flag does not establish causal specificity. Shared GWAS regions remain recorded in the source table (shared_gwas_region and gwas_locus_id), without figure labels. Different candidate genes from the same GWAS region do not constitute independent GWAS signals.",
      "Molecular PPA/CS retain the original PLASMA implementation output. Encoding dependence identified in another block remains unresolved for these candidates; this figure does not establish model stability. No new model, significance test or composite score was calculated.",
      "Layout adapted from the aligned gene-label and separate evidence-panel organization in Schwartzentruber et al. (2021), doi:10.1038/s41588-020-00776-w, author function scoreBarPlotSeparated5(); the author code offers both bar and point variants. Candidate grouping also informed by Weeks et al. (2023), Figure 5, doi:10.1038/s41588-023-01443-6. This is an independent R/grid implementation, not reproduction of either published figure, and does not import their composite scores or ranking.",
      "Reference code inspected: https://github.com/Jeremy37/AD_finemap/blob/master/src/plot_gene_evidence.Rmd (scoreBarPlotSeparated5, cowplot::plot_grid with aligned panels). Figure reference: https://www.nature.com/articles/s41588-023-01443-6/figures/5"
    )

    evidence_notes <- c(evidence_notes, "Matched-variant GWAS P values (not regional lead P):",
                 paste(d$display_candidate, d$normalized_variant_key,
                       paste0("P=", d$gwas_overlap_snp_p_lrt), sep = "; "))
    writeLines(evidence_notes, file.path(run_dir, "Figure5A_evidence_notes.txt"), useBytes = TRUE)
    caption <- c(
      "Figure 5A. Candidate new loci supported through molecular or positional GWAS integration.",
      sprintf("The panel summarises %d integration records from %d candidate loci across %d trait-associated GWAS regions. Molecular-route records combine molecular and GWAS credible-set (CS) membership of the matched variant, molecular PPA >= 0.1, and candidate/block-level Allelome LINK to protein-coding partners. Positional-route records combine a GWAS CS variant within the candidate span with LINK support.",
              nrow(d), length(unique(d$new_gene_id)), length(unique(d$gwas_locus_id))),
      sprintf("Primary regions have regional lead SNP P < %g; exploratory records are retained separately. Under the grouped GWAS P header, the Matched SNP subcolumn reports the GWAS P_LRT of the same variant whose posterior inclusion probability (PIP) is shown as a blue point. Full SNP identifiers remain in the source table; P and PIP are related summaries, not independent validation. Regional lead P may refer to a different SNP. Matched SNP does not imply that the variant lies within the candidate span for the molecular route. Mol. PPA denotes molecular posterior probability of association; a dash means molecular evidence was not used for nomination.", LEAD_P_THRESHOLD),
      "Under PC-link support, Partners and Dogs count distinct protein-coding partner genes and dogs supporting these LINK records, respectively. Molecular records use candidate/block tissue support; positional records pool selected tissues, so these counts are not comparable replication rates. Candidate labels share the DOG5_NEWGENE_ prefix. CB, cerebellum; FC, frontal cortex.",
      "Shared CS membership does not establish formal colocalisation, and LINK does not establish causal regulatory direction. Molecular estimates retain the original PLASMA output; the relevance of encoding dependence identified in another block remains unresolved for these candidates. Multigene attribution flags, matched-SNP P values and full provenance are retained in the source tables and Figure5A_evidence_notes.txt."
    )
    writeLines(caption, file.path(run_dir, "Figure5A_caption.txt"), useBytes = TRUE)
    manifest <- data.frame(item = c("script", "input", "partner_input", "input_producer", "run_dir", "backend",
                                     "selection_rule", "filter_field", "filter_operator", "filter_threshold", "module8_input_records", "displayed_records", "module8_excluded_records",
                                     "molecular_result_status", "width_mm", "height_mm", "font_pt",
                                     "concurrent_tasks", "threads_per_task", "rendered_QA"),
      value = c(SCRIPT, INPUT, PARTNER_INPUT, PRODUCER, run_dir, "R grid; R 4.5.1 documentation",
                "module8_primary_bothCS_LINK;module5_original_focus_relations_with_explicit_inclusion_status;module7_annotation_only",
                "gwas_locus_lead_p_lrt", "<", LEAD_P_THRESHOLD, prepared$input_count, nrow(d), nrow(excluded),
                "original_implementation_output_not_finalized", WIDTH_MM, height_mm, FONT_PT,
                "1", "1", "not_performed; review PDF text, collisions and source data separately"))
    manifest <- rbind(manifest, data.frame(
      item = c("position_input", "raw_input", "locus_lead_input", "position_producer", "raw_producer",
               "position_input_expanded_rows", "raw_input_records", "primary_records", "exploratory_records", "module4_included"),
      value = c(POSITION_INPUT, RAW_INPUT, LOCUS_INPUT, POSITION_PRODUCER, RAW_PRODUCER,
                nrow(prepared$position_rows), nrow(prepared$raw_rows), sum(d$inclusion_status == "Primary"),
                sum(d$inclusion_status == "Exploratory"), "no")))
    write.table(manifest, file.path(run_dir, "Figure5A_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
    writeLines(c("execution_status=success", paste0("run_id=", run_id)),
               file.path(run_dir, "execution_completion.marker"))
    log_stage("task completed / done; execution_status=success")
    log_stage("please separately check output completeness, source table corresponding, final font size, text overlap and biological interpretation; execution success does not imply these check already passed.")
    0L
  }, error = function(e) {log_stage(paste("execution_status=failed;", conditionMessage(e))); 1L},
     interrupt = function(e) {log_stage("execution_status=failed; interrupted"); 130L})
  status
}

quit(save = "no", status = main(), runLast = FALSE)
