#!/usr/bin/env Rscript


options(stringsAsFactors = FALSE)

project_dir <- "/path/to/project/legacy/five-dogs/5.cnv"
ase_root <- "/path/to/project/legacy/five-dogs/4.ASE-process/batch-v7/4.6bulk-ase-count"

input_tag <- "keep_locus_genes.all_genes"

casper_expression <- file.path(
  project_dir,
  "5.casper_input",
  paste0("bulk_rna.cnv_expression.", input_tag, ".log2_norm.casper_ordered.tsv")
)

sample_meta <- file.path(
  project_dir,
  "5.casper_input",
  paste0("bulk_rna.samples.", input_tag, ".for_casper.tsv")
)

out_dir <- file.path(project_dir, "6.baf_signal")

min_total_count <- 1L
include_sex_in_primary <- TRUE
include_mito_in_primary <- TRUE

full_dir <- file.path(out_dir, input_tag, "full")
primary_dir <- file.path(out_dir, input_tag, "primary")
dir.create(file.path(full_dir, "per_sample_baf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(full_dir, "per_sample_snp"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(primary_dir, "per_sample_baf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(primary_dir, "per_sample_snp"), recursive = TRUE, showWarnings = FALSE)

out_mapping <- file.path(out_dir, paste0("bulk_rna.", input_tag, ".casper_loh_name_mapping.tsv"))
out_qc <- file.path(out_dir, paste0("bulk_rna.", input_tag, ".baf_from_ase_counts.qc.tsv"))
out_missing_ase <- file.path(out_dir, paste0("bulk_rna.", input_tag, ".missing_ase_tables.tsv"))
out_check <- file.path(out_dir, paste0("prepare_baf_from_ase_counts.", input_tag, ".check.log"))

assert_nonempty_file <- function(path, label) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop(label, " not found or empty: ", path)
  }
}

normalize_chr <- function(chr) {
  x <- sub("^chr", "", chr)
  x[x == "M"] <- "MT"
  x
}

is_autosome <- function(chr) {
  x <- normalize_chr(chr)
  suppressWarnings(num <- as.integer(x))
  !is.na(num) & num >= 1L & num <= 38L
}

is_sex_chr <- function(chr) {
  normalize_chr(chr) %in% c("X", "Y")
}

is_mito_chr <- function(chr) {
  normalize_chr(chr) == "MT"
}

is_primary_chr <- function(chr) {
  primary <- is_autosome(chr)
  if (include_sex_in_primary) {
    primary <- primary | is_sex_chr(chr)
  }
  if (include_mito_in_primary) {
    primary <- primary | is_mito_chr(chr)
  }
  primary
}

chr_order_key <- function(chr) {
  x <- normalize_chr(chr)
  suppressWarnings(num <- as.integer(x))
  ifelse(
    !is.na(num),
    num,
    ifelse(x == "X", 100L, ifelse(x == "Y", 101L, ifelse(x == "MT", 102L, 999L)))
  )
}

find_ase_table <- function(animal_id, term, sample_id) {
  expected <- file.path(
    ase_root,
    animal_id,
    term,
    paste0(sample_id, ".ASEReadCounter.table")
  )
  if (file.exists(expected) && file.info(expected)$size > 0) {
    return(expected)
  }

  sample_dir <- file.path(ase_root, animal_id, term)
  if (!dir.exists(sample_dir)) {
    return(NA_character_)
  }

  candidates <- list.files(
    sample_dir,
    pattern = "\\.ASEReadCounter\\.table$",
    full.names = TRUE
  )
  candidates <- candidates[file.info(candidates)$size > 0]

  if (length(candidates) == 1L) {
    return(candidates[[1]])
  }

  NA_character_
}

read_ase_table <- function(path) {
  x <- read.delim(path, check.names = FALSE)
  required <- c(
    "contig", "position", "refAllele", "altAllele",
    "refCount", "altCount", "totalCount"
  )
  missing <- setdiff(required, colnames(x))
  if (length(missing) > 0) {
    stop(
      "ASEReadCounter table is missing required columns: ",
      paste(missing, collapse = ", "),
      " in ",
      path
    )
  }

  x$position <- as.integer(x$position)
  x$refCount <- as.integer(x$refCount)
  x$altCount <- as.integer(x$altCount)
  x$totalCount <- as.integer(x$totalCount)

  x <- x[
    !is.na(x$position) &
      !is.na(x$totalCount) &
      !is.na(x$altCount) &
      x$totalCount >= min_total_count,
    ,
    drop = FALSE
  ]

  if (nrow(x) == 0) {
    return(x)
  }

  x$baf <- x$altCount / x$totalCount
  x$dev <- abs(x$baf - 0.5)
  x$chr_order <- chr_order_key(x$contig)
  x <- x[order(x$chr_order, x$contig, x$position), , drop = FALSE]

  x
}

write_baf_outputs <- function(x, sample_id, mode_dir, mode_name) {
  baf <- data.frame(
    chr = x$contig,
    position = x$position,
    refAllele = x$refAllele,
    altAllele = x$altAllele,
    ref_count = x$refCount,
    alt_count = x$altCount,
    total_count = x$totalCount,
    baf = round(x$baf, 6),
    dev = round(x$dev, 6),
    stringsAsFactors = FALSE
  )

  baf_file <- file.path(
    mode_dir,
    "per_sample_baf",
    paste0(sample_id, ".baf.", mode_name, ".tsv")
  )
  snp_file <- file.path(
    mode_dir,
    "per_sample_snp",
    paste0(sample_id, ".", mode_name, ".snp")
  )

  write.table(baf, baf_file, sep = "\t", quote = FALSE, row.names = FALSE)


  snp <- data.frame(
    chr = baf$chr,
    position = baf$position,
    refAllele = baf$refAllele,
    altAllele = baf$altAllele,
    alt_count = baf$alt_count,
    total_count = baf$total_count,
    stringsAsFactors = FALSE
  )
  write.table(snp, snp_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  c(baf_file = baf_file, snp_file = snp_file)
}

assert_nonempty_file(casper_expression, "CaSpER expression matrix")
assert_nonempty_file(sample_meta, "CaSpER sample metadata")

expr_header <- read.delim(casper_expression, nrows = 1, check.names = FALSE)
if (ncol(expr_header) < 2 || colnames(expr_header)[1] != "gene_id") {
  stop("CaSpER expression matrix must start with gene_id and contain sample columns.")
}
expr_samples <- colnames(expr_header)[-1]

meta <- read.delim(sample_meta, check.names = FALSE)
required_meta <- c("sample_id", "animal_id", "term")
missing_meta <- setdiff(required_meta, colnames(meta))
if (length(missing_meta) > 0) {
  stop("Sample metadata is missing columns: ", paste(missing_meta, collapse = ", "))
}

if (anyDuplicated(meta$sample_id) > 0) {
  stop("Duplicated sample_id in sample metadata: ", meta$sample_id[duplicated(meta$sample_id)][1])
}

meta <- meta[match(expr_samples, meta$sample_id), , drop = FALSE]
if (!identical(expr_samples, meta$sample_id)) {
  stop("Sample metadata sample_id values do not match CaSpER expression columns.")
}

missing_ase <- data.frame(
  sample_id = character(),
  animal_id = character(),
  term = character(),
  expected_dir = character(),
  stringsAsFactors = FALSE
)

ase_paths <- character(nrow(meta))
for (i in seq_len(nrow(meta))) {
  ase_paths[[i]] <- find_ase_table(meta$animal_id[[i]], meta$term[[i]], meta$sample_id[[i]])
  if (is.na(ase_paths[[i]])) {
    missing_ase <- rbind(
      missing_ase,
      data.frame(
        sample_id = meta$sample_id[[i]],
        animal_id = meta$animal_id[[i]],
        term = meta$term[[i]],
        expected_dir = file.path(ase_root, meta$animal_id[[i]], meta$term[[i]]),
        stringsAsFactors = FALSE
      )
    )
  }
}

write.table(missing_ase, out_missing_ase, sep = "\t", quote = FALSE, row.names = FALSE)
if (nrow(missing_ase) > 0) {
  stop("Missing ASEReadCounter tables were found. See: ", out_missing_ase)
}

qc <- vector("list", nrow(meta))
mapping <- data.frame(loh.name = meta$sample_id, sample.name = meta$sample_id)

for (i in seq_len(nrow(meta))) {
  sample_id <- meta$sample_id[[i]]
  animal_id <- meta$animal_id[[i]]
  term <- meta$term[[i]]
  ase_table <- ase_paths[[i]]

  message("[prepare] ", sample_id)
  baf_all <- read_ase_table(ase_table)
  baf_primary <- baf_all[is_primary_chr(baf_all$contig), , drop = FALSE]

  full_files <- write_baf_outputs(baf_all, sample_id, full_dir, "full")
  primary_files <- write_baf_outputs(baf_primary, sample_id, primary_dir, "primary")

  qc[[i]] <- data.frame(
    sample_id = sample_id,
    animal_id = animal_id,
    term = term,
    ase_table = ase_table,
    full_sites = nrow(baf_all),
    full_sites_total_ge10 = sum(baf_all$totalCount >= 10L),
    full_sites_total_ge20 = sum(baf_all$totalCount >= 20L),
    full_sites_ge20_baf_0.2_0.8 = sum(baf_all$totalCount >= 20L & baf_all$baf > 0.2 & baf_all$baf < 0.8),
    primary_sites = nrow(baf_primary),
    primary_autosome_sites = sum(is_autosome(baf_primary$contig)),
    primary_sex_sites = sum(is_sex_chr(baf_primary$contig)),
    primary_mito_sites = sum(is_mito_chr(baf_primary$contig)),
    primary_sites_total_ge10 = sum(baf_primary$totalCount >= 10L),
    primary_sites_total_ge20 = sum(baf_primary$totalCount >= 20L),
    primary_sites_ge20_baf_0.2_0.8 = sum(baf_primary$totalCount >= 20L & baf_primary$baf > 0.2 & baf_primary$baf < 0.8),
    full_baf_file = unname(full_files[["baf_file"]]),
    primary_baf_file = unname(primary_files[["baf_file"]]),
    full_snp_file = unname(full_files[["snp_file"]]),
    primary_snp_file = unname(primary_files[["snp_file"]]),
    stringsAsFactors = FALSE
  )
}

qc <- do.call(rbind, qc)
write.table(mapping, out_mapping, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(qc, out_qc, sep = "\t", quote = FALSE, row.names = FALSE)

log_lines <- c(
  paste0("status=OK"),
  paste0("input_tag=", input_tag),
  paste0("casper_expression=", casper_expression),
  paste0("sample_meta=", sample_meta),
  paste0("ase_root=", ase_root),
  paste0("min_total_count=", min_total_count),
  paste0("include_sex_in_primary=", include_sex_in_primary),
  paste0("include_mito_in_primary=", include_mito_in_primary),
  paste0("samples=", nrow(meta)),
  paste0("missing_ase_tables=", nrow(missing_ase)),
  paste0("total_full_sites=", sum(qc$full_sites)),
  paste0("total_primary_sites=", sum(qc$primary_sites)),
  paste0("output_full_dir=", full_dir),
  paste0("output_primary_dir=", primary_dir),
  paste0("output_mapping=", out_mapping),
  paste0("output_qc=", out_qc),
  paste0("output_missing_ase=", out_missing_ase)
)
writeLines(log_lines, out_check)
writeLines(log_lines)
