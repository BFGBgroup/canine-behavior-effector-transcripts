#!/usr/bin/env Rscript


options(stringsAsFactors = FALSE, warn = 1)

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_message <- function(stage, message) {
  cat(sprintf("[%s ] [%s] %s\n", timestamp(), stage, message))
  flush.console()
}

abort <- function(message) stop(message, call. = FALSE)
assert_true <- function(condition, message) {
  if (!isTRUE(condition)) abort(message)
}

parse_args <- function(args) {
  assert_true(length(args) %% 2L == 0L, "parameter must follow --key value provide as a pair")
  parsed <- list()
  if (length(args) == 0L) return(parsed)
  for (i in seq.int(1L, length(args), by = 2L)) {
    key <- sub("^--", "", args[[i]])
    assert_true(startsWith(args[[i]], "--"), paste0("cannot identify parameter: ", args[[i]]))
    parsed[[key]] <- args[[i + 1L]]
  }
  parsed
}

require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  assert_true(length(missing) == 0L,
              paste0(label, " missing fields: ", paste(missing, collapse = ", ")))
}

read_tsv <- function(path, label) {
  assert_true(file.exists(path), paste0(label, " does not exist: ", path))
  read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "",
             check.names = FALSE, stringsAsFactors = FALSE,
             na.strings = c("NA", "NaN", "nan"))
}

read_tsv_gz <- function(path, label) {
  assert_true(file.exists(path), paste0(label, " does not exist: ", path))
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  read.delim(connection, sep = "\t", header = TRUE, quote = "", comment.char = "",
             check.names = FALSE, stringsAsFactors = FALSE,
             na.strings = c("NA", "NaN", "nan"))
}

write_tsv <- function(data, path) {
  write.table(data, file = path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE, na = "NA")
}

write_tsv_gz <- function(data, path) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  write.table(data, file = connection, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE, na = "NA")
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required_args <- c("locus-id", "manifest", "ld-prefix", "plink2", "output-dir")
missing_args <- required_args[!required_args %in% names(args)]
assert_true(length(missing_args) == 0L,
            paste0("missing required parameter: ", paste(missing_args, collapse = ", ")))

locus_id <- args[["locus-id"]]
manifest_path <- args[["manifest"]]
ld_prefix <- args[["ld-prefix"]]
plink2_bin <- args[["plink2"]]
output_dir <- args[["output-dir"]]
plink2_threads <- as.integer(ifelse(is.null(args[["plink2-threads"]]), "1", args[["plink2-threads"]]))
frequency_tolerance <- as.numeric(ifelse(is.null(args[["frequency-tolerance"]]), "1e-5", args[["frequency-tolerance"]]))

assert_true(grepl("^/", manifest_path), "manifest must be an absolute path")
assert_true(grepl("^/", ld_prefix), "ld-prefix must be an absolute path")
assert_true(grepl("^/", plink2_bin), "plink2 must be an absolute path")
assert_true(grepl("^/", output_dir), "output-dir must be an absolute path")
assert_true(file.exists(plink2_bin), paste0("PLINK2 does not exist: ", plink2_bin))
assert_true(plink2_threads >= 1L, "plink2-threads must be at least 1")
assert_true(is.finite(frequency_tolerance) && frequency_tolerance > 0,
            "frequency-tolerance must be positive number")
assert_true(!dir.exists(output_dir) && !file.exists(output_dir),
            paste0("refuse to overwrite existing locus output: ", output_dir))

log_message("INPUT", paste0("read locus list: ", manifest_path))
manifest <- read_tsv(manifest_path, "locus input list")
require_columns(
  manifest,
  c("locus_id", "trait", "analysis_scope", "locus_input_status",
    "harmonized_variant_count", "locus_directory", "harmonized_variants",
    "variant_order", "allele_harmonization_audit", "dropped_variants",
    "locus_metadata"),
  "locus input list"
)
locus_row <- manifest[manifest$locus_id == locus_id, , drop = FALSE]
assert_true(nrow(locus_row) == 1L,
            paste0("locus_id must be in list in unique encountered once: ", locus_id))
assert_true(locus_row$locus_input_status[[1L]] == "READY_FOR_GENOTYPE_MATRIX",
            paste0("locus not yet ready: ", locus_row$locus_input_status[[1L]]))

variant_order_path <- locus_row$variant_order[[1L]]
harmonized_path <- locus_row$harmonized_variants[[1L]]
variant_order <- read_tsv(variant_order_path, "variant_order")
harmonized <- read_tsv_gz(harmonized_path, "harmonized_variants")

require_columns(
  variant_order,
  c("variant_index", "locus_id", "analysis_scope", "trait", "chrom", "pos",
    "ld_variant_id", "ld_a1", "ld_a2", "z_for_ld_a1", "n_obs",
    "ambiguous_at_cg", "strict_no_ambiguous_include"),
  "variant_order"
)
require_columns(
  harmonized,
  c("variant_index", "locus_id", "trait", "analysis_scope", "chrom", "pos",
    "ld_variant_id", "ld_bim_order", "ld_a1", "ld_a2", "ld_a1_frequency",
    "beta", "se", "z_for_ld_a1", "p_lrt", "n_obs", "ambiguous_at_cg",
    "harmonization_status", "main_analysis_include", "strict_no_ambiguous_include"),
  "harmonized_variants"
)

assert_true(nrow(variant_order) >= 2L, "locus must contain at least 2 variants")
assert_true(!anyDuplicated(variant_order$variant_index), "variant_order in variant_index duplicate")
assert_true(!anyDuplicated(variant_order$ld_variant_id), "variant_order in ld_variant_id duplicate")
assert_true(all(variant_order$locus_id == locus_id), "variant_order unexpected inclusion of other locus")
assert_true(all(variant_order$trait == locus_row$trait[[1L]]) &&
                all(variant_order$analysis_scope == locus_row$analysis_scope[[1L]]),
            "trait or analysis_scope in variant_order differs from the manifest")
expected_variant_count <- suppressWarnings(as.integer(locus_row$harmonized_variant_count[[1L]]))
assert_true(is.finite(expected_variant_count) && expected_variant_count >= 2L,
            "harmonized_variant_count in the manifest is not a valid integer")
assert_true(nrow(variant_order) == expected_variant_count,
            "variant_order row count differs from harmonized_variant_count in the manifest")
assert_true(!anyDuplicated(harmonized$variant_index),
            "harmonized_variants in variant_index duplicate")

match_index <- match(variant_order$variant_index, harmonized$variant_index)
assert_true(!anyNA(match_index), "variant_order in exists cannot trace back to harmonized_variants records")
harmonized <- harmonized[match_index, , drop = FALSE]
assert_true(nrow(harmonized) == nrow(variant_order), "two copies variant table row counts differ")
assert_true(all(harmonized$locus_id == locus_id) &&
                all(harmonized$trait == locus_row$trait[[1L]]) &&
                all(harmonized$analysis_scope == locus_row$analysis_scope[[1L]]),
            "locus_id, trait or analysis_scope in harmonized_variants differs from the manifest")
assert_true(all(harmonized$ld_variant_id == variant_order$ld_variant_id),
            "Variant order or ID differs from harmonized_variants")
assert_true(all(harmonized$ld_a1 == variant_order$ld_a1 &
                harmonized$ld_a2 == variant_order$ld_a2),
            "LD alleles differ from variant_order")
assert_true(all(harmonized$main_analysis_include == "yes"),
            "variant_order contains variants excluded from the main analysis")
assert_true(all(variant_order$strict_no_ambiguous_include %in% c("yes", "no")) &&
                all(harmonized$strict_no_ambiguous_include %in% c("yes", "no")),
            "strict_no_ambiguous_include contains non-yes/no value")
assert_true(all(is.finite(as.numeric(harmonized$ld_bim_order)) &
                as.numeric(harmonized$ld_bim_order) >= 1),
            "ld_bim_order contains invalid value")
assert_true(all(is.finite(as.numeric(harmonized$ld_a1_frequency)) &
                as.numeric(harmonized$ld_a1_frequency) >= 0 &
                as.numeric(harmonized$ld_a1_frequency) <= 1),
            "ld_a1_frequency exists exceeds 0 to 1 value")

bed_path <- paste0(ld_prefix, ".bed")
bim_path <- paste0(ld_prefix, ".bim")
fam_path <- paste0(ld_prefix, ".fam")
for (path in c(bed_path, bim_path, fam_path)) {
  assert_true(file.exists(path), paste0("PLINK reference file does not exist: ", path))
}

fam <- read.table(fam_path, header = FALSE, stringsAsFactors = FALSE,
                  col.names = c("FID", "IID", "father", "mother", "sex", "phenotype"),
                  colClasses = "character", quote = "", comment.char = "")
assert_true(nrow(fam) >= 2L, "FAM sample count insufficient")
fam$sample_key <- paste(fam$FID, fam$IID, sep = ":")
assert_true(!anyDuplicated(fam$sample_key), "FAM in FID+IID is not unique")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
intermediate_dir <- file.path(output_dir, "intermediate")
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
extract_path <- file.path(intermediate_dir, "extract_snps.txt")
writeLines(variant_order$ld_variant_id, extract_path, useBytes = TRUE)

plink_out_prefix <- file.path(intermediate_dir, "plink_export")
plink_args <- c(
  "--bfile", ld_prefix,
  "--dog",
  "--allow-extra-chr",
  "--extract", extract_path,
  "--threads", as.character(plink2_threads),
  "--export", "A", "include-alt",
  "--out", plink_out_prefix
)

log_message("PLINK2", paste0("export ", nrow(variant_order), " items locus additive genotype for variant"))
plink_output <- system2(plink2_bin, args = plink_args, stdout = TRUE, stderr = TRUE)
if (length(plink_output) > 0L) cat(paste(plink_output, collapse = "\n"), "\n")
plink_status <- attr(plink_output, "status")
if (is.null(plink_status)) plink_status <- 0L
assert_true(plink_status == 0L, paste0("PLINK2 export failed, exit code: ", plink_status))

raw_path <- paste0(plink_out_prefix, ".raw")
assert_true(file.exists(raw_path), paste0("PLINK2 not generated .raw: ", raw_path))
raw_data <- read.table(raw_path, header = TRUE, quote = "", comment.char = "",
                       check.names = FALSE, stringsAsFactors = FALSE,
                       na.strings = c("NA", "nan", "NaN"))

fid_column <- if ("#FID" %in% names(raw_data)) "#FID" else if ("FID" %in% names(raw_data)) "FID" else NA_character_
assert_true(!is.na(fid_column) && "IID" %in% names(raw_data),
            "PLINK2 .raw missing FID/IID field")
raw_sample_key <- paste(raw_data[[fid_column]], raw_data[["IID"]], sep = ":")
assert_true(!anyDuplicated(raw_sample_key), "PLINK2 .raw in FID+IID is not unique")
sample_match <- match(fam$sample_key, raw_sample_key)
assert_true(!anyNA(sample_match) && nrow(raw_data) == nrow(fam),
            "PLINK2 exported samples differ from the FAM")
raw_data <- raw_data[sample_match, , drop = FALSE]


raw_names <- names(raw_data)
genotype_columns <- vapply(variant_order$ld_variant_id, function(variant_id) {
  column_prefix <- paste0(variant_id, "_")
  candidates <- which(startsWith(raw_names, column_prefix))
  assert_true(length(candidates) == 1L,
              paste0("cannot unique mapping PLINK2 genotype column: ", variant_id))
  header_suffix <- sub(paste0("^", variant_id, "_"), "", raw_names[candidates[[1L]]])
  allele_match <- regexec("^([ACGTN]+)\\(/([ACGTN]+)\\)$",
                          header_suffix, perl = TRUE)
  allele_tokens <- regmatches(header_suffix, allele_match)[[1L]]
  current_variant <- harmonized[harmonized$ld_variant_id == variant_id, , drop = FALSE]
  header_ref <- if (length(allele_tokens) == 3L) allele_tokens[[2L]] else NA_character_
  header_alt <- if (length(allele_tokens) == 3L) allele_tokens[[3L]] else NA_character_
  assert_true(nrow(current_variant) == 1L &&
                length(allele_tokens) == 3L &&
                header_ref == current_variant$ld_a2[[1L]] &&
                header_alt == current_variant$ld_a1[[1L]],
              paste0("PLINK2 .raw column header REF/ALT with BIM A2/A1 does not match: ",
                     raw_names[candidates[[1L]]]))
  candidates[[1L]]
}, integer(1L))

ref_dosage <- as.matrix(data.frame(
  lapply(raw_data[, genotype_columns, drop = FALSE], as.numeric),
  check.names = FALSE
))
colnames(ref_dosage) <- variant_order$ld_variant_id
assert_true(all(is.na(ref_dosage) | ref_dosage %in% c(0, 1, 2)),
            "PLINK2 additive genotypes contain 0/1/2/NA outside value")
a1_dosage <- 2 - ref_dosage

expected_a1_frequency <- as.numeric(harmonized$ld_a1_frequency)
observed_a1_frequency <- colMeans(a1_dosage, na.rm = TRUE) / 2
frequency_difference <- abs(observed_a1_frequency - expected_a1_frequency)
assert_true(all(is.finite(frequency_difference)), "A1 frequency validation found non-finite values")
assert_true(max(frequency_difference) <= frequency_tolerance,
            paste0("A1 dosage orientation validate failed, maximum frequency difference = ",
                   format(max(frequency_difference), scientific = TRUE),
                   ", tolerance = ", frequency_tolerance))

missing_before <- is.na(a1_dosage)
column_means <- colMeans(a1_dosage, na.rm = TRUE)
observed_variance <- vapply(seq_len(ncol(a1_dosage)), function(index) {
  stats::var(a1_dosage[, index], na.rm = TRUE)
}, numeric(1L))
keep_nonzero <- is.finite(column_means) & is.finite(observed_variance) & observed_variance > 0
original_main_variant_count <- nrow(variant_order)
original_strict_variant_count <- sum(variant_order$strict_no_ambiguous_include == "yes")
strict_zero_variance_count <- sum(!keep_nonzero &
                                  variant_order$strict_no_ambiguous_include == "yes")

zero_variance_index <- which(!keep_nonzero)
zero_variance <- data.frame(
  locus_id = rep(locus_id, length(zero_variance_index)),
  variant_index = variant_order$variant_index[zero_variance_index],
  ld_variant_id = variant_order$ld_variant_id[zero_variance_index],
  reason = rep("ZERO_OR_UNDEFINED_VARIANCE", length(zero_variance_index)),
  stringsAsFactors = FALSE
)
write_tsv(zero_variance, file.path(output_dir, "zero_variance_variants.tsv"))
assert_true(sum(keep_nonzero) >= 2L, "remove zero-variance variant after insufficient 2 variants")

a1_dosage <- a1_dosage[, keep_nonzero, drop = FALSE]
missing_before <- missing_before[, keep_nonzero, drop = FALSE]
column_means <- column_means[keep_nonzero]
variant_order <- variant_order[keep_nonzero, , drop = FALSE]
harmonized <- harmonized[keep_nonzero, , drop = FALSE]
frequency_difference <- frequency_difference[keep_nonzero]

for (column_index in seq_len(ncol(a1_dosage))) {
  missing_index <- is.na(a1_dosage[, column_index])
  if (any(missing_index)) {
    a1_dosage[missing_index, column_index] <- column_means[[column_index]]
  }
}
assert_true(!anyNA(a1_dosage), "all value impute after still contains missing genotype")
rownames(a1_dosage) <- fam$sample_key
colnames(a1_dosage) <- variant_order$ld_variant_id

sample_order <- data.frame(
  sample_index = seq_len(nrow(fam)),
  FID = fam$FID,
  IID = fam$IID,
  sample_key = fam$sample_key,
  stringsAsFactors = FALSE
)
write_tsv(sample_order, file.path(output_dir, "sample_order.tsv"))

build_branch <- function(branch_name, include, original_input_count,
                         branch_zero_variance_count) {
  branch_dir <- file.path(output_dir, branch_name)
  dir.create(branch_dir, recursive = TRUE, showWarnings = FALSE)
  assert_true(sum(include) >= 2L,
              paste0(branch_name, " branch less than 2 variants"))

  x_matrix <- a1_dosage[, include, drop = FALSE]
  branch_order <- variant_order[include, , drop = FALSE]
  branch_harmonized <- harmonized[include, , drop = FALSE]
  branch_frequency_difference <- frequency_difference[include]
  branch_missing_before <- missing_before[, include, drop = FALSE]

  r_matrix <- stats::cor(x_matrix)
  assert_true(all(is.finite(r_matrix)), paste0(branch_name, " LD matrix contains non-finite values"))
  symmetry_error <- max(abs(r_matrix - t(r_matrix)))
  diagonal_error <- max(abs(diag(r_matrix) - 1))
  minimum_eigenvalue <- min(eigen(r_matrix, symmetric = TRUE, only.values = TRUE)$values)
  assert_true(symmetry_error <= 1e-10,
              paste0(branch_name, " LD matrix asymmetric: ", symmetry_error))
  assert_true(diagonal_error <= 1e-10,
              paste0(branch_name, " LD matrix diagonal deviation 1: ", diagonal_error))
  assert_true(minimum_eigenvalue >= -1e-8,
              paste0(branch_name, " LD matrix exists substantially negative eigenvalue: ", minimum_eigenvalue))

  final_order <- data.frame(
    matrix_variant_index = seq_len(nrow(branch_order)),
    branch_order,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  z_beta_se <- data.frame(
    matrix_variant_index = seq_len(nrow(branch_harmonized)),
    locus_id = branch_harmonized$locus_id,
    trait = branch_harmonized$trait,
    analysis_scope = branch_harmonized$analysis_scope,
    chrom = branch_harmonized$chrom,
    pos = branch_harmonized$pos,
    ld_variant_id = branch_harmonized$ld_variant_id,
    ld_a1 = branch_harmonized$ld_a1,
    ld_a2 = branch_harmonized$ld_a2,
    beta_for_gwas_allele1 = branch_harmonized$beta,
    beta_for_ld_a1 = as.numeric(branch_harmonized$z_for_ld_a1) *
      as.numeric(branch_harmonized$se),
    se = branch_harmonized$se,
    z_for_ld_a1 = branch_harmonized$z_for_ld_a1,
    p_lrt = branch_harmonized$p_lrt,
    n_obs = branch_harmonized$n_obs,
    ambiguous_at_cg = branch_harmonized$ambiguous_at_cg,
    harmonization_status = branch_harmonized$harmonization_status,
    stringsAsFactors = FALSE
  )

  saveRDS(x_matrix, file.path(branch_dir, "X_ref.rds"), compress = "gzip")
  saveRDS(r_matrix, file.path(branch_dir, "R_signed.rds"), compress = "gzip")
  write_tsv(final_order, file.path(branch_dir, "variant_order.tsv"))
  write_tsv_gz(z_beta_se, file.path(branch_dir, "z_beta_se.tsv.gz"))

  data.frame(
    locus_id = locus_id,
    trait = locus_row$trait[[1L]],
    analysis_scope = locus_row$analysis_scope[[1L]],
    branch = branch_name,
    status = "PASS",
    sample_count = nrow(x_matrix),
    input_variant_count = original_input_count,
    zero_variance_dropped_count = branch_zero_variance_count,
    output_variant_count = ncol(x_matrix),
    missing_genotype_count_before_imputation = sum(branch_missing_before),
    missing_genotype_fraction_before_imputation = sum(branch_missing_before) / length(branch_missing_before),
    missing_genotype_count_after_imputation = sum(is.na(x_matrix)),
    max_a1_frequency_difference = max(branch_frequency_difference),
    r_symmetry_max_abs = symmetry_error,
    r_diagonal_max_abs_difference = diagonal_error,
    r_min_eigenvalue = minimum_eigenvalue,
    n_locus_median = stats::median(as.numeric(branch_harmonized$n_obs)),
    n_obs_min = min(as.numeric(branch_harmonized$n_obs)),
    n_obs_max = max(as.numeric(branch_harmonized$n_obs)),
    stringsAsFactors = FALSE
  )
}

log_message("MATRIX", "build the main analysis and the sensitivity branch without ambiguous variants")
main_include <- rep(TRUE, nrow(variant_order))
strict_include <- variant_order$strict_no_ambiguous_include == "yes"
qc_main <- build_branch(
  "main", main_include, original_main_variant_count, nrow(zero_variance)
)
qc_strict <- build_branch(
  "strict_no_ambiguous", strict_include, original_strict_variant_count,
  strict_zero_variance_count
)
matrix_qc <- rbind(qc_main, qc_strict)
write_tsv(matrix_qc, file.path(output_dir, "matrix_qc.tsv"))

locus_manifest <- data.frame(
  locus_id = locus_id,
  trait = locus_row$trait[[1L]],
  analysis_scope = locus_row$analysis_scope[[1L]],
  locus_status = "PASS",
  sample_count = nrow(fam),
  main_variant_count = qc_main$output_variant_count,
  strict_variant_count = qc_strict$output_variant_count,
  locus_output = output_dir,
  main_x_ref = file.path(output_dir, "main", "X_ref.rds"),
  main_r_signed = file.path(output_dir, "main", "R_signed.rds"),
  main_z_beta_se = file.path(output_dir, "main", "z_beta_se.tsv.gz"),
  strict_x_ref = file.path(output_dir, "strict_no_ambiguous", "X_ref.rds"),
  strict_r_signed = file.path(output_dir, "strict_no_ambiguous", "R_signed.rds"),
  strict_z_beta_se = file.path(output_dir, "strict_no_ambiguous", "z_beta_se.tsv.gz"),
  stringsAsFactors = FALSE
)
write_tsv(locus_manifest, file.path(output_dir, "locus_matrix_manifest.tsv"))

completion <- c(
  "overall_status\tPASS",
  paste0("locus_id\t", locus_id),
  paste0("trait\t", locus_row$trait[[1L]]),
  paste0("analysis_scope\t", locus_row$analysis_scope[[1L]]),
  paste0("sample_count\t", nrow(fam)),
  paste0("main_variant_count\t", qc_main$output_variant_count),
  paste0("strict_variant_count\t", qc_strict$output_variant_count),
  paste0("completed_at\t", timestamp())
)
writeLines(completion, file.path(output_dir, "COMPLETED.ok"), useBytes = TRUE)
log_message("DONE", paste0("locus completed: ", locus_id))
