#!/usr/bin/env Rscript


options(stringsAsFactors = FALSE, warn = 1)

EXPECTED_SUSIER_VERSION <- "0.16.5"
SUSIE_COVERAGE <- 0.95
SUSIE_MIN_ABS_CORR <- 0.5
SUSIE_MAX_ITER <- 50L
MATRIX_TOLERANCE <- 1e-8

log_message <- function(stage, ...) {
  message(sprintf(
    "[%s ] [%s] %s",
    format(Sys.time(), "%F %T"),
    stage,
    paste0(..., collapse = "")
  ))
}

stop_message <- function(...) {
  stop(paste0(..., collapse = ""), call. = FALSE)
}

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop_message(message)
  }
}

parse_arguments <- function(arguments) {
  if (length(arguments) %% 2L != 0L) {
    stop_message("parameter must follow --name value provide as a pair")
  }
  parsed <- list()
  if (length(arguments) > 0L) {
    for (index in seq.int(1L, length(arguments), by = 2L)) {
      key <- sub("^--", "", arguments[[index]])
      if (!startsWith(arguments[[index]], "--") || !nzchar(key)) {
        stop_message("unrecognized argument name: ", arguments[[index]])
      }
      parsed[[key]] <- arguments[[index + 1L]]
    }
  }
  required <- c(
    "locus-id", "manifest", "qc-summary", "metadata-manifest", "branch",
    "output-dir"
  )
  missing <- required[!vapply(required, function(key) {
    !is.null(parsed[[key]]) && nzchar(parsed[[key]])
  }, logical(1))]
  if (length(missing) > 0L) {
    stop_message("missing required parameter: ", paste(missing, collapse = ", "))
  }
  parsed
}

read_tsv <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(con), add = TRUE)
  read.delim(
    con,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

write_tsv <- function(table, path) {
  write.table(
    table,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}

write_tsv_gz <- function(table, path) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  write.table(
    table,
    file = connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}

check_columns <- function(table, required, label) {
  missing <- setdiff(required, names(table))
  if (length(missing) > 0L) {
    stop_message(label, " missing required columns: ", paste(missing, collapse = ", "))
  }
}

read_key_value_metadata <- function(path, required_keys) {
  metadata_table <- read_tsv(path)
  check_columns(metadata_table, c("key", "value"), path)
  assert_true(!anyDuplicated(metadata_table$key),
              paste0(path, " in key must be unique"))
  missing_keys <- setdiff(required_keys, metadata_table$key)
  assert_true(length(missing_keys) == 0L, paste0(
    path, " missing required key: ", paste(missing_keys, collapse = ", ")
  ))
  stats::setNames(
    as.character(metadata_table$value),
    as.character(metadata_table$key)
  )
}

check_nonempty_file <- function(path, label) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0) {
    stop_message("missing or is empty of input file: ", label, " -> ", path)
  }
}

capture_stage_warnings <- function(stage, expression, warning_store) {
  withCallingHandlers(
    expression,
    warning = function(warning_condition) {
      warning_store$stage <- c(warning_store$stage, stage)
      warning_store$message <- c(
        warning_store$message,
        conditionMessage(warning_condition)
      )
      log_message("WARNING", stage, ": ", conditionMessage(warning_condition))
      invokeRestart("muffleWarning")
    }
  )
}

empty_warning_table <- function() {
  data.frame(
    warning_index = integer(),
    stage = character(),
    message = character(),
    stringsAsFactors = FALSE
  )
}

warning_table <- function(warning_store) {
  if (length(warning_store$message) == 0L) {
    return(empty_warning_table())
  }
  data.frame(
    warning_index = seq_along(warning_store$message),
    stage = warning_store$stage,
    message = warning_store$message,
    stringsAsFactors = FALSE
  )
}

diag_scalar <- function(diagnostics, field, default = NA) {
  value <- diagnostics[[field]]
  if (is.null(value) || length(value) == 0L) {
    return(default)
  }
  value[[1L]]
}

numeric_summary <- function(value) {
  value <- as.numeric(value)
  value <- value[is.finite(value)]
  if (length(value) == 0L) {
    return(c(min = NA_real_, mean = NA_real_, max = NA_real_))
  }
  c(min = min(value), mean = mean(value), max = max(value))
}

purity_value <- function(purity, cs_name, cs_index, column_name) {
  if (is.null(purity) || !(is.matrix(purity) || is.data.frame(purity))) {
    return(NA_real_)
  }
  row_index <- cs_index
  if (!is.null(rownames(purity)) && cs_name %in% rownames(purity)) {
    row_index <- match(cs_name, rownames(purity))
  }
  if (row_index > nrow(purity) || !(column_name %in% colnames(purity))) {
    return(NA_real_)
  }
  as.numeric(purity[row_index, column_name])
}

credible_set_names <- function(cs_list) {
  if (is.null(cs_list) || length(cs_list) == 0L) {
    return(character())
  }
  names_out <- names(cs_list)
  if (is.null(names_out)) {
    names_out <- rep("", length(cs_list))
  }
  missing <- !nzchar(names_out)
  names_out[missing] <- paste0("CS", which(missing))
  names_out
}

credible_set_labels <- function(cs_list, cs_names, variant_count) {
  labels <- rep("", variant_count)
  if (is.null(cs_list) || length(cs_list) == 0L) {
    return(labels)
  }
  for (index in seq_along(cs_list)) {
    members <- as.integer(cs_list[[index]])
    current <- labels[members]
    labels[members] <- ifelse(
      nzchar(current),
      paste(current, cs_names[[index]], sep = ","),
      cs_names[[index]]
    )
  }
  labels
}

empty_credible_set_table <- function() {
  data.frame(
    locus_id = character(), trait = character(), analysis_scope = character(),
    branch = character(), cs_index = integer(), cs_id = character(),
    susie_effect_index = integer(), cs_size = integer(), coverage = numeric(),
    min_abs_corr = numeric(), mean_abs_corr = numeric(),
    median_abs_corr = numeric(), cs_min_pos = integer(), cs_max_pos = integer(),
    lead_ld_variant_id = character(), lead_pos = integer(), lead_pip = numeric(),
    lead_z_for_ld_a1 = numeric(), lead_p_lrt = numeric(), cs_pip_sum = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_credible_set_member_table <- function() {
  data.frame(
    locus_id = character(), trait = character(), analysis_scope = character(),
    branch = character(), cs_index = integer(), cs_id = character(),
    susie_effect_index = integer(), member_index = integer(),
    matrix_variant_index = integer(), chrom = character(), pos = integer(),
    ld_variant_id = character(), ld_a1 = character(), ld_a2 = character(),
    pip = numeric(), z_for_ld_a1 = numeric(), p_lrt = numeric(),
    stringsAsFactors = FALSE
  )
}

extract_credible_sets <- function(fit, z_table, metadata) {
  cs_list <- fit$sets$cs
  if (is.null(cs_list) || length(cs_list) == 0L) {
    return(list(
      summary = empty_credible_set_table(),
      members = empty_credible_set_member_table()
    ))
  }

  cs_names <- credible_set_names(cs_list)
  effect_indices <- fit$sets$cs_index
  if (is.null(effect_indices) || length(effect_indices) != length(cs_list)) {
    effect_indices <- rep(NA_integer_, length(cs_list))
  }
  summary_rows <- vector("list", length(cs_list))
  member_rows <- list()
  member_row_index <- 1L

  for (cs_index in seq_along(cs_list)) {
    members <- as.integer(cs_list[[cs_index]])
    assert_true(
      length(members) > 0L && all(members >= 1L) && all(members <= nrow(z_table)),
      paste0("credible set ", cs_names[[cs_index]], " contains out of bounds member")
    )
    cs_name <- cs_names[[cs_index]]
    lead_index <- members[[which.max(fit$pip[members])]]
    member_table <- z_table[members, , drop = FALSE]
    summary_rows[[cs_index]] <- data.frame(
      locus_id = metadata$locus_id,
      trait = metadata$trait,
      analysis_scope = metadata$analysis_scope,
      branch = metadata$branch,
      cs_index = cs_index,
      cs_id = cs_name,
      susie_effect_index = as.integer(effect_indices[[cs_index]]),
      cs_size = length(members),
      coverage = SUSIE_COVERAGE,
      min_abs_corr = purity_value(
        fit$sets$purity, cs_name, cs_index, "min.abs.corr"
      ),
      mean_abs_corr = purity_value(
        fit$sets$purity, cs_name, cs_index, "mean.abs.corr"
      ),
      median_abs_corr = purity_value(
        fit$sets$purity, cs_name, cs_index, "median.abs.corr"
      ),
      cs_min_pos = min(as.integer(member_table$pos)),
      cs_max_pos = max(as.integer(member_table$pos)),
      lead_ld_variant_id = z_table$ld_variant_id[[lead_index]],
      lead_pos = as.integer(z_table$pos[[lead_index]]),
      lead_pip = as.numeric(fit$pip[[lead_index]]),
      lead_z_for_ld_a1 = as.numeric(z_table$z_for_ld_a1[[lead_index]]),
      lead_p_lrt = as.numeric(z_table$p_lrt[[lead_index]]),
      cs_pip_sum = sum(as.numeric(fit$pip[members])),
      stringsAsFactors = FALSE
    )

    for (member_index in seq_along(members)) {
      variant_index <- members[[member_index]]
      member_rows[[member_row_index]] <- data.frame(
        locus_id = metadata$locus_id,
        trait = metadata$trait,
        analysis_scope = metadata$analysis_scope,
        branch = metadata$branch,
        cs_index = cs_index,
        cs_id = cs_name,
        susie_effect_index = as.integer(effect_indices[[cs_index]]),
        member_index = member_index,
        matrix_variant_index = as.integer(
          z_table$matrix_variant_index[[variant_index]]
        ),
        chrom = as.character(z_table$chrom[[variant_index]]),
        pos = as.integer(z_table$pos[[variant_index]]),
        ld_variant_id = z_table$ld_variant_id[[variant_index]],
        ld_a1 = z_table$ld_a1[[variant_index]],
        ld_a2 = z_table$ld_a2[[variant_index]],
        pip = as.numeric(fit$pip[[variant_index]]),
        z_for_ld_a1 = as.numeric(z_table$z_for_ld_a1[[variant_index]]),
        p_lrt = as.numeric(z_table$p_lrt[[variant_index]]),
        stringsAsFactors = FALSE
      )
      member_row_index <- member_row_index + 1L
    }
  }

  list(
    summary = do.call(rbind, summary_rows),
    members = do.call(rbind, member_rows)
  )
}

write_failure_status <- function(output_dir, locus_id, branch, reason) {
  if (!dir.exists(output_dir)) {
    return(invisible(NULL))
  }
  status_path <- file.path(output_dir, "run_status.tsv")
  if (!file.exists(status_path)) {
    write_tsv(data.frame(
      locus_id = locus_id,
      branch = branch,
      status = "FAILED",
      reason = reason,
      stringsAsFactors = FALSE
    ), status_path)
  }
  write_tsv(data.frame(
    overall_status = "FAILED",
    reason = reason,
    failed_at = format(Sys.time(), "%F %T"),
    stringsAsFactors = FALSE
  ), file.path(output_dir, "COMPLETED.failed"))
}

main <- function(arguments) {
  locus_id <- arguments[["locus-id"]]
  manifest_path <- arguments[["manifest"]]
  qc_summary_path <- arguments[["qc-summary"]]
  metadata_manifest_path <- arguments[["metadata-manifest"]]
  branch <- arguments[["branch"]]
  output_dir <- arguments[["output-dir"]]

  assert_true(
    branch %in% c("main", "strict_no_ambiguous"),
    paste0("branch must be one of main or strict_no_ambiguous, current value: ", branch)
  )
  check_nonempty_file(manifest_path, "Step 1.9 locus manifest")
  check_nonempty_file(qc_summary_path, "Step 1.9 matrix QC summary")
  check_nonempty_file(metadata_manifest_path, "Step 1.8 locus metadata manifest")
  if (dir.exists(output_dir) || file.exists(output_dir)) {
    stop_message("output directory already exists, refuse to overwrite: ", output_dir)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  log_message("START", "read locus=", locus_id, "，branch=", branch)
  if (!requireNamespace("susieR", quietly = TRUE)) {
    stop_message("current R environment missing susieR")
  }
  installed_version <- as.character(utils::packageVersion("susieR"))
  assert_true(
    identical(installed_version, EXPECTED_SUSIER_VERSION),
    paste0(
      "susieR version does not match the required version: actual=", installed_version,
      "; expected=", EXPECTED_SUSIER_VERSION
    )
  )

  manifest <- read_tsv(manifest_path)
  manifest_columns <- c(
    "locus_id", "trait", "analysis_scope", "locus_status", "sample_count",
    "main_variant_count", "strict_variant_count", "locus_output",
    "main_x_ref", "main_r_signed", "main_z_beta_se",
    "strict_x_ref", "strict_r_signed", "strict_z_beta_se"
  )
  check_columns(manifest, manifest_columns, manifest_path)
  locus_rows <- manifest[manifest$locus_id == locus_id, , drop = FALSE]
  assert_true(nrow(locus_rows) == 1L, paste0(
    "manifest in locus_id must be unique and exists: ", locus_id,
    "; match row count=", nrow(locus_rows)
  ))
  locus_row <- locus_rows[1L, , drop = FALSE]
  assert_true(
    identical(as.character(locus_row$locus_status[[1L]]), "PASS"),
    paste0("Step 1.9 locus_status is not PASS: ", locus_row$locus_status[[1L]])
  )

  metadata_manifest <- read_tsv(metadata_manifest_path)
  metadata_manifest_columns <- c(
    "locus_id", "trait", "analysis_scope", "locus_input_status",
    "locus_metadata"
  )
  check_columns(
    metadata_manifest,
    metadata_manifest_columns,
    metadata_manifest_path
  )
  metadata_manifest_rows <- metadata_manifest[
    metadata_manifest$locus_id == locus_id,
    , drop = FALSE
  ]
  assert_true(nrow(metadata_manifest_rows) == 1L, paste0(
    "Step 1.8 manifest in locus_id must be unique and exists: ", locus_id,
    "; match row count=", nrow(metadata_manifest_rows)
  ))
  metadata_manifest_row <- metadata_manifest_rows[1L, , drop = FALSE]
  assert_true(
    identical(
      as.character(metadata_manifest_row$locus_input_status[[1L]]),
      "READY_FOR_GENOTYPE_MATRIX"
    ),
    paste0(
      "Step 1.8 locus_input_status unavailable in modeling: ",
      metadata_manifest_row$locus_input_status[[1L]]
    )
  )
  assert_true(
    identical(
      as.character(metadata_manifest_row$trait[[1L]]),
      as.character(locus_row$trait[[1L]])
    ) && identical(
      as.character(metadata_manifest_row$analysis_scope[[1L]]),
      as.character(locus_row$analysis_scope[[1L]])
    ),
    "Step 1.8 and Step 1.9 differ in trait or analysis_scope"
  )

  locus_metadata_path <- metadata_manifest_row$locus_metadata[[1L]]
  check_nonempty_file(locus_metadata_path, "Step 1.8 locus_metadata.tsv")
  lead_metadata <- read_key_value_metadata(
    locus_metadata_path,
    c(
      "locus_id", "trait", "analysis_scope", "chrom", "lead_variant_id",
      "lead_pos", "lead_p_lrt", "lead_p_lrt_tier",
      "lead_in_harmonized_variants"
    )
  )
  assert_true(
    identical(unname(lead_metadata[["locus_id"]]), locus_id) &&
      identical(
        unname(lead_metadata[["trait"]]),
        as.character(locus_row$trait[[1L]])
      ) && identical(
        unname(lead_metadata[["analysis_scope"]]),
        as.character(locus_row$analysis_scope[[1L]])
      ),
    "locus_metadata.tsv of locus_id, trait or analysis_scope mismatch"
  )
  gwas_lead_pos <- suppressWarnings(as.integer(lead_metadata[["lead_pos"]]))
  gwas_lead_p_lrt <- suppressWarnings(as.numeric(lead_metadata[["lead_p_lrt"]]))
  assert_true(nzchar(unname(lead_metadata[["chrom"]])),
              "locus_metadata.tsv of chrom is empty")
  assert_true(nzchar(unname(lead_metadata[["lead_variant_id"]])),
              "locus_metadata.tsv of lead_variant_id is empty")
  assert_true(nzchar(unname(lead_metadata[["lead_p_lrt_tier"]])),
              "locus_metadata.tsv of lead_p_lrt_tier is empty")
  assert_true(is.finite(gwas_lead_pos) && gwas_lead_pos >= 1L,
              "locus_metadata.tsv of lead_pos invalid")
  assert_true(
    is.finite(gwas_lead_p_lrt) && gwas_lead_p_lrt >= 0 && gwas_lead_p_lrt <= 1,
    "locus_metadata.tsv of lead_p_lrt must be located in [0,1]"
  )
  assert_true(
    lead_metadata[["lead_in_harmonized_variants"]] %in% c("yes", "no"),
    "locus_metadata.tsv of lead_in_harmonized_variants must be yes/no"
  )

  qc_summary <- read_tsv(qc_summary_path)
  qc_columns <- c(
    "locus_id", "trait", "analysis_scope", "branch", "status",
    "sample_count", "output_variant_count", "n_locus_median",
    "n_obs_min", "n_obs_max"
  )
  check_columns(qc_summary, qc_columns, qc_summary_path)
  qc_rows <- qc_summary[
    qc_summary$locus_id == locus_id & qc_summary$branch == branch,
    , drop = FALSE
  ]
  assert_true(nrow(qc_rows) == 1L, paste0(
    "matrix QC in locus_id+branch must be unique and exists: ", locus_id, "/", branch,
    "; match row count=", nrow(qc_rows)
  ))
  qc_row <- qc_rows[1L, , drop = FALSE]
  assert_true(
    identical(as.character(qc_row$status[[1L]]), "PASS"),
    paste0("Step 1.9 matrix QC status is not PASS: ", qc_row$status[[1L]])
  )

  if (branch == "main") {
    x_path <- locus_row$main_x_ref[[1L]]
    r_path <- locus_row$main_r_signed[[1L]]
    z_path <- locus_row$main_z_beta_se[[1L]]
    manifest_variant_count <- as.integer(locus_row$main_variant_count[[1L]])
  } else {
    x_path <- locus_row$strict_x_ref[[1L]]
    r_path <- locus_row$strict_r_signed[[1L]]
    z_path <- locus_row$strict_z_beta_se[[1L]]
    manifest_variant_count <- as.integer(locus_row$strict_variant_count[[1L]])
  }
  check_nonempty_file(x_path, paste0(branch, " X_ref.rds"))
  check_nonempty_file(r_path, paste0(branch, " R_signed.rds"))
  check_nonempty_file(z_path, paste0(branch, " z_beta_se.tsv.gz"))

  log_message("INPUT", "read X, R and already harmonize z-score")
  x_reference <- readRDS(x_path)
  r_signed <- readRDS(r_path)
  z_table <- read_tsv(z_path)
  z_columns <- c(
    "matrix_variant_index", "locus_id", "trait", "analysis_scope", "chrom",
    "pos", "ld_variant_id", "ld_a1", "ld_a2", "beta_for_gwas_allele1",
    "beta_for_ld_a1", "se", "z_for_ld_a1", "p_lrt", "n_obs",
    "ambiguous_at_cg", "harmonization_status"
  )
  check_columns(z_table, z_columns, z_path)

  variant_count <- nrow(z_table)
  assert_true(is.matrix(x_reference) && is.numeric(x_reference),
              "X_ref.rds must be numeric value matrix")
  assert_true(is.matrix(r_signed) && is.numeric(r_signed),
              "R_signed.rds must be numeric value matrix")
  assert_true(variant_count >= 2L, "current locus less than 2 items eligible for modeling variant")
  assert_true(ncol(x_reference) == variant_count,
              "X column count and z table row counts differ")
  assert_true(nrow(r_signed) == variant_count && ncol(r_signed) == variant_count,
              "R dimensions and z table row counts differ")
  assert_true(variant_count == manifest_variant_count,
              "Actual variant count differs from the locus manifest")
  assert_true(variant_count == as.integer(qc_row$output_variant_count[[1L]]),
              "Actual variant count differs from the matrix QC record")
  assert_true(nrow(x_reference) == as.integer(locus_row$sample_count[[1L]]),
              "X sample count differs from the locus manifest")
  assert_true(nrow(x_reference) == as.integer(qc_row$sample_count[[1L]]),
              "X sample count differs from the matrix QC record")
  assert_true(all(is.finite(x_reference)), "X contains non-finite values")
  assert_true(all(is.finite(r_signed)), "R contains non-finite values")
  assert_true(all(is.finite(as.numeric(z_table$z_for_ld_a1))),
              "z_for_ld_a1 contains non-finite values")
  assert_true(all(is.finite(as.numeric(z_table$se))) && all(z_table$se > 0),
              "se must be all are finite positive number")
  assert_true(all(is.finite(as.numeric(z_table$n_obs))) && all(z_table$n_obs > 1),
              "n_obs must be all greater than 1")
  assert_true(all(is.finite(as.numeric(z_table$p_lrt))) &&
                all(z_table$p_lrt >= 0) && all(z_table$p_lrt <= 1),
              "p_lrt must be all located in [0,1]")

  expected_indices <- seq_len(variant_count)
  assert_true(identical(as.integer(z_table$matrix_variant_index), expected_indices),
              "matrix_variant_index not from 1 start of continuous order")
  assert_true(!anyDuplicated(z_table$ld_variant_id), "ld_variant_id contains duplicate")
  assert_true(identical(as.character(z_table$ld_variant_id), colnames(x_reference)),
              "Variant order differs between the z table and X")
  assert_true(identical(colnames(x_reference), rownames(r_signed)) &&
                identical(colnames(x_reference), colnames(r_signed)),
              "Variant order differs between X and R")
  assert_true(all(z_table$locus_id == locus_id), "z table unexpected inclusion of other locus_id")
  assert_true(all(z_table$trait == locus_row$trait[[1L]]), "z table trait mismatch")
  assert_true(all(z_table$analysis_scope == locus_row$analysis_scope[[1L]]),
              "z table analysis_scope mismatch")
  assert_true(
    all(as.character(z_table$chrom) == unname(lead_metadata[["chrom"]])),
    "Lead chromosome in locus_metadata.tsv differs from the z table"
  )
  assert_true(all(z_table$ambiguous_at_cg %in% c("yes", "no")),
              "ambiguous_at_cg unsupported values encountered")
  supported_harmonization <- c(
    "HARMONIZED_NONAMBIGUOUS",
    "HARMONIZED_AMBIGUOUS_FREQUENCY_RESOLVED"
  )
  assert_true(all(z_table$harmonization_status %in% supported_harmonization),
              "harmonization_status values unsupported by the current input were encountered")
  if (branch == "strict_no_ambiguous") {
    assert_true(all(z_table$ambiguous_at_cg == "no"),
                "strict_no_ambiguous branch still contains ambiguous A/T or C/G variant")
  }

  reconstructed_z <- as.numeric(z_table$beta_for_ld_a1) / as.numeric(z_table$se)
  z_difference <- max(abs(reconstructed_z - as.numeric(z_table$z_for_ld_a1)))
  assert_true(is.finite(z_difference) && z_difference <= MATRIX_TOLERANCE,
              paste0("beta_for_ld_a1/se differs from z_for_ld_a1: ", z_difference))
  symmetry_error <- max(abs(r_signed - t(r_signed)))
  diagonal_error <- max(abs(diag(r_signed) - 1))
  assert_true(symmetry_error <= MATRIX_TOLERANCE,
              paste0("R asymmetric: ", symmetry_error))
  assert_true(diagonal_error <= MATRIX_TOLERANCE,
              paste0("R diagonal deviation 1: ", diagonal_error))

  n_locus <- as.numeric(qc_row$n_locus_median[[1L]])
  n_from_z <- stats::median(as.numeric(z_table$n_obs))
  assert_true(is.finite(n_locus) && n_locus > 1, "n_locus_median invalid")
  assert_true(abs(n_locus - n_from_z) <= MATRIX_TOLERANCE,
              "matrix QC of n_locus_median with z table medians differ")

  metadata <- list(
    locus_id = locus_id,
    trait = as.character(locus_row$trait[[1L]]),
    analysis_scope = as.character(locus_row$analysis_scope[[1L]]),
    branch = branch
  )
  z_vector <- as.numeric(z_table$z_for_ld_a1)
  L_used <- min(10L, variant_count)
  warning_store <- new.env(parent = emptyenv())
  warning_store$stage <- character()
  warning_store$message <- character()

  log_message("DIAGNOSTIC", "compute estimate_s_rss with kriging_rss")
  estimate_s <- capture_stage_warnings(
    "estimate_s_rss",
    susieR::estimate_s_rss(z_vector, r_signed, n = n_locus),
    warning_store
  )
  kriging <- capture_stage_warnings(
    "kriging_rss",
    susieR::kriging_rss(
      z_vector,
      r_signed,
      n = n_locus,
      s = estimate_s
    ),
    warning_store
  )
  kriging_table <- as.data.frame(kriging$conditional_dist)
  assert_true(nrow(kriging_table) == variant_count,
              "Row count returned by kriging_rss differs from the variant count")
  kriging_table <- cbind(
    z_table[c("matrix_variant_index", "ld_variant_id", "chrom", "pos")],
    kriging_table
  )

  log_message(
    "MODEL",
    "run susie_rss: p=", variant_count,
    "，n_locus=", format(n_locus, scientific = FALSE),
    "，B=", nrow(x_reference), "，L=", L_used
  )
  fit <- capture_stage_warnings(
    "susie_rss",
    susieR::susie_rss(
      z = z_vector,
      X = x_reference,
      n = n_locus,
      z_method = "wald",
      L = L_used,
      coverage = SUSIE_COVERAGE,
      min_abs_corr = SUSIE_MIN_ABS_CORR,
      standardize = TRUE,
      estimate_residual_variance = FALSE,
      refine = TRUE,
      R_finite = TRUE,
      R_mismatch = "eb",
      max_iter = SUSIE_MAX_ITER,
      check_input = TRUE,
      verbose = FALSE
    ),
    warning_store
  )

  assert_true(!is.null(fit$pip) && length(fit$pip) == variant_count,
              "susie_rss did not return complete PIP")
  assert_true(all(is.finite(fit$pip)) && all(fit$pip >= 0) && all(fit$pip <= 1),
              "PIP contains non-finite values or exceeds [0,1]")

  cs_list <- fit$sets$cs
  cs_names <- credible_set_names(cs_list)
  cs_labels <- credible_set_labels(cs_list, cs_names, variant_count)
  variant_pip <- z_table
  variant_pip$pip <- as.numeric(fit$pip)
  variant_pip$pip_rank <- rank(-variant_pip$pip, ties.method = "min")
  variant_pip$credible_sets <- cs_labels
  variant_pip$in_credible_set <- ifelse(nzchar(cs_labels), "yes", "no")
  variant_pip <- variant_pip[order(-variant_pip$pip, variant_pip$matrix_variant_index), ]

  cs_tables <- extract_credible_sets(fit, z_table, metadata)
  diagnostics <- fit$R_finite_diagnostics
  if (is.null(diagnostics)) {
    diagnostics <- list()
  }
  lambda_summary <- numeric_summary(diagnostics$lambda_bias)
  corrected_summary <- numeric_summary(diagnostics$B_corrected)
  penalty_summary <- numeric_summary(diagnostics$per_variable_penalty)
  rhat_summary <- numeric_summary(diagnostics$Rhat_diag_deviation)
  reliability_flag <- isTRUE(diag_scalar(
    diagnostics, "R_reliability_flag", FALSE
  ))
  artifact_flag <- isTRUE(diag_scalar(diagnostics, "artifact_flag", FALSE))
  sensitivity_flag <- isTRUE(diag_scalar(
    diagnostics, "R_sensitivity_flag", FALSE
  ))
  kriging_outlier_count <- sum(
    kriging_table$logLR > 2 & abs(kriging_table$z) > 2,
    na.rm = TRUE
  )
  converged <- isTRUE(fit$converged)
  warning_rows <- warning_table(warning_store)
  model_status <- if (converged) "CONVERGED" else "FAILED_NONCONVERGED"
  diagnostic_warning <- reliability_flag || sensitivity_flag ||
    artifact_flag || kriging_outlier_count > 0L
  overall_status <- if (!converged) {
    "FAILED"
  } else if (diagnostic_warning || nrow(warning_rows) > 0L) {
    "PASS_WITH_WARNINGS"
  } else {
    "PASS"
  }

  rss_diagnostics <- data.frame(
    locus_id = locus_id,
    trait = metadata$trait,
    analysis_scope = metadata$analysis_scope,
    branch = branch,
    susieR_version = installed_version,
    variant_count = variant_count,
    reference_sample_count_B = nrow(x_reference),
    n_locus_median = n_locus,
    n_obs_min = min(as.numeric(z_table$n_obs)),
    n_obs_max = max(as.numeric(z_table$n_obs)),
    estimate_s_rss = as.numeric(estimate_s),
    kriging_outlier_count_logLR_gt2_absz_gt2 = kriging_outlier_count,
    R_finite_B = as.numeric(diag_scalar(diagnostics, "B")),
    R_finite_p = as.numeric(diag_scalar(diagnostics, "p")),
    effective_rank = as.numeric(diag_scalar(diagnostics, "effective_rank")),
    r_over_B = as.numeric(diag_scalar(diagnostics, "r_over_B")),
    Rhat_diag_deviation_max = unname(rhat_summary[["max"]]),
    Q_art = as.numeric(diag_scalar(diagnostics, "Q_art")),
    artifact_flag = artifact_flag,
    R_sensitivity_flag = sensitivity_flag,
    R_reliability_flag = reliability_flag,
    low_eigen_count = as.numeric(diag_scalar(diagnostics, "low_eigen_count")),
    low_eigen_fraction = as.numeric(diag_scalar(diagnostics, "low_eigen_fraction")),
    eig_delta = as.numeric(diag_scalar(diagnostics, "eig_delta")),
    mode_label = as.character(diag_scalar(diagnostics, "mode_label", NA_character_)),
    lambda_bias_min = unname(lambda_summary[["min"]]),
    lambda_bias_mean = unname(lambda_summary[["mean"]]),
    lambda_bias_max = unname(lambda_summary[["max"]]),
    B_corrected_min = unname(corrected_summary[["min"]]),
    B_corrected_mean = unname(corrected_summary[["mean"]]),
    B_corrected_max = unname(corrected_summary[["max"]]),
    per_variable_penalty_max = unname(penalty_summary[["max"]]),
    R_symmetry_max_abs = symmetry_error,
    R_diagonal_max_abs_difference = diagonal_error,
    z_beta_se_max_abs_difference = z_difference,
    warning_count = nrow(warning_rows),
    converged = converged,
    niter = as.integer(if (is.null(fit$niter)) NA_integer_ else fit$niter),
    final_elbo = as.numeric(if (is.null(fit$elbo)) NA_real_ else tail(fit$elbo, 1L)),
    credible_set_count = length(cs_names),
    stringsAsFactors = FALSE
  )

  top_index <- which.max(fit$pip)
  locus_summary <- data.frame(
    locus_id = locus_id,
    trait = metadata$trait,
    analysis_scope = metadata$analysis_scope,
    branch = branch,
    status = overall_status,
    model_status = model_status,
    gwas_lead_variant_id = unname(lead_metadata[["lead_variant_id"]]),
    gwas_lead_pos = gwas_lead_pos,
    gwas_lead_p_lrt = gwas_lead_p_lrt,
    gwas_lead_p_lrt_tier = unname(lead_metadata[["lead_p_lrt_tier"]]),
    gwas_lead_in_harmonized_variants = unname(
      lead_metadata[["lead_in_harmonized_variants"]]
    ),
    variant_count = variant_count,
    reference_sample_count = nrow(x_reference),
    n_locus_median = n_locus,
    min_observed_p_lrt = min(as.numeric(z_table$p_lrt)),
    converged = converged,
    niter = as.integer(if (is.null(fit$niter)) NA_integer_ else fit$niter),
    L = L_used,
    coverage = SUSIE_COVERAGE,
    estimate_s_rss = as.numeric(estimate_s),
    credible_set_count = length(cs_names),
    max_pip = max(as.numeric(fit$pip)),
    top_variant_id = z_table$ld_variant_id[[top_index]],
    top_variant_pos = as.integer(z_table$pos[[top_index]]),
    artifact_flag = artifact_flag,
    R_sensitivity_flag = sensitivity_flag,
    R_reliability_flag = reliability_flag,
    kriging_outlier_count = kriging_outlier_count,
    output_dir = output_dir,
    variant_pip = file.path(output_dir, "variant_pip.tsv.gz"),
    credible_sets = file.path(output_dir, "credible_sets.tsv"),
    credible_set_members = file.path(output_dir, "credible_set_members.tsv.gz"),
    rss_diagnostics = file.path(output_dir, "rss_diagnostics.tsv"),
    kriging_diagnostics = file.path(output_dir, "kriging_rss.tsv.gz"),
    susie_fit = file.path(output_dir, "susie_fit.rds"),
    stringsAsFactors = FALSE
  )

  run_status <- data.frame(
    locus_id = locus_id,
    trait = metadata$trait,
    analysis_scope = metadata$analysis_scope,
    branch = branch,
    status = if (converged) "SUCCESS" else "FAILED",
    overall_status = overall_status,
    model_status = model_status,
    gwas_lead_variant_id = unname(lead_metadata[["lead_variant_id"]]),
    gwas_lead_pos = gwas_lead_pos,
    gwas_lead_p_lrt = gwas_lead_p_lrt,
    gwas_lead_p_lrt_tier = unname(lead_metadata[["lead_p_lrt_tier"]]),
    gwas_lead_in_harmonized_variants = unname(
      lead_metadata[["lead_in_harmonized_variants"]]
    ),
    susieR_version = installed_version,
    z_method = "wald",
    L = L_used,
    coverage = SUSIE_COVERAGE,
    min_abs_corr = SUSIE_MIN_ABS_CORR,
    estimate_residual_variance = FALSE,
    refine = TRUE,
    R_finite = TRUE,
    R_mismatch = "eb",
    max_iter = SUSIE_MAX_ITER,
    check_input = TRUE,
    variant_count = variant_count,
    n_locus_median = n_locus,
    reference_sample_count_B = nrow(x_reference),
    warning_count = nrow(warning_rows),
    artifact_flag = artifact_flag,
    R_sensitivity_flag = sensitivity_flag,
    R_reliability_flag = reliability_flag,
    kriging_outlier_count = kriging_outlier_count,
    completed_at = format(Sys.time(), "%F %T"),
    stringsAsFactors = FALSE
  )

  log_message("OUTPUT", "write PIP, credible set, diagnostic and model object")
  saveRDS(fit, file.path(output_dir, "susie_fit.rds"), compress = "gzip")
  write_tsv_gz(variant_pip, file.path(output_dir, "variant_pip.tsv.gz"))
  write_tsv(cs_tables$summary, file.path(output_dir, "credible_sets.tsv"))
  write_tsv_gz(
    cs_tables$members,
    file.path(output_dir, "credible_set_members.tsv.gz")
  )
  write_tsv_gz(kriging_table, file.path(output_dir, "kriging_rss.tsv.gz"))
  write_tsv(rss_diagnostics, file.path(output_dir, "rss_diagnostics.tsv"))
  write_tsv(warning_rows, file.path(output_dir, "warnings.tsv"))
  write_tsv(locus_summary, file.path(output_dir, "locus_summary.tsv"))
  write_tsv(run_status, file.path(output_dir, "run_status.tsv"))

  if (!converged) {
    write_tsv(data.frame(
      overall_status = "FAILED_NONCONVERGED",
      locus_id = locus_id,
      branch = branch,
      completed_at = format(Sys.time(), "%F %T"),
      stringsAsFactors = FALSE
    ), file.path(output_dir, "COMPLETED.failed"))
    stop_message("susie_rss not converged, result already retain for diagnostic")
  }

  write_tsv(data.frame(
    overall_status = overall_status,
    locus_id = locus_id,
    branch = branch,
    completed_at = format(Sys.time(), "%F %T"),
    stringsAsFactors = FALSE
  ), file.path(output_dir, "COMPLETED.ok"))
  log_message(
    "DONE",
    "locus completed / done: ", locus_id,
    "; status=", overall_status,
    "；credible_sets=", length(cs_names)
  )
}

arguments <- NULL
exit_code <- tryCatch({
  arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
  main(arguments)
  0L
}, error = function(error_condition) {
  log_message("ERROR", conditionMessage(error_condition))
  if (!is.null(arguments)) {
    write_failure_status(
      arguments[["output-dir"]],
      arguments[["locus-id"]],
      arguments[["branch"]],
      conditionMessage(error_condition)
    )
  }
  1L
})

quit(save = "no", status = exit_code, runLast = FALSE)
