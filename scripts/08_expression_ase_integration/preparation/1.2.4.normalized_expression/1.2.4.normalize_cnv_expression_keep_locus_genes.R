#!/path/to/software/miniconda3/bin/Rscript


options(stringsAsFactors = FALSE)

log_msg <- function(stage, message) {
  cat(sprintf("[%s ] [%s] %s\n", format(Sys.time(), "%F %T"), stage, message))
  flush.console()
}

project_dir <- "/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation"
script_dir <- file.path(project_dir, "scripts", "1.2.4.normalized_expression")
input_123_run_id <- "20260821_141015"
input_122_run_id <- "20260821_134953"

run_id <- Sys.getenv("RUN_ID", unset = format(Sys.time(), "%Y%m%d_%H%M%S"))
if (!grepl("^[0-9]{8}_[0-9]{6}$", run_id)) {
  stop("RUN_ID must use YYYYMMDD_HHMMSS format, found: ", run_id)
}


expected_samples <- 19L
expected_genes <- 40825L
expected_gene_LOC_rows <- 18858L
expected_dog5_new_gene_rows <- 3828L
expected_anno_cols <- c(
  "gene_id", "gene_name", "seqname", "start", "end", "strand", "midpoint", "source_feature"
)
expected_sample_meta_cols <- c(
  "sample_id", "animal_id", "term", "raw_data_name", "counts_file", "bam_file", "bai_file"
)

in_matrix <- file.path(
  project_dir,
  "1.2.3.cnv_raw_counts",
  input_123_run_id,
  "bulk_rna.star.reverse_stranded.cnv_raw_counts.manual_new_gene.tsv"
)

in_gene_annotation <- file.path(
  project_dir,
  "1.2.3.cnv_raw_counts",
  input_123_run_id,
  "bulk_rna.star.reverse_stranded.cnv_raw_counts.manual_new_gene.gene_annotation.tsv"
)

input_123_completion_marker <- file.path(
  project_dir,
  "1.2.3.cnv_raw_counts",
  input_123_run_id,
  "execution_completion.marker"
)

sample_meta <- file.path(
  project_dir,
  "1.2.2.expression_matrix",
  input_122_run_id,
  "bulk_rna.star.samples.tsv"
)

input_122_completion_marker <- file.path(
  project_dir,
  "1.2.2.expression_matrix",
  input_122_run_id,
  "execution_completion.marker"
)

output_root <- file.path(project_dir, "1.2.4.normalized_expression")
out_dir <- file.path(output_root, run_id)
log_dir <- file.path(script_dir, "log")
master_log <- file.path(log_dir, paste0("1.2.4.normalized_expression.", run_id, ".log"))
completion_marker <- file.path(out_dir, "execution_completion.marker")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
if (dir.exists(out_dir)) {
  stop("Timestamped output directory already exists; refusing to overwrite: ", out_dir)
}
if (file.exists(master_log)) {
  stop("Run-specific master log already exists; refusing to overwrite: ", master_log)
}
if (!dir.create(out_dir, recursive = FALSE, showWarnings = FALSE)) {
  stop("Failed to create timestamped output directory: ", out_dir)
}

log_connection <- file(master_log, open = "wt")
sink(log_connection, type = "output")
sink(log_connection, type = "message")

out_sample_qc <- file.path(out_dir, "bulk_rna.cnv_expression.manual_new_gene.sample_qc.tsv")
out_gene_qc <- file.path(out_dir, "bulk_rna.cnv_expression.manual_new_gene.gene_qc.tsv")
out_norm_counts <- file.path(out_dir, "bulk_rna.cnv_expression.manual_new_gene.sizefactor_normalized_counts.tsv")
out_log2_norm <- file.path(out_dir, "bulk_rna.cnv_expression.manual_new_gene.log2_sizefactor_norm_counts.tsv")
out_log2_norm_filtered <- file.path(out_dir, "bulk_rna.cnv_expression.manual_new_gene.log2_sizefactor_norm_counts.cpm1_min2.tsv")
out_gene_annotation_all <- file.path(out_dir, "bulk_rna.cnv_expression.manual_new_gene.gene_annotation.tsv")
out_gene_annotation_filtered <- file.path(out_dir, "bulk_rna.cnv_expression.manual_new_gene.gene_annotation.cpm1_min2.tsv")
out_check_log <- file.path(out_dir, "bulk_rna.cnv_expression.manual_new_gene.normalization.check.log")

main <- function() {
log_msg("start", "start finemaping CNV expression normalization")
log_msg("config", paste0("run_id=", run_id))
log_msg("config", paste0("project_dir=", project_dir))
log_msg("config", paste0("input_matrix=", in_matrix))
log_msg("config", paste0("input_gene_annotation=", in_gene_annotation))
log_msg("config", paste0("sample_meta=", sample_meta))
log_msg("config", paste0("out_dir=", out_dir))

log_msg("check-inputs", "check input file")
for (marker in c(input_123_completion_marker, input_122_completion_marker)) {
  if (!file.exists(marker) || file.info(marker)$size == 0) {
    stop("Upstream completion marker not found or empty: ", marker)
  }
  marker_lines <- readLines(marker, warn = FALSE)
  if (!any(marker_lines == "execution_status\tsuccess")) {
    stop("Upstream completion marker does not report success: ", marker)
  }
}

if (!file.exists(in_matrix) || file.info(in_matrix)$size == 0) {
  stop("Input CNV raw-count matrix not found or empty: ", in_matrix)
}

if (!file.exists(in_gene_annotation) || file.info(in_gene_annotation)$size == 0) {
  stop("Input gene annotation not found or empty: ", in_gene_annotation)
}

if (!file.exists(sample_meta) || file.info(sample_meta)$size == 0) {
  stop("Input sample metadata not found or empty: ", sample_meta)
}

log_msg("load-inputs", "read count matrix, gene annotation and sample metadata")
counts_df <- read.delim(in_matrix, check.names = FALSE)
gene_anno <- read.delim(in_gene_annotation, check.names = FALSE)
meta <- read.delim(sample_meta, check.names = FALSE)

log_msg("validate-inputs", "check matrix and annotation table structure")
if (ncol(counts_df) < 2) {
  stop("Input matrix must contain gene_id plus at least one sample column.")
}

if (colnames(counts_df)[1] != "gene_id") {
  stop("First column of count matrix must be gene_id, found: ", colnames(counts_df)[1])
}

if (ncol(counts_df) - 1L != expected_samples) {
  stop("Unexpected sample count in count matrix: ", ncol(counts_df) - 1L, "; expected ", expected_samples)
}

if (nrow(counts_df) != expected_genes) {
  stop("Unexpected gene count in count matrix: ", nrow(counts_df), "; expected ", expected_genes)
}

if (!identical(colnames(gene_anno), expected_anno_cols)) {
  stop(
    "Unexpected gene annotation columns: ", paste(colnames(gene_anno), collapse = ","),
    "; expected ", paste(expected_anno_cols, collapse = ",")
  )
}

if (nrow(gene_anno) != expected_genes) {
  stop("Unexpected gene count in gene annotation: ", nrow(gene_anno), "; expected ", expected_genes)
}

if (!identical(colnames(meta), expected_sample_meta_cols)) {
  stop(
    "Unexpected sample metadata columns: ", paste(colnames(meta), collapse = ","),
    "; expected ", paste(expected_sample_meta_cols, collapse = ",")
  )
}

if (nrow(meta) != expected_samples) {
  stop("Unexpected sample metadata rows: ", nrow(meta), "; expected ", expected_samples)
}

if (anyDuplicated(meta$sample_id) > 0) {
  dup <- meta$sample_id[duplicated(meta$sample_id)][1]
  stop("Duplicated sample_id found in sample metadata: ", dup)
}

if (anyDuplicated(counts_df$gene_id) > 0) {
  dup <- counts_df$gene_id[duplicated(counts_df$gene_id)][1]
  stop("Duplicated gene_id found in count matrix: ", dup)
}

if (anyDuplicated(gene_anno$gene_id) > 0) {
  dup <- gene_anno$gene_id[duplicated(gene_anno$gene_id)][1]
  stop("Duplicated gene_id found in gene annotation: ", dup)
}

if (!identical(counts_df$gene_id, gene_anno$gene_id)) {
  stop("Count matrix and gene annotation are not in the same gene order.")
}

if (any(counts_df$gene_id == "MissingGeneID")) {
  stop("MissingGeneID is present in count matrix. It has no reliable genomic coordinate.")
}

if (sum(grepl("^gene-LOC", counts_df$gene_id)) != expected_gene_LOC_rows) {
  stop("Unexpected gene-LOC row count in count matrix; expected ", expected_gene_LOC_rows)
}

if (sum(grepl("^DOG5_NEWGENE_", counts_df$gene_id)) != expected_dog5_new_gene_rows) {
  stop("Unexpected DOG5_NEWGENE row count in count matrix; expected ", expected_dog5_new_gene_rows)
}

if (any(grepl("^UNIFIED_NEWGENE_", counts_df$gene_id))) {
  stop("Legacy UNIFIED_NEWGENE rows are present in the five-dog-only count matrix.")
}

gene_id <- counts_df$gene_id
counts <- as.matrix(counts_df[, -1, drop = FALSE])
suppressWarnings(storage.mode(counts) <- "numeric")
rownames(counts) <- gene_id

if (any(is.na(counts))) {
  stop("NA or non-numeric values found in count matrix.")
}

if (any(counts < 0)) {
  stop("Negative counts found in count matrix.")
}

if (any(counts != floor(counts))) {
  stop("Non-integer count values found in raw count matrix.")
}

log_msg(
  "validate-inputs",
  paste0("samples=", ncol(counts), "; genes=", nrow(counts))
)

sample_id <- colnames(counts)
if (!identical(sample_id, meta$sample_id)) {
  stop("Count matrix sample columns and sample metadata are not in the same sample order.")
}

lib_size <- colSums(counts)

if (any(lib_size <= 0)) {
  bad <- sample_id[lib_size <= 0]
  stop("Samples with zero library size found: ", paste(bad, collapse = ","))
}

detected_genes_count_gt0 <- colSums(counts > 0)
detected_genes_count_ge10 <- colSums(counts >= 10)

log_msg("normalize", "compute CPM, size factors and log2 normalized expression")
cpm <- t(t(counts) / lib_size * 1e6)
detected_genes_cpm_ge1 <- colSums(cpm >= 1)


positive_counts <- counts[rowSums(counts > 0) == ncol(counts), , drop = FALSE]
if (nrow(positive_counts) < 100) {
  stop(
    "Too few genes have positive counts in all samples for median-ratio size factor estimation: ",
    nrow(positive_counts)
  )
}

geomeans <- exp(rowMeans(log(positive_counts)))
ratios <- sweep(positive_counts, 1, geomeans, "/")
size_factors <- apply(ratios, 2, median, na.rm = TRUE)

if (any(!is.finite(size_factors)) || any(size_factors <= 0)) {
  stop("Invalid size factors detected.")
}

norm_counts <- sweep(counts, 2, size_factors, "/")
log2_norm <- log2(norm_counts + 1)

keep_cpm1_min2 <- rowSums(cpm >= 1) >= 2
log2_norm_filtered <- log2_norm[keep_cpm1_min2, , drop = FALSE]
gene_anno_filtered <- gene_anno[keep_cpm1_min2, , drop = FALSE]

sample_qc <- data.frame(
  sample_id = sample_id,
  library_size = as.numeric(lib_size),
  size_factor = as.numeric(size_factors[sample_id]),
  detected_genes_count_gt0 = as.numeric(detected_genes_count_gt0),
  detected_genes_count_ge10 = as.numeric(detected_genes_count_ge10),
  detected_genes_cpm_ge1 = as.numeric(detected_genes_cpm_ge1),
  stringsAsFactors = FALSE
)

sample_qc <- cbind(
  meta,
  sample_qc[, setdiff(colnames(sample_qc), "sample_id"), drop = FALSE]
)

gene_qc <- data.frame(
  gene_id = gene_id,
  gene_name = if ("gene_name" %in% colnames(gene_anno)) gene_anno$gene_name else gene_id,
  seqname = gene_anno$seqname,
  start = gene_anno$start,
  end = gene_anno$end,
  strand = gene_anno$strand,
  midpoint = gene_anno$midpoint,
  is_gene_LOC = grepl("^gene-LOC", gene_id),
  is_dog5_new_gene = grepl("^DOG5_NEWGENE_", gene_id),
  total_count = as.numeric(rowSums(counts)),
  mean_count = as.numeric(rowMeans(counts)),
  samples_count_gt0 = as.numeric(rowSums(counts > 0)),
  samples_count_ge10 = as.numeric(rowSums(counts >= 10)),
  samples_cpm_ge1 = as.numeric(rowSums(cpm >= 1)),
  keep_cpm1_min2 = keep_cpm1_min2,
  stringsAsFactors = FALSE
)

log_msg("write-outputs", "write sample QC, gene QC and normalization matrix")
write.table(sample_qc, out_sample_qc, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_qc, out_gene_qc, sep = "\t", quote = FALSE, row.names = FALSE)

write.table(
  data.frame(gene_id = rownames(norm_counts), norm_counts, check.names = FALSE),
  out_norm_counts,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(gene_id = rownames(log2_norm), log2_norm, check.names = FALSE),
  out_log2_norm,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(gene_id = rownames(log2_norm_filtered), log2_norm_filtered, check.names = FALSE),
  out_log2_norm_filtered,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(gene_anno, out_gene_annotation_all, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_anno_filtered, out_gene_annotation_filtered, sep = "\t", quote = FALSE, row.names = FALSE)

log_lines <- c(
  paste0("status=OK"),
  paste0("run_id=", run_id),
  paste0("input_matrix=", in_matrix),
  paste0("input_gene_annotation=", in_gene_annotation),
  paste0("sample_meta=", sample_meta),
  paste0("samples=", ncol(counts)),
  paste0("input_genes=", nrow(counts)),
  paste0("gene_LOC_rows_kept=", sum(grepl("^gene-LOC", gene_id))),
  paste0("DOG5_NEWGENE_rows=", sum(grepl("^DOG5_NEWGENE_", gene_id))),
  paste0("MissingGeneID_rows=", sum(gene_id == "MissingGeneID")),
  paste0("positive_genes_for_size_factors=", nrow(positive_counts)),
  paste0("output_genes_all=", nrow(log2_norm)),
  paste0("output_genes_cpm1_min2=", nrow(log2_norm_filtered)),
  paste0("gene_LOC_rows_cpm1_min2=", sum(grepl("^gene-LOC", rownames(log2_norm_filtered)))),
  paste0("DOG5_NEWGENE_rows_cpm1_min2=", sum(grepl("^DOG5_NEWGENE_", rownames(log2_norm_filtered)))),
  paste0("min_library_size=", min(lib_size)),
  paste0("max_library_size=", max(lib_size)),
  paste0("min_size_factor=", signif(min(size_factors), 6)),
  paste0("max_size_factor=", signif(max(size_factors), 6)),
  paste0("sample_qc=", out_sample_qc),
  paste0("gene_qc=", out_gene_qc),
  paste0("normalized_counts=", out_norm_counts),
  paste0("log2_normalized_matrix=", out_log2_norm),
  paste0("log2_normalized_matrix_filtered=", out_log2_norm_filtered),
  paste0("gene_annotation_all=", out_gene_annotation_all),
  paste0("gene_annotation_filtered=", out_gene_annotation_filtered)
)

log_msg("write-outputs", paste0("write check log: ", out_check_log))
writeLines(log_lines, out_check_log)
writeLines(log_lines)
writeLines(
  c(
    paste("run_id", run_id, sep = "\t"),
    "execution_status\tsuccess",
    paste("completed_at", format(Sys.time(), "%F %T"), sep = "\t")
  ),
  completion_marker
)
log_msg("success", "execution_status=success")
log_msg("success", "finemaping CNV expression normalization script execution completed; task completed / done")
log_msg("reminder", "execution completion does not establish output completeness, format, counts or biological quality have passed review check, still requires separately review result.")
}

exit_status <- tryCatch(
  {
    main()
    0L
  },
  error = function(e) {
    log_msg("error", conditionMessage(e))
    log_msg("error", "execution_status=failed")
    1L
  }
)

sink(type = "message")
sink(type = "output")
close(log_connection)

if (exit_status != 0L) {
  quit(save = "no", status = exit_status)
}
