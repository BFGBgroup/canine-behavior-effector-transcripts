#!/usr/bin/env Rscript


options(stringsAsFactors = FALSE)

project_dir <- "/path/to/project/legacy/five-dogs/5.cnv"

in_matrix <- file.path(
  project_dir,
  "3.cnv_raw_counts",
  "bulk_rna.star.reverse_stranded.cnv_raw_counts.keep_locus_genes.tsv"
)

in_gene_annotation <- file.path(
  project_dir,
  "3.cnv_raw_counts",
  "bulk_rna.star.reverse_stranded.cnv_raw_counts.keep_locus_genes.gene_annotation.tsv"
)

sample_meta <- file.path(
  project_dir,
  "2.expression_matrix",
  "bulk_rna.star.samples.tsv"
)

out_dir <- file.path(project_dir, "4.normalized_expression")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_sample_qc <- file.path(out_dir, "bulk_rna.cnv_expression.keep_locus_genes.sample_qc.tsv")
out_gene_qc <- file.path(out_dir, "bulk_rna.cnv_expression.keep_locus_genes.gene_qc.tsv")
out_norm_counts <- file.path(out_dir, "bulk_rna.cnv_expression.keep_locus_genes.sizefactor_normalized_counts.tsv")
out_log2_norm <- file.path(out_dir, "bulk_rna.cnv_expression.keep_locus_genes.log2_sizefactor_norm_counts.tsv")
out_log2_norm_filtered <- file.path(out_dir, "bulk_rna.cnv_expression.keep_locus_genes.log2_sizefactor_norm_counts.cpm1_min2.tsv")
out_gene_annotation_all <- file.path(out_dir, "bulk_rna.cnv_expression.keep_locus_genes.gene_annotation.tsv")
out_gene_annotation_filtered <- file.path(out_dir, "bulk_rna.cnv_expression.keep_locus_genes.gene_annotation.cpm1_min2.tsv")
out_check_log <- file.path(out_dir, "bulk_rna.cnv_expression.keep_locus_genes.normalization.check.log")

if (!file.exists(in_matrix) || file.info(in_matrix)$size == 0) {
  stop("Input CNV raw-count matrix not found or empty: ", in_matrix)
}

if (!file.exists(in_gene_annotation) || file.info(in_gene_annotation)$size == 0) {
  stop("Input gene annotation not found or empty: ", in_gene_annotation)
}

counts_df <- read.delim(in_matrix, check.names = FALSE)
gene_anno <- read.delim(in_gene_annotation, check.names = FALSE)

if (ncol(counts_df) < 2) {
  stop("Input matrix must contain gene_id plus at least one sample column.")
}

if (colnames(counts_df)[1] != "gene_id") {
  stop("First column of count matrix must be gene_id, found: ", colnames(counts_df)[1])
}

if (ncol(gene_anno) < 2 || colnames(gene_anno)[1] != "gene_id") {
  stop("First column of gene annotation must be gene_id.")
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

required_anno_cols <- c("gene_id", "seqname", "start", "end", "strand", "midpoint")
missing_anno_cols <- setdiff(required_anno_cols, colnames(gene_anno))
if (length(missing_anno_cols) > 0) {
  stop("Gene annotation missing required columns: ", paste(missing_anno_cols, collapse = ","))
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

sample_id <- colnames(counts)
lib_size <- colSums(counts)

if (any(lib_size <= 0)) {
  bad <- sample_id[lib_size <= 0]
  stop("Samples with zero library size found: ", paste(bad, collapse = ","))
}

detected_genes_count_gt0 <- colSums(counts > 0)
detected_genes_count_ge10 <- colSums(counts >= 10)

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

if (file.exists(sample_meta) && file.info(sample_meta)$size > 0) {
  meta <- read.delim(sample_meta, check.names = FALSE)
  keep_cols <- intersect(
    c("sample_id", "animal_id", "term", "raw_data_name", "counts_file", "bam_file", "bai_file"),
    colnames(meta)
  )
  if ("sample_id" %in% keep_cols) {
    sample_qc <- merge(
      meta[, keep_cols, drop = FALSE],
      sample_qc,
      by = "sample_id",
      all.y = TRUE,
      sort = FALSE
    )
    sample_qc <- sample_qc[match(sample_id, sample_qc$sample_id), , drop = FALSE]
  }
}

gene_qc <- data.frame(
  gene_id = gene_id,
  gene_name = if ("gene_name" %in% colnames(gene_anno)) gene_anno$gene_name else gene_id,
  seqname = gene_anno$seqname,
  start = gene_anno$start,
  end = gene_anno$end,
  strand = gene_anno$strand,
  midpoint = gene_anno$midpoint,
  is_gene_LOC = grepl("^gene-LOC", gene_id),
  total_count = as.numeric(rowSums(counts)),
  mean_count = as.numeric(rowMeans(counts)),
  samples_count_gt0 = as.numeric(rowSums(counts > 0)),
  samples_count_ge10 = as.numeric(rowSums(counts >= 10)),
  samples_cpm_ge1 = as.numeric(rowSums(cpm >= 1)),
  keep_cpm1_min2 = keep_cpm1_min2,
  stringsAsFactors = FALSE
)

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
  paste0("input_matrix=", in_matrix),
  paste0("input_gene_annotation=", in_gene_annotation),
  paste0("sample_meta=", sample_meta),
  paste0("samples=", ncol(counts)),
  paste0("input_genes=", nrow(counts)),
  paste0("gene_LOC_rows_kept=", sum(grepl("^gene-LOC", gene_id))),
  paste0("MissingGeneID_rows=", sum(gene_id == "MissingGeneID")),
  paste0("positive_genes_for_size_factors=", nrow(positive_counts)),
  paste0("output_genes_all=", nrow(log2_norm)),
  paste0("output_genes_cpm1_min2=", nrow(log2_norm_filtered)),
  paste0("gene_LOC_rows_cpm1_min2=", sum(grepl("^gene-LOC", rownames(log2_norm_filtered)))),
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

writeLines(log_lines, out_check_log)
writeLines(log_lines)
