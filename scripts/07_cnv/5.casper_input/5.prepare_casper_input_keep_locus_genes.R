#!/usr/bin/env Rscript


options(stringsAsFactors = FALSE)

project_dir <- "/path/to/project/legacy/five-dogs/5.cnv"


use_cpm1_min2_filtered_matrix <- FALSE

norm_dir <- file.path(project_dir, "4.normalized_expression")
out_dir <- file.path(project_dir, "5.casper_input")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (use_cpm1_min2_filtered_matrix) {
  expr_matrix <- file.path(
    norm_dir,
    "bulk_rna.cnv_expression.keep_locus_genes.log2_sizefactor_norm_counts.cpm1_min2.tsv"
  )
  gene_annotation <- file.path(
    norm_dir,
    "bulk_rna.cnv_expression.keep_locus_genes.gene_annotation.cpm1_min2.tsv"
  )
  output_tag <- "keep_locus_genes.cpm1_min2"
} else {
  expr_matrix <- file.path(
    norm_dir,
    "bulk_rna.cnv_expression.keep_locus_genes.log2_sizefactor_norm_counts.tsv"
  )
  gene_annotation <- file.path(
    norm_dir,
    "bulk_rna.cnv_expression.keep_locus_genes.gene_annotation.tsv"
  )
  output_tag <- "keep_locus_genes.all_genes"
}

sample_qc <- file.path(
  norm_dir,
  "bulk_rna.cnv_expression.keep_locus_genes.sample_qc.tsv"
)

out_expr <- file.path(out_dir, paste0("bulk_rna.cnv_expression.", output_tag, ".log2_norm.casper_ordered.tsv"))
out_annotation <- file.path(out_dir, paste0("CanFam3.1.", output_tag, ".gene_annotation.for_casper.tsv"))
out_expr_primary <- file.path(out_dir, paste0("bulk_rna.cnv_expression.", output_tag, ".log2_norm.primary_chr.casper_ordered.tsv"))
out_annotation_primary <- file.path(out_dir, paste0("CanFam3.1.", output_tag, ".gene_annotation.primary_chr.for_casper.tsv"))
out_missing <- file.path(out_dir, paste0("bulk_rna.cnv_expression.", output_tag, ".genes_missing_annotation.tsv"))
out_sample_meta <- file.path(out_dir, paste0("bulk_rna.samples.", output_tag, ".for_casper.tsv"))
out_check <- file.path(out_dir, paste0("prepare_casper_input.", output_tag, ".check.log"))

chr_order_key <- function(chr) {
  x <- sub("^chr", "", chr)
  suppressWarnings(num <- as.integer(x))
  ifelse(
    !is.na(num),
    num,
    ifelse(x == "X", 100L, ifelse(x == "Y", 101L, ifelse(x %in% c("M", "MT"), 102L, 999L)))
  )
}

is_primary_dog_chr <- function(chr) {
  x <- sub("^chr", "", chr)
  suppressWarnings(num <- as.integer(x))
  (!is.na(num) & num >= 1L & num <= 38L) | x %in% c("X", "Y", "M", "MT")
}

assert_nonempty_file <- function(path, label) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop(label, " not found or empty: ", path)
  }
}

assert_nonempty_file(expr_matrix, "Expression matrix")
assert_nonempty_file(gene_annotation, "Gene annotation")
assert_nonempty_file(sample_qc, "Sample QC")

expr <- read.delim(expr_matrix, check.names = FALSE)
anno <- read.delim(gene_annotation, check.names = FALSE)
samples <- read.delim(sample_qc, check.names = FALSE)

if (ncol(expr) < 2 || colnames(expr)[1] != "gene_id") {
  stop("Expression matrix must start with gene_id and contain sample columns.")
}

if (ncol(anno) < 2 || colnames(anno)[1] != "gene_id") {
  stop("Gene annotation must start with gene_id.")
}

if (!"sample_id" %in% colnames(samples)) {
  stop("Sample QC must contain sample_id column.")
}

if (anyDuplicated(expr$gene_id) > 0) {
  stop("Duplicated gene_id in expression matrix: ", expr$gene_id[duplicated(expr$gene_id)][1])
}

if (anyDuplicated(anno$gene_id) > 0) {
  stop("Duplicated gene_id in gene annotation: ", anno$gene_id[duplicated(anno$gene_id)][1])
}

if (!identical(expr$gene_id, anno$gene_id)) {
  missing_from_anno <- setdiff(expr$gene_id, anno$gene_id)
  write.table(
    data.frame(gene_id = missing_from_anno),
    out_missing,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  stop("Expression matrix and gene annotation are not in the same gene order.")
}

if (any(expr$gene_id == "MissingGeneID")) {
  stop("MissingGeneID is present in expression matrix. It has no reliable genomic coordinate.")
}

required_anno_cols <- c("gene_id", "gene_name", "seqname", "start", "end", "strand", "midpoint")
missing_anno_cols <- setdiff(required_anno_cols, colnames(anno))
if (length(missing_anno_cols) > 0) {
  stop("Gene annotation missing required columns: ", paste(missing_anno_cols, collapse = ","))
}

expr_values <- as.matrix(expr[, -1, drop = FALSE])
suppressWarnings(storage.mode(expr_values) <- "numeric")
if (any(!is.finite(expr_values))) {
  stop("Expression matrix contains NA/NaN/Inf or non-numeric values.")
}

sample_id <- colnames(expr)[-1]
if (!identical(sample_id, samples$sample_id)) {
  samples <- samples[match(sample_id, samples$sample_id), , drop = FALSE]
  if (!identical(sample_id, samples$sample_id)) {
    stop("Sample QC sample_id values do not match expression matrix columns.")
  }
}

anno$start <- as.integer(anno$start)
anno$end <- as.integer(anno$end)
anno$midpoint <- as.integer(anno$midpoint)

if (any(is.na(anno$seqname) | anno$seqname == "" | is.na(anno$start) | is.na(anno$end) | is.na(anno$midpoint))) {
  stop("Gene annotation contains missing chromosome/start/end/midpoint values.")
}

casper_annotation <- data.frame(
  Gene = anno$gene_id,
  GeneSymbol = ifelse(is.na(anno$gene_name) | anno$gene_name == "", anno$gene_id, anno$gene_name),
  Chr = anno$seqname,
  start = anno$start,
  end = anno$end,
  strand = anno$strand,
  Position = anno$midpoint,
  cytoband = anno$seqname,
  isCentromer = "no",
  source_feature = if ("source_feature" %in% colnames(anno)) anno$source_feature else NA_character_,
  is_gene_LOC = grepl("^gene-LOC", anno$gene_id),
  is_primary_dog_chr = is_primary_dog_chr(anno$seqname),
  stringsAsFactors = FALSE
)

missing_genes <- setdiff(expr$gene_id, casper_annotation$Gene)
write.table(
  data.frame(gene_id = missing_genes),
  out_missing,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ord <- order(
  chr_order_key(casper_annotation$Chr),
  casper_annotation$Chr,
  casper_annotation$start,
  casper_annotation$end,
  casper_annotation$Gene
)

expr_ordered <- expr[ord, , drop = FALSE]
annotation_ordered <- casper_annotation[ord, , drop = FALSE]

primary_keep <- annotation_ordered$is_primary_dog_chr
expr_primary <- expr_ordered[primary_keep, , drop = FALSE]
annotation_primary <- annotation_ordered[primary_keep, , drop = FALSE]

if (!identical(expr_ordered$gene_id, annotation_ordered$Gene)) {
  stop("Ordered expression and annotation gene order do not match.")
}

if (!identical(expr_primary$gene_id, annotation_primary$Gene)) {
  stop("Primary-chromosome expression and annotation gene order do not match.")
}

write.table(expr_ordered, out_expr, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(annotation_ordered, out_annotation, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(expr_primary, out_expr_primary, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(annotation_primary, out_annotation_primary, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(samples, out_sample_meta, sep = "\t", quote = FALSE, row.names = FALSE)

log_lines <- c(
  paste0("status=OK"),
  paste0("use_cpm1_min2_filtered_matrix=", use_cpm1_min2_filtered_matrix),
  paste0("expression_matrix=", expr_matrix),
  paste0("gene_annotation=", gene_annotation),
  paste0("sample_qc=", sample_qc),
  paste0("input_expression_genes=", nrow(expr)),
  paste0("input_samples=", length(sample_id)),
  paste0("gene_LOC_rows=", sum(casper_annotation$is_gene_LOC)),
  paste0("MissingGeneID_rows=", sum(expr$gene_id == "MissingGeneID")),
  paste0("missing_annotation_genes=", length(missing_genes)),
  paste0("primary_chromosome_genes=", sum(annotation_ordered$is_primary_dog_chr)),
  paste0("non_primary_scaffold_genes=", sum(!annotation_ordered$is_primary_dog_chr)),
  paste0("primary_chromosome_gene_LOC_rows=", sum(annotation_primary$is_gene_LOC)),
  paste0("all_casper_expression=", out_expr),
  paste0("all_casper_annotation=", out_annotation),
  paste0("primary_chr_casper_expression=", out_expr_primary),
  paste0("primary_chr_casper_annotation=", out_annotation_primary),
  paste0("missing_genes=", out_missing),
  paste0("sample_meta=", out_sample_meta)
)

writeLines(log_lines, out_check)
writeLines(log_lines)
