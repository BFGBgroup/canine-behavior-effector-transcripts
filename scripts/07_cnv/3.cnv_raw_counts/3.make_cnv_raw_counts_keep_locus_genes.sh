#!/usr/bin/env bash
set -euo pipefail


PROJECT_DIR="/path/to/project/legacy/five-dogs/5.cnv"
IN_DIR="${PROJECT_DIR}/2.expression_matrix"
OUT_DIR="${PROJECT_DIR}/3.cnv_raw_counts"

IN_MATRIX="${IN_DIR}/bulk_rna.star.reverse_stranded.raw_counts.keep_locus_genes.tsv"
IN_GENE_ANNOTATION="${IN_DIR}/CanFam3.1.patched_gene_annotation.for_expression_matrix.tsv"

OUT_MATRIX="${OUT_DIR}/bulk_rna.star.reverse_stranded.cnv_raw_counts.keep_locus_genes.tsv"
OUT_GENE_ANNOTATION="${OUT_DIR}/bulk_rna.star.reverse_stranded.cnv_raw_counts.keep_locus_genes.gene_annotation.tsv"
MISSING_ANNOTATION_GENES="${OUT_DIR}/bulk_rna.star.reverse_stranded.cnv_raw_counts.genes_missing_coordinate_annotation.tsv"
EXTRA_ANNOTATION_GENES="${OUT_DIR}/bulk_rna.star.reverse_stranded.cnv_raw_counts.annotation_genes_not_in_matrix.tsv"
CHECK_LOG="${OUT_DIR}/bulk_rna.star.reverse_stranded.cnv_raw_counts.keep_locus_genes.check.log"

mkdir -p "${OUT_DIR}"

if [[ ! -s "${IN_MATRIX}" ]]; then
    echo "[ERROR] Input expression matrix not found or empty: ${IN_MATRIX}" >&2
    echo "[ERROR] Please run scripts/2.expression_matrix/2.make_expression_matrix_keep_locus_genes.sh first." >&2
    exit 1
fi

if [[ ! -s "${IN_GENE_ANNOTATION}" ]]; then
    echo "[ERROR] Input gene annotation table not found or empty: ${IN_GENE_ANNOTATION}" >&2
    echo "[ERROR] Please run scripts/2.expression_matrix/2.make_expression_matrix_keep_locus_genes.sh first." >&2
    exit 1
fi

header_first_field="$(awk -F'\t' 'NR == 1 { print $1; exit }' "${IN_MATRIX}")"
if [[ "${header_first_field}" != "gene_id" ]]; then
    echo "[ERROR] First column of input matrix must be gene_id, got: ${header_first_field}" >&2
    exit 1
fi

annotation_header_first_field="$(awk -F'\t' 'NR == 1 { print $1; exit }' "${IN_GENE_ANNOTATION}")"
if [[ "${annotation_header_first_field}" != "gene_id" ]]; then
    echo "[ERROR] First column of gene annotation must be gene_id, got: ${annotation_header_first_field}" >&2
    exit 1
fi

input_genes=$(( $(wc -l < "${IN_MATRIX}") - 1 ))
sample_columns="$(awk -F'\t' 'NR == 1 { print NF - 1; exit }' "${IN_MATRIX}")"
annotation_genes=$(( $(wc -l < "${IN_GENE_ANNOTATION}") - 1 ))

if [[ "${input_genes}" -le 0 ]]; then
    echo "[ERROR] Input matrix has no gene rows: ${IN_MATRIX}" >&2
    exit 1
fi

if [[ "${sample_columns}" -le 0 ]]; then
    echo "[ERROR] Input matrix has no sample columns: ${IN_MATRIX}" >&2
    exit 1
fi

duplicate_matrix_genes="$(awk -F'\t' 'NR > 1 { seen[$1]++ } END { for (g in seen) if (seen[g] > 1) dup++; print dup + 0 }' "${IN_MATRIX}")"
duplicate_annotation_genes="$(awk -F'\t' 'NR > 1 { seen[$1]++ } END { for (g in seen) if (seen[g] > 1) dup++; print dup + 0 }' "${IN_GENE_ANNOTATION}")"
missing_geneid_rows="$(awk -F'\t' 'NR > 1 && $1 == "MissingGeneID" { n++ } END { print n + 0 }' "${IN_MATRIX}")"
gene_loc_rows="$(awk -F'\t' 'NR > 1 && $1 ~ /^gene-LOC/ { n++ } END { print n + 0 }' "${IN_MATRIX}")"

if [[ "${duplicate_matrix_genes}" -ne 0 ]]; then
    echo "[ERROR] Duplicate gene_id rows found in input matrix: ${duplicate_matrix_genes}" >&2
    exit 1
fi

if [[ "${duplicate_annotation_genes}" -ne 0 ]]; then
    echo "[ERROR] Duplicate gene_id rows found in gene annotation table: ${duplicate_annotation_genes}" >&2
    exit 1
fi

if [[ "${missing_geneid_rows}" -ne 0 ]]; then
    echo "[ERROR] Matrix contains MissingGeneID rows: ${missing_geneid_rows}" >&2
    echo "[ERROR] MissingGeneID has no reliable genomic coordinate, so it is not valid for CNV analysis." >&2
    echo "[ERROR] Please check whether STAR index and STAR GeneCounts used the patched GTF." >&2
    exit 1
fi


awk -F'\t' '
NR == FNR {
    if (NR > 1) {
        annotated[$1] = 1
    }
    next
}
FNR > 1 && !($1 in annotated) {
    print $1
}
' "${IN_GENE_ANNOTATION}" "${IN_MATRIX}" > "${MISSING_ANNOTATION_GENES}"

missing_annotation_genes="$(wc -l < "${MISSING_ANNOTATION_GENES}")"
if [[ "${missing_annotation_genes}" -ne 0 ]]; then
    echo "[ERROR] Matrix genes without coordinate annotation were found: ${missing_annotation_genes}" >&2
    echo "[ERROR] See: ${MISSING_ANNOTATION_GENES}" >&2
    exit 1
fi


awk -F'\t' '
NR == FNR {
    if (NR > 1) {
        in_matrix[$1] = 1
    }
    next
}
FNR > 1 && !($1 in in_matrix) {
    print $1
}
' "${IN_MATRIX}" "${IN_GENE_ANNOTATION}" > "${EXTRA_ANNOTATION_GENES}"

extra_annotation_genes="$(wc -l < "${EXTRA_ANNOTATION_GENES}")"


cp "${IN_MATRIX}" "${OUT_MATRIX}"


awk -F'\t' '
BEGIN { OFS = "\t" }
NR == FNR {
    if (NR == 1) {
        annotation_header = $0
        next
    }
    annotation[$1] = $0
    next
}
FNR == 1 {
    print annotation_header
    next
}
{
    print annotation[$1]
}
' "${IN_GENE_ANNOTATION}" "${OUT_MATRIX}" > "${OUT_GENE_ANNOTATION}"

output_genes=$(( $(wc -l < "${OUT_MATRIX}") - 1 ))
output_annotation_genes=$(( $(wc -l < "${OUT_GENE_ANNOTATION}") - 1 ))

if [[ "${output_genes}" -ne "${input_genes}" ]]; then
    echo "[ERROR] Output matrix gene count differs from input matrix." >&2
    echo "[ERROR] input_genes=${input_genes}, output_genes=${output_genes}" >&2
    exit 1
fi

if [[ "${output_annotation_genes}" -ne "${output_genes}" ]]; then
    echo "[ERROR] Output annotation gene count differs from output matrix." >&2
    echo "[ERROR] output_genes=${output_genes}, output_annotation_genes=${output_annotation_genes}" >&2
    exit 1
fi

if ! diff -q \
    <(awk -F'\t' 'NR > 1 { print $1 }' "${OUT_MATRIX}") \
    <(awk -F'\t' 'NR > 1 { print $1 }' "${OUT_GENE_ANNOTATION}") \
    >/dev/null; then
    echo "[ERROR] Output matrix and output gene annotation are not in the same gene order." >&2
    exit 1
fi

{
    echo "status=OK"
    echo "project_dir=${PROJECT_DIR}"
    echo "input_matrix=${IN_MATRIX}"
    echo "input_gene_annotation=${IN_GENE_ANNOTATION}"
    echo "output_matrix=${OUT_MATRIX}"
    echo "output_gene_annotation=${OUT_GENE_ANNOTATION}"
    echo "sample_columns=${sample_columns}"
    echo "input_genes=${input_genes}"
    echo "output_genes=${output_genes}"
    echo "gene_rows_removed=0"
    echo "gene_LOC_rows_kept=${gene_loc_rows}"
    echo "MissingGeneID_rows=${missing_geneid_rows}"
    echo "matrix_genes_missing_coordinate_annotation=${missing_annotation_genes}"
    echo "annotation_genes_not_in_matrix=${extra_annotation_genes}"
    echo "duplicate_matrix_gene_ids=${duplicate_matrix_genes}"
    echo "duplicate_annotation_gene_ids=${duplicate_annotation_genes}"
    echo "missing_annotation_genes_file=${MISSING_ANNOTATION_GENES}"
    echo "extra_annotation_genes_file=${EXTRA_ANNOTATION_GENES}"
} > "${CHECK_LOG}"

echo "[OK] CNV raw count matrix generated:"
echo "     ${OUT_MATRIX}"
echo "[OK] Ordered gene coordinate annotation generated:"
echo "     ${OUT_GENE_ANNOTATION}"
echo "[OK] Check log:"
echo "     ${CHECK_LOG}"
