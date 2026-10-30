#!/usr/bin/bash
set -euo pipefail


PROJECT_DIR="/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation"
INPUT_RUN_ID="20260821_134953"
IN_DIR="${PROJECT_DIR}/1.2.2.expression_matrix/${INPUT_RUN_ID}"
INPUT_COMPLETION_MARKER="${IN_DIR}/execution_completion.marker"
OUTPUT_ROOT="${PROJECT_DIR}/1.2.3.cnv_raw_counts"

SCRIPT_DIR="/path/to/repository/scripts/08_expression_ase_integration/preparation/1.2.3.cnv_raw_counts"
LOG_DIR="${SCRIPT_DIR}/log"
RUN_ID="${RUN_ID:-$(/usr/bin/date '+%Y%m%d_%H%M%S')}"
if [[ ! "${RUN_ID}" =~ ^[0-9]{8}_[0-9]{6}$ ]]; then
    echo "[ERROR] RUN_ID must use YYYYMMDD_HHMMSS format: ${RUN_ID}" >&2
    exit 1
fi
OUT_DIR="${OUTPUT_ROOT}/${RUN_ID}"
MASTER_LOG="${LOG_DIR}/1.2.3.cnv_raw_counts.${RUN_ID}.log"
COMPLETION_MARKER="${OUT_DIR}/execution_completion.marker"


EXPECTED_SAMPLE_COLUMNS=19
EXPECTED_INPUT_GENES=40825
EXPECTED_ANNOTATION_GENES=40924
EXPECTED_GENE_LOC_ROWS=18858
EXPECTED_DOG5_NEW_GENE_ROWS=3828
EXPECTED_EXTRA_ANNOTATION_GENES=99
EXPECTED_ANNOTATION_HEADER=$'gene_id\tgene_name\tseqname\tstart\tend\tstrand\tmidpoint\tsource_feature'

AWK_BIN="/usr/bin/awk"
CP_BIN="/usr/bin/cp"
DATE_BIN="/usr/bin/date"
DIFF_BIN="/usr/bin/diff"
HEAD_BIN="/usr/bin/head"
MKDIR_BIN="/usr/bin/mkdir"
TR_BIN="/usr/bin/tr"
WC_BIN="/usr/bin/wc"

IN_MATRIX="${IN_DIR}/bulk_rna.star.reverse_stranded.raw_counts.manual_new_gene.tsv"
IN_GENE_ANNOTATION="${IN_DIR}/CanFam3.1.manual_new_gene_annotation.for_expression_matrix.tsv"

OUT_MATRIX="${OUT_DIR}/bulk_rna.star.reverse_stranded.cnv_raw_counts.manual_new_gene.tsv"
OUT_GENE_ANNOTATION="${OUT_DIR}/bulk_rna.star.reverse_stranded.cnv_raw_counts.manual_new_gene.gene_annotation.tsv"
MISSING_ANNOTATION_GENES="${OUT_DIR}/bulk_rna.star.reverse_stranded.cnv_raw_counts.manual_new_gene.genes_missing_coordinate_annotation.tsv"
EXTRA_ANNOTATION_GENES="${OUT_DIR}/bulk_rna.star.reverse_stranded.cnv_raw_counts.manual_new_gene.annotation_genes_not_in_matrix.tsv"
CHECK_LOG="${OUT_DIR}/bulk_rna.star.reverse_stranded.cnv_raw_counts.manual_new_gene.check.log"

log() {
    echo "[$(${DATE_BIN} '+%F %T') ] $*"
}

${MKDIR_BIN} -p "${OUTPUT_ROOT}" "${LOG_DIR}"
if [[ -e "${OUT_DIR}" || -e "${MASTER_LOG}" ]]; then
    echo "[ERROR] this run run directory or log already exists, refuse to overwrite: ${RUN_ID}" >&2
    exit 1
fi
${MKDIR_BIN} "${OUT_DIR}"
: > "${MASTER_LOG}"
exec >> "${MASTER_LOG}" 2>&1

on_error() {
    local rc=$?
    trap - ERR
    log "[ERROR] Script failed at line ${BASH_LINENO[0]} with exit code ${rc}"
    exit "${rc}"
}
trap on_error ERR

on_exit() {
    local rc=$?
    trap - EXIT
    if [[ "${rc}" -ne 0 ]]; then
        log "execution_status=failed"
    fi
    exit "${rc}"
}
trap on_exit EXIT

log "[start] prepare finemaping CNV raw-count input"
log "[config] RUN_ID=${RUN_ID}"
log "[config] PROJECT_DIR=${PROJECT_DIR}"
log "[config] INPUT_RUN_ID=${INPUT_RUN_ID}"
log "[config] IN_DIR=${IN_DIR}"
log "[config] OUT_DIR=${OUT_DIR}"
log "[config] MASTER_LOG=${MASTER_LOG}"
log "[config] IN_MATRIX=${IN_MATRIX}"
log "[config] IN_GENE_ANNOTATION=${IN_GENE_ANNOTATION}"

log "[stage:validate-inputs] check input matrix and gene annotation table"
if [[ ! -s "${IN_MATRIX}" ]]; then
    echo "[ERROR] Input expression matrix not found or empty: ${IN_MATRIX}" >&2
    echo "[ERROR] Please generate 1.2.2.expression_matrix results first." >&2
    exit 1
fi

if [[ ! -s "${IN_GENE_ANNOTATION}" ]]; then
    echo "[ERROR] Input gene annotation table not found or empty: ${IN_GENE_ANNOTATION}" >&2
    echo "[ERROR] Please generate 1.2.2.expression_matrix results first." >&2
    exit 1
fi

if [[ ! -s "${INPUT_COMPLETION_MARKER}" ]] || ! ${AWK_BIN} -F'\t' '$1 == "execution_status" && $2 == "success" {ok=1} END {exit(ok ? 0 : 1)}' "${INPUT_COMPLETION_MARKER}"; then
    echo "[ERROR] upstream 1.2.2 run no valid of success completion marker: ${INPUT_COMPLETION_MARKER}" >&2
    exit 1
fi

header_first_field="$(${AWK_BIN} -F'\t' 'NR == 1 { print $1; exit }' "${IN_MATRIX}")"
if [[ "${header_first_field}" != "gene_id" ]]; then
    echo "[ERROR] First column of input matrix must be gene_id, got: ${header_first_field}" >&2
    exit 1
fi

annotation_header="$(${HEAD_BIN} -n 1 "${IN_GENE_ANNOTATION}" | ${TR_BIN} -d '\r')"
if [[ "${annotation_header}" != "${EXPECTED_ANNOTATION_HEADER}" ]]; then
    echo "[ERROR] Unexpected gene annotation header: ${annotation_header}" >&2
    echo "[ERROR] Expected header: ${EXPECTED_ANNOTATION_HEADER}" >&2
    exit 1
fi

input_genes=$(( $(${WC_BIN} -l < "${IN_MATRIX}") - 1 ))
sample_columns="$(${AWK_BIN} -F'\t' 'NR == 1 { print NF - 1; exit }' "${IN_MATRIX}")"
annotation_genes=$(( $(${WC_BIN} -l < "${IN_GENE_ANNOTATION}") - 1 ))
log "[stage:validate-inputs] input_genes=${input_genes}; sample_columns=${sample_columns}; annotation_genes=${annotation_genes}"

if [[ "${input_genes}" -ne "${EXPECTED_INPUT_GENES}" ]]; then
    echo "[ERROR] Unexpected input gene count: ${input_genes}; expected ${EXPECTED_INPUT_GENES}" >&2
    exit 1
fi

if [[ "${sample_columns}" -ne "${EXPECTED_SAMPLE_COLUMNS}" ]]; then
    echo "[ERROR] Unexpected sample column count: ${sample_columns}; expected ${EXPECTED_SAMPLE_COLUMNS}" >&2
    exit 1
fi

if [[ "${annotation_genes}" -ne "${EXPECTED_ANNOTATION_GENES}" ]]; then
    echo "[ERROR] Unexpected annotation gene count: ${annotation_genes}; expected ${EXPECTED_ANNOTATION_GENES}" >&2
    exit 1
fi

duplicate_matrix_genes="$(${AWK_BIN} -F'\t' 'NR > 1 { seen[$1]++ } END { for (g in seen) if (seen[g] > 1) dup++; print dup + 0 }' "${IN_MATRIX}")"
duplicate_annotation_genes="$(${AWK_BIN} -F'\t' 'NR > 1 { seen[$1]++ } END { for (g in seen) if (seen[g] > 1) dup++; print dup + 0 }' "${IN_GENE_ANNOTATION}")"
missing_geneid_rows="$(${AWK_BIN} -F'\t' 'NR > 1 && $1 == "MissingGeneID" { n++ } END { print n + 0 }' "${IN_MATRIX}")"
gene_loc_rows="$(${AWK_BIN} -F'\t' 'NR > 1 && $1 ~ /^gene-LOC/ { n++ } END { print n + 0 }' "${IN_MATRIX}")"
dog5_new_gene_rows="$(${AWK_BIN} -F'\t' 'NR > 1 && $1 ~ /^DOG5_NEWGENE_/ { n++ } END { print n + 0 }' "${IN_MATRIX}")"
unified_new_gene_rows="$(${AWK_BIN} -F'\t' 'NR > 1 && $1 ~ /^UNIFIED_NEWGENE_/ { n++ } END { print n + 0 }' "${IN_MATRIX}")"

if [[ "${gene_loc_rows}" -ne "${EXPECTED_GENE_LOC_ROWS}" ]]; then
    echo "[ERROR] Unexpected gene-LOC row count: ${gene_loc_rows}; expected ${EXPECTED_GENE_LOC_ROWS}" >&2
    exit 1
fi

if [[ "${dog5_new_gene_rows}" -ne "${EXPECTED_DOG5_NEW_GENE_ROWS}" ]]; then
    echo "[ERROR] Unexpected DOG5_NEWGENE row count: ${dog5_new_gene_rows}; expected ${EXPECTED_DOG5_NEW_GENE_ROWS}" >&2
    exit 1
fi

if [[ "${unified_new_gene_rows}" -ne 0 ]]; then
    echo "[ERROR] five-only the matrix must not contain UNIFIED_NEWGENE row: ${unified_new_gene_rows}" >&2
    exit 1
fi

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
    echo "[ERROR] Please check whether STAR index and STAR GeneCounts used the Manual-annotation GTF." >&2
    exit 1
fi
log "[stage:validate-inputs] duplicate_matrix_gene_ids=${duplicate_matrix_genes}; duplicate_annotation_gene_ids=${duplicate_annotation_genes}; MissingGeneID_rows=${missing_geneid_rows}"
log "[stage:validate-inputs] gene_LOC_rows_kept=${gene_loc_rows}; dog5_new_gene_rows_kept=${dog5_new_gene_rows}; legacy_unified_new_gene_rows=${unified_new_gene_rows}"


log "[stage:annotation-check] check matrix gene whether all has coordinates annotation"
${AWK_BIN} -F'\t' '
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

missing_annotation_genes="$(${WC_BIN} -l < "${MISSING_ANNOTATION_GENES}")"
if [[ "${missing_annotation_genes}" -ne 0 ]]; then
    echo "[ERROR] Matrix genes without coordinate annotation were found: ${missing_annotation_genes}" >&2
    echo "[ERROR] See: ${MISSING_ANNOTATION_GENES}" >&2
    exit 1
fi
log "[stage:annotation-check] matrix_genes_missing_coordinate_annotation=${missing_annotation_genes}"


log "[stage:audit] record annotation not included in the matrix from gene"
${AWK_BIN} -F'\t' '
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

extra_annotation_genes="$(${WC_BIN} -l < "${EXTRA_ANNOTATION_GENES}")"
if [[ "${extra_annotation_genes}" -ne "${EXPECTED_EXTRA_ANNOTATION_GENES}" ]]; then
    echo "[ERROR] Unexpected annotation-only gene count: ${extra_annotation_genes}; expected ${EXPECTED_EXTRA_ANNOTATION_GENES}" >&2
    exit 1
fi
log "[stage:audit] annotation_genes_not_in_matrix=${extra_annotation_genes}"


log "[stage:write] write CNV raw count matrix"
${CP_BIN} "${IN_MATRIX}" "${OUT_MATRIX}"


log "[stage:write] by matrix order reorder and write gene coordinates annotation table"
${AWK_BIN} -F'\t' '
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

output_genes=$(( $(${WC_BIN} -l < "${OUT_MATRIX}") - 1 ))
output_annotation_genes=$(( $(${WC_BIN} -l < "${OUT_GENE_ANNOTATION}") - 1 ))

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

if ! ${DIFF_BIN} -q \
    <(${AWK_BIN} -F'\t' 'NR > 1 { print $1 }' "${OUT_MATRIX}") \
    <(${AWK_BIN} -F'\t' 'NR > 1 { print $1 }' "${OUT_GENE_ANNOTATION}") \
    >/dev/null; then
    echo "[ERROR] Output matrix and output gene annotation are not in the same gene order." >&2
    exit 1
fi
log "[stage:final-check] output_genes=${output_genes}; output_annotation_genes=${output_annotation_genes}; gene_order=matched"

{
    echo "status=OK"
    echo "run_id=${RUN_ID}"
    echo "project_dir=${PROJECT_DIR}"
    echo "input_run_id=${INPUT_RUN_ID}"
    echo "input_completion_marker=${INPUT_COMPLETION_MARKER}"
    echo "output_dir=${OUT_DIR}"
    echo "input_matrix=${IN_MATRIX}"
    echo "input_gene_annotation=${IN_GENE_ANNOTATION}"
    echo "output_matrix=${OUT_MATRIX}"
    echo "output_gene_annotation=${OUT_GENE_ANNOTATION}"
    echo "sample_columns=${sample_columns}"
    echo "input_genes=${input_genes}"
    echo "output_genes=${output_genes}"
    echo "gene_rows_removed=0"
    echo "gene_LOC_rows_kept=${gene_loc_rows}"
    echo "dog5_new_gene_rows_kept=${dog5_new_gene_rows}"
    echo "legacy_unified_new_gene_rows=${unified_new_gene_rows}"
    echo "MissingGeneID_rows=${missing_geneid_rows}"
    echo "matrix_genes_missing_coordinate_annotation=${missing_annotation_genes}"
    echo "annotation_genes_not_in_matrix=${extra_annotation_genes}"
    echo "duplicate_matrix_gene_ids=${duplicate_matrix_genes}"
    echo "duplicate_annotation_gene_ids=${duplicate_annotation_genes}"
    echo "missing_annotation_genes_file=${MISSING_ANNOTATION_GENES}"
    echo "extra_annotation_genes_file=${EXTRA_ANNOTATION_GENES}"
} > "${CHECK_LOG}"

log "[summary] output_genes=${output_genes}; gene_rows_removed=0; gene_LOC_rows_kept=${gene_loc_rows}; dog5_new_gene_rows_kept=${dog5_new_gene_rows}; legacy_unified_new_gene_rows=${unified_new_gene_rows}"
log "[summary] matrix_genes_missing_coordinate_annotation=${missing_annotation_genes}; annotation_genes_not_in_matrix=${extra_annotation_genes}"
printf 'run_id\t%s\nexecution_status\tsuccess\ncompleted_at\t%s\n' \
    "${RUN_ID}" "$(${DATE_BIN} '+%F %T')" > "${COMPLETION_MARKER}"
log "[done] COMPLETION_MARKER=${COMPLETION_MARKER}"
log "[done] finemaping CNV raw-count input prepare completed; task completed / done"
log "execution_status=success"
log "reminder: execution completion does not establish output completeness, format, count or biological quality already passed; please in run after separately check result."

echo "[OK] CNV raw count matrix generated:"
echo "     ${OUT_MATRIX}"
echo "[OK] Ordered gene coordinate annotation generated:"
echo "     ${OUT_GENE_ANNOTATION}"
echo "[OK] Check log:"
echo "     ${CHECK_LOG}"
