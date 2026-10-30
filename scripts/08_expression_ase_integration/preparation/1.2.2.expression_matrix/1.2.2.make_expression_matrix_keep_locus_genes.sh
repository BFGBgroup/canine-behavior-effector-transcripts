#!/usr/bin/bash
set -euo pipefail


PROJECT_DIR="/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation"
STAR_RUN_ID="20260820_225724"
STAR_DIR="${PROJECT_DIR}/1.2.1star/${STAR_RUN_ID}"
STAR_COMPLETION_MARKER="${STAR_DIR}/execution_completion.marker"
OUTPUT_ROOT="${PROJECT_DIR}/1.2.2.expression_matrix"

SCRIPT_DIR="/path/to/repository/scripts/08_expression_ase_integration/preparation/1.2.2.expression_matrix"
LOG_DIR="${SCRIPT_DIR}/log"
RUN_ID="${RUN_ID:-$(/usr/bin/date '+%Y%m%d_%H%M%S')}"
if [[ ! "${RUN_ID}" =~ ^[0-9]{8}_[0-9]{6}$ ]]; then
    echo "[ERROR] RUN_ID must use YYYYMMDD_HHMMSS format: ${RUN_ID}" >&2
    exit 1
fi
OUT_DIR="${OUTPUT_ROOT}/${RUN_ID}"
MASTER_LOG="${LOG_DIR}/1.2.2.expression_matrix.${RUN_ID}.log"
COMPLETION_MARKER="${OUT_DIR}/execution_completion.marker"

SAMPLE_TSV="/path/to/project/legacy/five-dogs/data/statistics/five_dogs_fastq_pair_stats.tsv"
MANUAL_GTF="/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/GCF_000002285.3_CanFam3.1_genomic_chr.plus_five_only_new_gene.gtf"


EXPECTED_SAMPLES=19
EXPECTED_MATRIX_GENES=40825
EXPECTED_ANNOTATION_GENES=40924
EXPECTED_MATRIX_GENE_LOC_ROWS=18858
EXPECTED_ANNOTATION_GENE_LOC_ROWS=18940
EXPECTED_DOG5_NEW_GENE_ROWS=3828
EXPECTED_ANNOTATION_ONLY_GENES=99
EXPECTED_SAMPLE_HEADER=$'animal_id\tterm\traw_data_name'

AWK_BIN="/usr/bin/awk"
CAT_BIN="/usr/bin/cat"
CMP_BIN="/usr/bin/cmp"
CP_BIN="/usr/bin/cp"
CUT_BIN="/usr/bin/cut"
DATE_BIN="/usr/bin/date"
HEAD_BIN="/usr/bin/head"
MKDIR_BIN="/usr/bin/mkdir"
MKTEMP_BIN="/usr/bin/mktemp"
MV_BIN="/usr/bin/mv"
PASTE_BIN="/usr/bin/paste"
RM_BIN="/usr/bin/rm"
TAIL_BIN="/usr/bin/tail"
TR_BIN="/usr/bin/tr"
WC_BIN="/usr/bin/wc"

SAMPLES_FILE="${OUT_DIR}/bulk_rna.star.samples.tsv"
MISSING_OUTPUTS_FILE="${OUT_DIR}/bulk_rna.star.missing_outputs.tsv"
STRANDNESS_SUMMARY="${OUT_DIR}/star_gene_counts_strandness_summary.tsv"
STRANDNESS_OVERALL="${OUT_DIR}/star_gene_counts_strandness_overall.txt"
GENE_ANNOTATION="${OUT_DIR}/CanFam3.1.manual_new_gene_annotation.for_expression_matrix.tsv"
MATRIX_FILE="${OUT_DIR}/bulk_rna.star.reverse_stranded.raw_counts.manual_new_gene.tsv"
MISSING_ANNOTATION_GENES="${OUT_DIR}/bulk_rna.star.genes_missing_from_manual_new_gtf.tsv"
CHECK_LOG="${OUT_DIR}/bulk_rna.star.reverse_stranded.merge.manual_new_gene.check.log"


COUNT_COLUMN=4
COUNT_COLUMN_NAME="reverse_stranded"

log() {
    echo "[$(${DATE_BIN} '+%F %T ')] $*"
}

${MKDIR_BIN} -p "${OUTPUT_ROOT}" "${LOG_DIR}"
if [[ -e "${OUT_DIR}" || -e "${MASTER_LOG}" ]]; then
    echo "[ERROR] this run run directory or log already exists, refuse to overwrite: ${RUN_ID}" >&2
    exit 1
fi
${MKDIR_BIN} "${OUT_DIR}"
: > "${MASTER_LOG}"
exec >> "${MASTER_LOG}" 2>&1

gene_list_tmp=""
matrix_body_tmp=""
header_tmp=""

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

    [[ -z "${gene_list_tmp}" ]] || ${RM_BIN} -f "${gene_list_tmp}"
    [[ -z "${matrix_body_tmp}" ]] || ${RM_BIN} -f "${matrix_body_tmp}" "${matrix_body_tmp}.new"
    [[ -z "${header_tmp}" ]] || ${RM_BIN} -f "${header_tmp}"

    if [[ "${rc}" -ne 0 ]]; then
        log "execution_status=failed"
    fi
    exit "${rc}"
}
trap on_exit EXIT

log "[start] build finemaping expression raw-count matrix"
log "[config] RUN_ID=${RUN_ID}"
log "[config] PROJECT_DIR=${PROJECT_DIR}"
log "[config] STAR_RUN_ID=${STAR_RUN_ID}"
log "[config] STAR_DIR=${STAR_DIR}"
log "[config] OUT_DIR=${OUT_DIR}"
log "[config] MASTER_LOG=${MASTER_LOG}"
log "[config] SAMPLE_TSV=${SAMPLE_TSV}"
log "[config] MANUAL_GTF=${MANUAL_GTF}"
log "[config] COUNT_COLUMN=${COUNT_COLUMN} (${COUNT_COLUMN_NAME})"

if [[ ! -s "${SAMPLE_TSV}" ]]; then
    echo "[ERROR] sample TSV not found or empty: ${SAMPLE_TSV}" >&2
    exit 1
fi

sample_header="$(${HEAD_BIN} -n 1 "${SAMPLE_TSV}" | ${TR_BIN} -d '\r')"
if [[ "${sample_header}" != "${EXPECTED_SAMPLE_HEADER}" ]]; then
    echo "[ERROR] Unexpected sample TSV header: ${sample_header}" >&2
    echo "[ERROR] Expected header: ${EXPECTED_SAMPLE_HEADER}" >&2
    exit 1
fi

if [[ ! -d "${STAR_DIR}" ]]; then
    echo "[ERROR] STAR output directory not found: ${STAR_DIR}" >&2
    exit 1
fi

if [[ ! -s "${STAR_COMPLETION_MARKER}" ]] || ! ${AWK_BIN} -F'\t' '$1 == "execution_status" && $2 == "success" {ok=1} END {exit(ok ? 0 : 1)}' "${STAR_COMPLETION_MARKER}"; then
    echo "[ERROR] upstream STAR run no valid of success completion marker: ${STAR_COMPLETION_MARKER}" >&2
    exit 1
fi

if [[ ! -s "${MANUAL_GTF}" ]]; then
    echo "[ERROR] Manual-annotation GTF not found or empty: ${MANUAL_GTF}" >&2
    echo "[ERROR] This script intentionally does not fall back to the original GTF." >&2
    exit 1
fi


log "[stage:sample-list] read sample table and construct sample list"
${AWK_BIN} -v star_dir="${STAR_DIR}" -F'\t' '
BEGIN { OFS="\t"; print "sample_id", "animal_id", "term", "raw_data_name", "counts_file", "bam_file", "bai_file" }
NR == 1 { next }
{
    animal_id = $1
    term = $2
    raw_data_name = $3
    gsub(/\r/, "", animal_id)
    gsub(/\r/, "", term)
    gsub(/\r/, "", raw_data_name)

    if (animal_id != "" && term != "") {
        sample_id = animal_id "_" term
        sample_dir = star_dir "/" animal_id "/" term "/align_star"
        counts_file = sample_dir "/" sample_id ".ReadsPerGene.out.tab"
        bam_file = sample_dir "/" sample_id ".Aligned.sortedByCoord.out.bam"
        bai_file = bam_file ".bai"
        print sample_id, animal_id, term, raw_data_name, counts_file, bam_file, bai_file
    }
}
' "${SAMPLE_TSV}" > "${SAMPLES_FILE}"

expected_samples=$(( $(${WC_BIN} -l < "${SAMPLES_FILE}") - 1 ))
if [[ "${expected_samples}" -ne "${EXPECTED_SAMPLES}" ]]; then
    echo "[ERROR] Unexpected sample count: ${expected_samples}; expected ${EXPECTED_SAMPLES}" >&2
    exit 1
fi

duplicate_sample_ids="$(${AWK_BIN} -F'\t' 'NR > 1 { seen[$1]++ } END { for (s in seen) if (seen[s] > 1) dup++; print dup + 0 }' "${SAMPLES_FILE}")"
if [[ "${duplicate_sample_ids}" -ne 0 ]]; then
    echo "[ERROR] Duplicate sample_id rows found: ${duplicate_sample_ids}" >&2
    exit 1
fi
log "[stage:sample-list] expected_samples=${expected_samples}; duplicate_sample_ids=${duplicate_sample_ids}; samples_file=${SAMPLES_FILE}"


log "[stage:check-inputs] check STAR output whether complete"
{
    echo -e "sample_id\tanimal_id\tterm\tmissing_item\tpath"
    ${TAIL_BIN} -n +2 "${SAMPLES_FILE}" | while IFS=$'\t' read -r sample_id animal_id term _raw_data_name counts_file bam_file bai_file; do

        if [[ ! -s "${counts_file}" ]]; then
            echo -e "${sample_id}\t${animal_id}\t${term}\tReadsPerGene.out.tab\t${counts_file}"
        else
            bad_column_rows="$(${AWK_BIN} -F'\t' 'NF != 4 { n++ } END { print n + 0 }' "${counts_file}")"
            [[ "${bad_column_rows}" -eq 0 ]] || echo -e "${sample_id}\t${animal_id}\t${term}\tinvalid_ReadsPerGene_column_count\t${counts_file}"
        fi
        [[ -s "${bam_file}" ]] || echo -e "${sample_id}\t${animal_id}\t${term}\tsorted_bam\t${bam_file}"
        [[ -s "${bai_file}" ]] || echo -e "${sample_id}\t${animal_id}\t${term}\tbam_index\t${bai_file}"
    done
} > "${MISSING_OUTPUTS_FILE}"

missing_outputs=$(( $(${WC_BIN} -l < "${MISSING_OUTPUTS_FILE}") - 1 ))
if [[ "${missing_outputs}" -ne 0 ]]; then
    echo "[ERROR] Missing STAR outputs were found. See: ${MISSING_OUTPUTS_FILE}" >&2
    exit 1
fi
log "[stage:check-inputs] missing_star_outputs=${missing_outputs}"


log "[stage:strandness] summary STAR GeneCounts strand specificity signal"
{
    echo -e "sample_id\tgene_rows\tmissing_geneid_rows\tunstranded_sum\tforward_sum\treverse_sum\treverse_over_forward"
    ${TAIL_BIN} -n +2 "${SAMPLES_FILE}" | while IFS=$'\t' read -r sample_id _animal_id _term _raw_data_name counts_file _bam_file _bai_file; do
        ${AWK_BIN} -v sample_id="${sample_id}" -F'\t' '
        BEGIN { OFS="\t" }
        $1 !~ /^N_/ {
            gene_rows++
            if ($1 == "MissingGeneID") {
                missing_geneid_rows++
            }
            unstranded += $2
            forward += $3
            reverse += $4
        }
        END {
            ratio = (forward > 0 ? reverse / forward : "NA")
            print sample_id, gene_rows + 0, missing_geneid_rows + 0, unstranded + 0, forward + 0, reverse + 0, ratio
        }
        ' "${counts_file}"
    done
} > "${STRANDNESS_SUMMARY}"

${AWK_BIN} -F'\t' '
BEGIN { OFS="\t" }
NR == 1 { next }
{
    samples++
    gene_rows[$2]++
    missing_geneid += $3
    unstranded += $4
    forward += $5
    reverse += $6
}
END {
    print "samples", samples + 0
    print "selected_count_column", "4"
    print "selected_count_type", "reverse_stranded"
    print "total_unstranded_counts", unstranded + 0
    print "total_forward_counts", forward + 0
    print "total_reverse_counts", reverse + 0
    print "reverse_over_forward", (forward > 0 ? reverse / forward : "NA")
    print "total_MissingGeneID_rows", missing_geneid + 0
    print "gene_row_count_classes", length(gene_rows)
    for (n in gene_rows) {
        print "gene_rows=" n, "sample_count=" gene_rows[n]
    }
}
' "${STRANDNESS_SUMMARY}" > "${STRANDNESS_OVERALL}"

strandness_samples=$(( $(${WC_BIN} -l < "${STRANDNESS_SUMMARY}") - 1 ))
unexpected_gene_row_samples="$(${AWK_BIN} -v expected="${EXPECTED_MATRIX_GENES}" -F'\t' 'NR > 1 && $2 != expected { n++ } END { print n + 0 }' "${STRANDNESS_SUMMARY}")"
strandness_missing_geneid_rows="$(${AWK_BIN} -F'\t' 'NR > 1 { n += $3 } END { print n + 0 }' "${STRANDNESS_SUMMARY}")"

if [[ "${strandness_samples}" -ne "${EXPECTED_SAMPLES}" ]]; then
    echo "[ERROR] Unexpected samples in strandness summary: ${strandness_samples}; expected ${EXPECTED_SAMPLES}" >&2
    exit 1
fi

if [[ "${unexpected_gene_row_samples}" -ne 0 ]]; then
    echo "[ERROR] Samples with unexpected GeneCounts gene-row count: ${unexpected_gene_row_samples}" >&2
    exit 1
fi

if [[ "${strandness_missing_geneid_rows}" -ne 0 ]]; then
    echo "[ERROR] MissingGeneID rows found across GeneCounts: ${strandness_missing_geneid_rows}" >&2
    exit 1
fi

log "[stage:strandness] strandness_summary=${STRANDNESS_SUMMARY}"
log "[stage:strandness] strandness_overall=${STRANDNESS_OVERALL}"


log "[stage:annotation] from GTF parse gene coordinates annotation"
${AWK_BIN} -F'\t' '
BEGIN {
    OFS = "\t"
    print "gene_id", "gene_name", "seqname", "start", "end", "strand", "midpoint", "source_feature"
}
$0 ~ /^#/ { next }
NF < 9 { next }
{
    gid = get_attr($9, "gene_id")
    if (gid == "") {
        next
    }

    gname = get_attr($9, "gene_name")
    if (gname == "") {
        gname = gid
    }

    if (!(gid in seen)) {
        seen[gid] = 1
        order[++n] = gid
        gene_name[gid] = gname
        seqname[gid] = $1
        start[gid] = $4
        end[gid] = $5
        strand[gid] = $7
        source_feature[gid] = $3
    } else {
        if ($4 < start[gid]) {
            start[gid] = $4
        }
        if ($5 > end[gid]) {
            end[gid] = $5
        }
        if (gene_name[gid] == gid && gname != gid) {
            gene_name[gid] = gname
        }
        if (source_feature[gid] != "gene" && $3 == "gene") {
            source_feature[gid] = $3
        }
    }
}
END {
    for (i = 1; i <= n; i++) {
        gid = order[i]
        midpoint = int((start[gid] + end[gid]) / 2)
        print gid, gene_name[gid], seqname[gid], start[gid], end[gid], strand[gid], midpoint, source_feature[gid]
    }
}
function get_attr(attr, key, pieces, value_parts, n) {
    n = split(attr, pieces, key " \"")
    if (n < 2) {
        return ""
    }
    split(pieces[2], value_parts, "\"")
    return value_parts[1]
}
' "${MANUAL_GTF}" > "${GENE_ANNOTATION}"
log "[stage:annotation] gene_annotation=${GENE_ANNOTATION}"


log "[stage:merge] merge ${expected_samples} samples of STAR reverse-stranded counts"
first_sample="$(${AWK_BIN} -F'\t' 'NR == 2 { print $1 }' "${SAMPLES_FILE}")"
first_counts="$(${AWK_BIN} -F'\t' 'NR == 2 { print $5 }' "${SAMPLES_FILE}")"
gene_list_tmp="$(${MKTEMP_BIN})"
matrix_body_tmp="$(${MKTEMP_BIN})"
header_tmp="$(${MKTEMP_BIN})"

${AWK_BIN} -F'\t' '$1 !~ /^N_/ { print $1 }' "${first_counts}" > "${gene_list_tmp}"

${TAIL_BIN} -n +2 "${SAMPLES_FILE}" | ${CUT_BIN} -f1 | ${AWK_BIN} '
BEGIN {
    printf "gene_id"
}
{
    printf "\t%s", $1
}
END {
    printf "\n"
}
' > "${header_tmp}"

${CP_BIN} "${gene_list_tmp}" "${matrix_body_tmp}"

${TAIL_BIN} -n +2 "${SAMPLES_FILE}" | while IFS=$'\t' read -r sample_id _animal_id _term _raw_data_name counts_file _bam_file _bai_file; do
    log "[stage:merge] process sample ${sample_id}: ${counts_file}"
    sample_genes_tmp="$(${MKTEMP_BIN})"
    sample_counts_tmp="$(${MKTEMP_BIN})"

    ${AWK_BIN} -F'\t' '$1 !~ /^N_/ { print $1 }' "${counts_file}" > "${sample_genes_tmp}"
    if ! ${CMP_BIN} -s "${gene_list_tmp}" "${sample_genes_tmp}"; then
        echo "[ERROR] Gene order/content differs for sample: ${sample_id}" >&2
        echo "[ERROR] Compare against first sample: ${first_sample}" >&2
        ${RM_BIN} -f "${sample_genes_tmp}" "${sample_counts_tmp}"
        exit 1
    fi

    ${AWK_BIN} -v col="${COUNT_COLUMN}" -F'\t' '$1 !~ /^N_/ { print $col }' "${counts_file}" > "${sample_counts_tmp}"
    ${PASTE_BIN} "${matrix_body_tmp}" "${sample_counts_tmp}" > "${matrix_body_tmp}.new"
    ${MV_BIN} "${matrix_body_tmp}.new" "${matrix_body_tmp}"
    ${RM_BIN} -f "${sample_genes_tmp}" "${sample_counts_tmp}"
done

${CAT_BIN} "${header_tmp}" "${matrix_body_tmp}" > "${MATRIX_FILE}"
log "[stage:merge] matrix_file=${MATRIX_FILE}"


log "[stage:final-check] check matrix gene and GTF coordinates annotation of consistency"
${AWK_BIN} -F'\t' 'NR == FNR { annotated[$1] = 1; next } NR > 1 && !($1 in annotated) { print $1 }' \
    "${GENE_ANNOTATION}" "${MATRIX_FILE}" > "${MISSING_ANNOTATION_GENES}"

matrix_genes=$(( $(${WC_BIN} -l < "${MATRIX_FILE}") - 1 ))
annotation_genes=$(( $(${WC_BIN} -l < "${GENE_ANNOTATION}") - 1 ))
missing_annotation_genes=$(${WC_BIN} -l < "${MISSING_ANNOTATION_GENES}")
duplicate_matrix_genes=$(${AWK_BIN} -F'\t' 'NR > 1 { seen[$1]++ } END { for (g in seen) if (seen[g] > 1) dup++; print dup + 0 }' "${MATRIX_FILE}")
matrix_missing_geneid_rows=$(${AWK_BIN} -F'\t' 'NR > 1 && $1 == "MissingGeneID" { n++ } END { print n + 0 }' "${MATRIX_FILE}")
matrix_gene_loc_rows=$(${AWK_BIN} -F'\t' 'NR > 1 && $1 ~ /^gene-LOC/ { n++ } END { print n + 0 }' "${MATRIX_FILE}")
matrix_dog5_new_gene_rows=$(${AWK_BIN} -F'\t' 'NR > 1 && $1 ~ /^DOG5_NEWGENE_/ { n++ } END { print n + 0 }' "${MATRIX_FILE}")
matrix_unified_new_gene_rows=$(${AWK_BIN} -F'\t' 'NR > 1 && $1 ~ /^UNIFIED_NEWGENE_/ { n++ } END { print n + 0 }' "${MATRIX_FILE}")
annotation_gene_loc_rows=$(${AWK_BIN} -F'\t' 'NR > 1 && $1 ~ /^gene-LOC/ { n++ } END { print n + 0 }' "${GENE_ANNOTATION}")
annotation_dog5_new_gene_rows=$(${AWK_BIN} -F'\t' 'NR > 1 && $1 ~ /^DOG5_NEWGENE_/ { n++ } END { print n + 0 }' "${GENE_ANNOTATION}")
annotation_unified_new_gene_rows=$(${AWK_BIN} -F'\t' 'NR > 1 && $1 ~ /^UNIFIED_NEWGENE_/ { n++ } END { print n + 0 }' "${GENE_ANNOTATION}")
annotation_only_genes=$(( annotation_genes - matrix_genes ))

if [[ "${matrix_genes}" -ne "${EXPECTED_MATRIX_GENES}" ]]; then
    echo "[ERROR] Unexpected matrix gene count: ${matrix_genes}; expected ${EXPECTED_MATRIX_GENES}" >&2
    exit 1
fi

if [[ "${annotation_genes}" -ne "${EXPECTED_ANNOTATION_GENES}" ]]; then
    echo "[ERROR] Unexpected annotation gene count: ${annotation_genes}; expected ${EXPECTED_ANNOTATION_GENES}" >&2
    exit 1
fi

if [[ "${matrix_gene_loc_rows}" -ne "${EXPECTED_MATRIX_GENE_LOC_ROWS}" ]]; then
    echo "[ERROR] Unexpected matrix gene-LOC rows: ${matrix_gene_loc_rows}; expected ${EXPECTED_MATRIX_GENE_LOC_ROWS}" >&2
    exit 1
fi

if [[ "${annotation_gene_loc_rows}" -ne "${EXPECTED_ANNOTATION_GENE_LOC_ROWS}" ]]; then
    echo "[ERROR] Unexpected annotation gene-LOC rows: ${annotation_gene_loc_rows}; expected ${EXPECTED_ANNOTATION_GENE_LOC_ROWS}" >&2
    exit 1
fi

if [[ "${matrix_dog5_new_gene_rows}" -ne "${EXPECTED_DOG5_NEW_GENE_ROWS}" || "${annotation_dog5_new_gene_rows}" -ne "${EXPECTED_DOG5_NEW_GENE_ROWS}" ]]; then
    echo "[ERROR] Unexpected DOG5_NEWGENE rows: matrix=${matrix_dog5_new_gene_rows}; annotation=${annotation_dog5_new_gene_rows}; expected ${EXPECTED_DOG5_NEW_GENE_ROWS}" >&2
    exit 1
fi

if [[ "${matrix_unified_new_gene_rows}" -ne 0 || "${annotation_unified_new_gene_rows}" -ne 0 ]]; then
    echo "[ERROR] five-only unexpected occurrence in the results UNIFIED_NEWGENE row: matrix=${matrix_unified_new_gene_rows}; annotation=${annotation_unified_new_gene_rows}" >&2
    exit 1
fi

if [[ "${annotation_only_genes}" -ne "${EXPECTED_ANNOTATION_ONLY_GENES}" ]]; then
    echo "[ERROR] Unexpected annotation-only gene count: ${annotation_only_genes}; expected ${EXPECTED_ANNOTATION_ONLY_GENES}" >&2
    exit 1
fi

if [[ "${duplicate_matrix_genes}" -ne 0 ]]; then
    echo "[ERROR] Duplicate gene_id rows found in matrix: ${duplicate_matrix_genes}" >&2
    exit 1
fi

if [[ "${missing_annotation_genes}" -ne 0 ]]; then
    echo "[ERROR] Some matrix genes have no coordinate annotation in Manual-annotation GTF." >&2
    echo "[ERROR] See: ${MISSING_ANNOTATION_GENES}" >&2
    exit 1
fi

if [[ "${matrix_missing_geneid_rows}" -ne 0 ]]; then
    echo "[ERROR] Matrix still contains MissingGeneID rows, which have no reliable locus coordinate." >&2
    echo "[ERROR] Please check whether STAR index was built with the Manual-annotation GTF." >&2
    exit 1
fi

{
    echo "status=OK"
    echo "run_id=${RUN_ID}"
    echo "project_dir=${PROJECT_DIR}"
    echo "star_run_id=${STAR_RUN_ID}"
    echo "star_dir=${STAR_DIR}"
    echo "out_dir=${OUT_DIR}"
    echo "sample_tsv=${SAMPLE_TSV}"
    echo "manual_gtf=${MANUAL_GTF}"
    echo "expected_samples=${expected_samples}"
    echo "duplicate_sample_ids=${duplicate_sample_ids}"
    echo "missing_star_outputs=${missing_outputs}"
    echo "selected_count_column=${COUNT_COLUMN}"
    echo "selected_count_type=${COUNT_COLUMN_NAME}"
    echo "matrix_file=${MATRIX_FILE}"
    echo "matrix_genes=${matrix_genes}"
    echo "gene_LOC_rows_kept=${matrix_gene_loc_rows}"
    echo "dog5_new_gene_rows_kept=${matrix_dog5_new_gene_rows}"
    echo "legacy_unified_new_gene_rows=${matrix_unified_new_gene_rows}"
    echo "MissingGeneID_rows=${matrix_missing_geneid_rows}"
    echo "annotation_file=${GENE_ANNOTATION}"
    echo "annotation_genes=${annotation_genes}"
    echo "annotation_gene_LOC_rows=${annotation_gene_loc_rows}"
    echo "annotation_dog5_new_gene_rows=${annotation_dog5_new_gene_rows}"
    echo "annotation_legacy_unified_new_gene_rows=${annotation_unified_new_gene_rows}"
    echo "annotation_only_genes=${annotation_only_genes}"
    echo "matrix_genes_missing_coordinate_annotation=${missing_annotation_genes}"
    echo "duplicate_matrix_gene_ids=${duplicate_matrix_genes}"
    echo "strandness_summary=${STRANDNESS_SUMMARY}"
    echo "strandness_overall=${STRANDNESS_OVERALL}"
} > "${CHECK_LOG}"

log "[summary] matrix_genes=${matrix_genes}; annotation_genes=${annotation_genes}; annotation_only_genes=${annotation_only_genes}"
log "[summary] gene_LOC_rows_kept=${matrix_gene_loc_rows}; dog5_new_gene_rows_kept=${matrix_dog5_new_gene_rows}; legacy_unified_new_gene_rows=${matrix_unified_new_gene_rows}"
log "[summary] MissingGeneID_rows=${matrix_missing_geneid_rows}; matrix_genes_missing_coordinate_annotation=${missing_annotation_genes}; duplicate_matrix_gene_ids=${duplicate_matrix_genes}"
printf 'run_id\t%s\nexecution_status\tsuccess\ncompleted_at\t%s\n' \
    "${RUN_ID}" "$(${DATE_BIN} '+%F %T')" > "${COMPLETION_MARKER}"
log "[done] COMPLETION_MARKER=${COMPLETION_MARKER}"
log "[done] finemaping expression raw-count matrix build completed; task completed / done"
log "execution_status=success"
log "reminder: execution completion does not establish output completeness, format, count or biological quality already passed; please in run after separately check result."

echo "[OK] Expression matrix generated:"
echo "     ${MATRIX_FILE}"
echo "[OK] Gene coordinate annotation generated:"
echo "     ${GENE_ANNOTATION}"
echo "[OK] Check log:"
echo "     ${CHECK_LOG}"
