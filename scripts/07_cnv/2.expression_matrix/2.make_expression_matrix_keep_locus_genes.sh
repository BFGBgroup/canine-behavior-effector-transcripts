#!/usr/bin/env bash
set -euo pipefail


PROJECT_DIR="/path/to/project/legacy/five-dogs/5.cnv"
STAR_DIR="${PROJECT_DIR}/1.star"
OUT_DIR="${PROJECT_DIR}/2.expression_matrix"

SAMPLE_TSV="/path/to/project/legacy/five-dogs/data/statistics/five_dogs_fastq_pair_stats.tsv"
PATCHED_GTF="${PROJECT_DIR}/0.gtf/GCF_000002285.3_CanFam3.1_genomic_chr.gene_id_patched.gtf"

SAMPLES_FILE="${OUT_DIR}/bulk_rna.star.samples.tsv"
MISSING_OUTPUTS_FILE="${OUT_DIR}/bulk_rna.star.missing_outputs.tsv"
STRANDNESS_SUMMARY="${OUT_DIR}/star_gene_counts_strandness_summary.tsv"
STRANDNESS_OVERALL="${OUT_DIR}/star_gene_counts_strandness_overall.txt"
GENE_ANNOTATION="${OUT_DIR}/CanFam3.1.patched_gene_annotation.for_expression_matrix.tsv"
MATRIX_FILE="${OUT_DIR}/bulk_rna.star.reverse_stranded.raw_counts.keep_locus_genes.tsv"
MISSING_ANNOTATION_GENES="${OUT_DIR}/bulk_rna.star.genes_missing_from_patched_gtf.tsv"
CHECK_LOG="${OUT_DIR}/bulk_rna.star.reverse_stranded.merge.keep_locus_genes.check.log"


COUNT_COLUMN=4
COUNT_COLUMN_NAME="reverse_stranded"

mkdir -p "${OUT_DIR}"

if [[ ! -s "${SAMPLE_TSV}" ]]; then
    echo "[ERROR] sample TSV not found or empty: ${SAMPLE_TSV}" >&2
    exit 1
fi

if [[ ! -d "${STAR_DIR}" ]]; then
    echo "[ERROR] STAR output directory not found: ${STAR_DIR}" >&2
    exit 1
fi

if [[ ! -s "${PATCHED_GTF}" ]]; then
    echo "[ERROR] patched GTF not found or empty: ${PATCHED_GTF}" >&2
    echo "[ERROR] This script intentionally does not fall back to the original GTF." >&2
    exit 1
fi


awk -v star_dir="${STAR_DIR}" -F'\t' '
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

expected_samples=$(( $(wc -l < "${SAMPLES_FILE}") - 1 ))
if [[ "${expected_samples}" -le 0 ]]; then
    echo "[ERROR] no samples parsed from ${SAMPLE_TSV}" >&2
    exit 1
fi


{
    echo -e "sample_id\tanimal_id\tterm\tmissing_item\tpath"
    tail -n +2 "${SAMPLES_FILE}" | while IFS=$'\t' read -r sample_id animal_id term _raw_data_name counts_file bam_file bai_file; do

        [[ -s "${counts_file}" ]] || echo -e "${sample_id}\t${animal_id}\t${term}\tReadsPerGene.out.tab\t${counts_file}"
        [[ -s "${bam_file}" ]] || echo -e "${sample_id}\t${animal_id}\t${term}\tsorted_bam\t${bam_file}"
        [[ -s "${bai_file}" ]] || echo -e "${sample_id}\t${animal_id}\t${term}\tbam_index\t${bai_file}"
    done
} > "${MISSING_OUTPUTS_FILE}"

missing_outputs=$(( $(wc -l < "${MISSING_OUTPUTS_FILE}") - 1 ))
if [[ "${missing_outputs}" -ne 0 ]]; then
    echo "[ERROR] Missing STAR outputs were found. See: ${MISSING_OUTPUTS_FILE}" >&2
    exit 1
fi


{
    echo -e "sample_id\tgene_rows\tmissing_geneid_rows\tunstranded_sum\tforward_sum\treverse_sum\treverse_over_forward"
    tail -n +2 "${SAMPLES_FILE}" | while IFS=$'\t' read -r sample_id _animal_id _term _raw_data_name counts_file _bam_file _bai_file; do
        awk -v sample_id="${sample_id}" -F'\t' '
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

awk -F'\t' '
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


awk -F'\t' '
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
' "${PATCHED_GTF}" > "${GENE_ANNOTATION}"


first_sample="$(awk -F'\t' 'NR == 2 { print $1 }' "${SAMPLES_FILE}")"
first_counts="$(awk -F'\t' 'NR == 2 { print $5 }' "${SAMPLES_FILE}")"
gene_list_tmp="$(mktemp)"
matrix_body_tmp="$(mktemp)"
header_tmp="$(mktemp)"
trap 'rm -f "${gene_list_tmp}" "${matrix_body_tmp}" "${header_tmp}"' EXIT

awk -F'\t' '$1 !~ /^N_/ { print $1 }' "${first_counts}" > "${gene_list_tmp}"

tail -n +2 "${SAMPLES_FILE}" | cut -f1 | awk '
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

cp "${gene_list_tmp}" "${matrix_body_tmp}"

tail -n +2 "${SAMPLES_FILE}" | while IFS=$'\t' read -r sample_id _animal_id _term _raw_data_name counts_file _bam_file _bai_file; do
    sample_genes_tmp="$(mktemp)"
    sample_counts_tmp="$(mktemp)"

    awk -F'\t' '$1 !~ /^N_/ { print $1 }' "${counts_file}" > "${sample_genes_tmp}"
    if ! cmp -s "${gene_list_tmp}" "${sample_genes_tmp}"; then
        echo "[ERROR] Gene order/content differs for sample: ${sample_id}" >&2
        echo "[ERROR] Compare against first sample: ${first_sample}" >&2
        rm -f "${sample_genes_tmp}" "${sample_counts_tmp}"
        exit 1
    fi

    awk -v col="${COUNT_COLUMN}" -F'\t' '$1 !~ /^N_/ { print $col }' "${counts_file}" > "${sample_counts_tmp}"
    paste "${matrix_body_tmp}" "${sample_counts_tmp}" > "${matrix_body_tmp}.new"
    mv "${matrix_body_tmp}.new" "${matrix_body_tmp}"
    rm -f "${sample_genes_tmp}" "${sample_counts_tmp}"
done

cat "${header_tmp}" "${matrix_body_tmp}" > "${MATRIX_FILE}"


awk -F'\t' 'NR == FNR { annotated[$1] = 1; next } NR > 1 && !($1 in annotated) { print $1 }' \
    "${GENE_ANNOTATION}" "${MATRIX_FILE}" > "${MISSING_ANNOTATION_GENES}"

matrix_genes=$(( $(wc -l < "${MATRIX_FILE}") - 1 ))
annotation_genes=$(( $(wc -l < "${GENE_ANNOTATION}") - 1 ))
missing_annotation_genes=$(wc -l < "${MISSING_ANNOTATION_GENES}")
duplicate_matrix_genes=$(awk -F'\t' 'NR > 1 { seen[$1]++ } END { for (g in seen) if (seen[g] > 1) dup++; print dup + 0 }' "${MATRIX_FILE}")
matrix_missing_geneid_rows=$(awk -F'\t' 'NR > 1 && $1 == "MissingGeneID" { n++ } END { print n + 0 }' "${MATRIX_FILE}")
matrix_gene_loc_rows=$(awk -F'\t' 'NR > 1 && $1 ~ /^gene-LOC/ { n++ } END { print n + 0 }' "${MATRIX_FILE}")

if [[ "${duplicate_matrix_genes}" -ne 0 ]]; then
    echo "[ERROR] Duplicate gene_id rows found in matrix: ${duplicate_matrix_genes}" >&2
    exit 1
fi

if [[ "${missing_annotation_genes}" -ne 0 ]]; then
    echo "[ERROR] Some matrix genes have no coordinate annotation in patched GTF." >&2
    echo "[ERROR] See: ${MISSING_ANNOTATION_GENES}" >&2
    exit 1
fi

if [[ "${matrix_missing_geneid_rows}" -ne 0 ]]; then
    echo "[ERROR] Matrix still contains MissingGeneID rows, which have no reliable locus coordinate." >&2
    echo "[ERROR] Please check whether STAR index was built with the patched GTF." >&2
    exit 1
fi

{
    echo "status=OK"
    echo "project_dir=${PROJECT_DIR}"
    echo "star_dir=${STAR_DIR}"
    echo "out_dir=${OUT_DIR}"
    echo "sample_tsv=${SAMPLE_TSV}"
    echo "patched_gtf=${PATCHED_GTF}"
    echo "expected_samples=${expected_samples}"
    echo "missing_star_outputs=${missing_outputs}"
    echo "selected_count_column=${COUNT_COLUMN}"
    echo "selected_count_type=${COUNT_COLUMN_NAME}"
    echo "matrix_file=${MATRIX_FILE}"
    echo "matrix_genes=${matrix_genes}"
    echo "gene_LOC_rows_kept=${matrix_gene_loc_rows}"
    echo "MissingGeneID_rows=${matrix_missing_geneid_rows}"
    echo "annotation_file=${GENE_ANNOTATION}"
    echo "annotation_genes=${annotation_genes}"
    echo "matrix_genes_missing_coordinate_annotation=${missing_annotation_genes}"
    echo "duplicate_matrix_gene_ids=${duplicate_matrix_genes}"
    echo "strandness_summary=${STRANDNESS_SUMMARY}"
    echo "strandness_overall=${STRANDNESS_OVERALL}"
} > "${CHECK_LOG}"

echo "[OK] Expression matrix generated:"
echo "     ${MATRIX_FILE}"
echo "[OK] Gene coordinate annotation generated:"
echo "     ${GENE_ANNOTATION}"
echo "[OK] Check log:"
echo "     ${CHECK_LOG}"
