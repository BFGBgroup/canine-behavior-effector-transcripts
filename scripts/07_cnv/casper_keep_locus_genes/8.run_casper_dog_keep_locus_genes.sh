#!/usr/bin/env bash
set -euo pipefail


TISSUE="${1:-amygdala}"
CNV_SCALE="${2:-${CNV_SCALE:-3}}"
LOH_SCALE="${3:-${LOH_SCALE:-3}}"
BAF_MODE="${4:-${BAF_MODE:-hybrid}}"
EXTRA_ARGS=("${@:5}")

INPUT_PROJECT_DIR="/path/to/project/legacy/five-dogs/5.cnv"
OUTPUT_ROOT="/path/to/project/analysis_20260609/5-dogs/4.cnv-rna/8.casper_dog_keep_locus_genes-v2"
SCRIPT_DIR="/path/to/repository/scripts/07_cnv/casper_keep_locus_genes"

RENDERED_R="${SCRIPT_DIR}/8.run_casper_dog_keep_locus_genes.rendered.R"
CONDA_ENV="CaSpER"
CHECK_ONLY="${CHECK_ONLY:-0}"

EXPR_FILE="${INPUT_PROJECT_DIR}/5.casper_input/bulk_rna.cnv_expression.keep_locus_genes.all_genes.log2_norm.primary_chr.casper_ordered.tsv"
ANNOT_FILE="${INPUT_PROJECT_DIR}/5.casper_input/CanFam3.1.keep_locus_genes.all_genes.gene_annotation.primary_chr.for_casper.tsv"
SAMPLE_META="${INPUT_PROJECT_DIR}/5.casper_input/bulk_rna.samples.keep_locus_genes.all_genes.for_casper.tsv"
BAF_DIR="${INPUT_PROJECT_DIR}/6.baf_signal/keep_locus_genes.all_genes/primary/per_sample_baf"
BAF_QC="${INPUT_PROJECT_DIR}/6.baf_signal/bulk_rna.keep_locus_genes.all_genes.baf_from_ase_counts.qc.tsv"
LOH_MAPPING="${INPUT_PROJECT_DIR}/6.baf_signal/bulk_rna.keep_locus_genes.all_genes.casper_loh_name_mapping.tsv"
OUT_ROOT="${OUTPUT_ROOT}"
OUT_TISSUE="${OUT_ROOT}/${TISSUE}"
SAMPLE_MODE_FILE="${OUT_TISSUE}/sample_baf_mode.tsv"

die() {
  echo "[ERROR] $*" >&2
  exit 1
}

require_file() {
  local path="$1"
  local label="$2"
  [[ -s "${path}" ]] || die "${label} not found or empty: ${path}"
}

require_dir() {
  local path="$1"
  local label="$2"
  [[ -d "${path}" ]] || die "${label} not found: ${path}"
}

require_file "${RENDERED_R}" "Local CaSpER R script"
require_file "${EXPR_FILE}" "CaSpER expression matrix"
require_file "${ANNOT_FILE}" "CaSpER gene annotation"
require_file "${SAMPLE_META}" "CaSpER sample metadata"
require_file "${BAF_QC}" "BAF QC table"
require_file "${LOH_MAPPING}" "LOH name mapping"
require_dir "${BAF_DIR}" "BAF primary directory"

expr_samples="$(awk -F '\t' 'NR == 1 { print NF - 1; exit }' "${EXPR_FILE}")"
expr_genes="$(awk 'END { n = NR - 1; if (n < 0) n = 0; print n }' "${EXPR_FILE}")"
annot_genes="$(awk 'END { n = NR - 1; if (n < 0) n = 0; print n }' "${ANNOT_FILE}")"
meta_samples="$(awk 'END { n = NR - 1; if (n < 0) n = 0; print n }' "${SAMPLE_META}")"
baf_files="$(find "${BAF_DIR}" -maxdepth 1 -name '*.baf.primary.tsv' -type f | wc -l)"

[[ "${expr_samples}" -gt 0 ]] || die "No sample columns found in expression matrix"
[[ "${expr_genes}" -gt 0 ]] || die "No genes found in expression matrix"
[[ "${expr_genes}" -eq "${annot_genes}" ]] || die "Expression gene count (${expr_genes}) != annotation gene count (${annot_genes})"
[[ "${expr_samples}" -eq "${meta_samples}" ]] || die "Expression sample count (${expr_samples}) != sample metadata count (${meta_samples})"

if ! diff -q \
  <(awk -F '\t' 'NR > 1 { print $1 }' "${EXPR_FILE}") \
  <(awk -F '\t' 'NR > 1 { print $1 }' "${ANNOT_FILE}") \
  >/dev/null; then
  die "Expression matrix and annotation are not in the same gene order"
fi

if ! diff -q \
  <(awk -F '\t' 'NR == 1 { for (i = 2; i <= NF; i++) print $i }' "${EXPR_FILE}") \
  <(awk -F '\t' 'NR > 1 { print $1 }' "${SAMPLE_META}") \
  >/dev/null; then
  die "Expression sample columns and sample metadata sample_id are not in the same order"
fi

missing_baf_all="$(
  awk -F '\t' 'NR == 1 { for (i = 2; i <= NF; i++) print $i }' "${EXPR_FILE}" |
    while read -r sample_id; do
      [[ -s "${BAF_DIR}/${sample_id}.baf.primary.tsv" ]] || echo "${sample_id}"
    done
)"
if [[ -n "${missing_baf_all}" ]]; then
  echo "${missing_baf_all}" >&2
  die "Missing BAF file(s) for expression samples"
fi

tissue_samples="$(
  awk -F '\t' -v tissue="${TISSUE}" 'NR == 1 {
    suffix = "_" tissue
    for (i = 2; i <= NF; i++) {
      if (length($i) >= length(suffix) &&
          substr($i, length($i) - length(suffix) + 1) == suffix) {
        print $i
      }
    }
  }' "${EXPR_FILE}"
)"
tissue_sample_count="$(printf '%s\n' "${tissue_samples}" | awk 'NF > 0 { n++ } END { print n + 0 }')"
[[ "${tissue_sample_count}" -gt 0 ]] || die "No expression samples found for tissue=${TISSUE}; five-dogs expects sample_id suffix _${TISSUE}"

missing_baf_for_tissue="$(
  printf '%s\n' "${tissue_samples}" |
    awk 'NF > 0' |
    while read -r sample_id; do
      [[ -s "${BAF_DIR}/${sample_id}.baf.primary.tsv" ]] || echo "${sample_id}"
    done
)"
if [[ -n "${missing_baf_for_tissue}" ]]; then
  echo "${missing_baf_for_tissue}" >&2
  die "Missing BAF file(s) for tissue=${TISSUE}"
fi

first_baf="$(find "${BAF_DIR}" -maxdepth 1 -name '*.baf.primary.tsv' -type f -print -quit)"
require_file "${first_baf}" "Example BAF file"
missing_baf_cols="$(
  awk -F '\t' 'NR == 1 {
    need["chr"] = 1
    need["position"] = 1
    need["total_count"] = 1
    need["baf"] = 1
    for (i = 1; i <= NF; i++) delete need[$i]
    for (k in need) print k
  }' "${first_baf}"
)"
[[ -z "${missing_baf_cols}" ]] || die "BAF file missing required column(s): ${missing_baf_cols}"

unusable_baf_for_tissue="$(
  printf '%s\n' "${tissue_samples}" |
    awk 'NF > 0' |
    while read -r sample_id; do
      baf_file="${BAF_DIR}/${sample_id}.baf.primary.tsv"
      usable_sites="$(
        awk -F '\t' '
          NR == 1 {
            for (i = 1; i <= NF; i++) h[$i] = i
            next
          }
          $h["chr"] ~ /^chr([1-9]|[1-2][0-9]|3[0-8])$/ &&
            ($(h["total_count"]) + 0) >= 20 &&
            ($(h["baf"]) + 0) > 0.2 &&
            ($(h["baf"]) + 0) < 0.8 {
            n++
          }
          END { print n + 0 }
        ' "${baf_file}"
      )"
      [[ "${usable_sites}" -gt 0 ]] || echo "${sample_id}"
    done
)"
if [[ -n "${unusable_baf_for_tissue}" && "${BAF_MODE}" != "off" ]]; then
  echo "[WARN] sample-level BAF fallback will be used for tissue=${TISSUE}; unusable sample(s):" >&2
  echo "${unusable_baf_for_tissue}" >&2
fi

mkdir -p "${OUT_TISSUE}"

{
  printf 'sample_id\trequested_baf_mode\teffective_sample_mode\tusable_autosomal_baf_sites\tnote\n'
  printf '%s\n' "${tissue_samples}" |
    awk 'NF > 0' |
    while read -r sample_id; do
      baf_file="${BAF_DIR}/${sample_id}.baf.primary.tsv"
      usable_sites="$(
        awk -F '\t' '
          NR == 1 {
            for (i = 1; i <= NF; i++) h[$i] = i
            next
          }
          $h["chr"] ~ /^chr([1-9]|[1-2][0-9]|3[0-8])$/ &&
            ($(h["total_count"]) + 0) >= 20 &&
            ($(h["baf"]) + 0) > 0.2 &&
            ($(h["baf"]) + 0) < 0.8 {
            n++
          }
          END { print n + 0 }
        ' "${baf_file}"
      )"

      effective_mode="hybrid"
      note="baf_used"

      if [[ "${BAF_MODE}" == "off" ]]; then
        effective_mode="expression_only"
        note="baf_mode_off"
      elif [[ "${usable_sites}" -eq 0 ]]; then
        effective_mode="expression_only_sample_fallback"
        note="no_usable_baf_sites; sample_fallback_to_expression_only"
      fi

      printf '%s\t%s\t%s\t%s\t%s\n' \
        "${sample_id}" "${BAF_MODE}" "${effective_mode}" "${usable_sites}" "${note}"
    done
} > "${SAMPLE_MODE_FILE}"

chmod +x "${RENDERED_R}"

echo "=========================================="
echo "Dog CaSpER CNV Workflow: five-dogs keep_locus_genes"
echo "=========================================="
echo "Tissue:        ${TISSUE}"
echo "CNV scale:     ${CNV_SCALE}"
echo "LOH scale:     ${LOH_SCALE}"
echo "BAF mode:      ${BAF_MODE}"
echo "Conda env:     ${CONDA_ENV}"
echo "Expression:    ${EXPR_FILE}"
echo "Annotation:    ${ANNOT_FILE}"
echo "Sample meta:   ${SAMPLE_META}"
echo "BAF dir:       ${BAF_DIR}"
echo "BAF QC:        ${BAF_QC}"
echo "LOH mapping:   ${LOH_MAPPING}"
echo "Output root:   ${OUT_ROOT}"
echo "Sample modes:  ${SAMPLE_MODE_FILE}"
echo "Rendered R:    ${RENDERED_R}"
echo "Expr samples:  ${expr_samples}"
echo "Expr genes:    ${expr_genes}"
echo "BAF files:     ${baf_files}"
echo "Tissue samples:${tissue_sample_count}"
if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
  echo "Extra args:    ${EXTRA_ARGS[*]}"
fi
echo ""

if [[ "${CHECK_ONLY}" == "1" ]]; then
  echo "[CHECK_ONLY] input checks passed; CaSpER was not run."
  exit 0
fi

if ! command -v conda >/dev/null 2>&1; then
  die "conda not found in PATH"
fi

mkdir -p "${OUT_ROOT}"

set +e
conda run -n "${CONDA_ENV}" Rscript "${RENDERED_R}" \
  --tissue "${TISSUE}" \
  --cnv-scale "${CNV_SCALE}" \
  --loh-scale "${LOH_SCALE}" \
  --baf-mode "${BAF_MODE}" \
  "${EXTRA_ARGS[@]}"
EXIT_CODE=$?
set -e

if [[ "${EXIT_CODE}" -eq 0 ]]; then
  echo ""
  echo "=========================================="
  echo "Analysis completed successfully"
  echo "Output:"
  echo "${OUT_ROOT}/${TISSUE}"
  echo "=========================================="
else
  echo ""
  echo "=========================================="
  echo "Analysis failed with exit code ${EXIT_CODE}"
  echo "Log:"
  echo "${OUT_ROOT}/${TISSUE}/7.run_casper_dog.${TISSUE}.log"
  echo "=========================================="
  exit "${EXIT_CODE}"
fi
