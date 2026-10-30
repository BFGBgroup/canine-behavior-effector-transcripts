#!/usr/bin/bash
set -euo pipefail


SCRIPT_DIR="/path/to/repository/scripts/12_evidence_integration/6.allelome-compare-eqtl-finemap/scripts/1.1_compare_dog5_allelome_eqtl_finemap"
LOG_DIR="${SCRIPT_DIR}/log"
MASTER_LOG="${LOG_DIR}/1.1_compare_dog5_allelome_eqtl_finemap.log"
LOCK_FILE="${LOG_DIR}/1.1_compare_dog5_allelome_eqtl_finemap.lock"
PYTHON_BIN="/usr/bin/python3.10"
PYTHON_SCRIPT="${SCRIPT_DIR}/1.1.compare_dog5_allelome_eqtl_finemap.py"
ALLELOME_LINKS="/path/to/project/analysis_20260817/5-dogs/7.Allelome.LINK/1.4_annotate_and_merge_dog5_allelome_links/20260907_213429/dog5_allelome_links.annotated.tsv"
ALLELOME_SUMMARY="/path/to/project/analysis_20260817/5-dogs/7.Allelome.LINK/1.4_annotate_and_merge_dog5_allelome_links/20260907_213429/dog5_allelome_annotation.summary.tsv"
CANDIDATE_BLOCKS="/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.2_no_gwas_prepare_expression_ase_plasma_inputs/20260821_181104/no_gwas_plasma_candidate_blocks.tsv"
VALID_BLOCKS="/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.5_no_gwas_build_expression_ase_regulatory_catalogue/20260827_201055/no_gwas_expression_ase_valid_block_catalogue.tsv"
CREDIBLE_VARIANTS="/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.6_no_gwas_export_expression_ase_catalogue_deliverables/20260827_224912/no_gwas_expression_ase_credible_set_variants.tsv.gz"
FINAL_QC_SUMMARY="/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.7_no_gwas_final_catalogue_qc_report/20260828_122551/no_gwas_expression_ase_final_qc_summary.tsv"
NEW_GENE_GTF="/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.new_only.gtf"
NEW_GENE_TRANSCRIPT_MAP="/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.transcript_map.tsv"
NEW_GENE_BUILD_SUMMARY="/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.summary.tsv"
OUTPUT_ROOT="/path/to/project/analysis_20260817/results/5-dogs/6.allelome-compare-eqtl-finemap/1.1_compare_dog5_allelome_eqtl_finemap"
RUN_ID="${RUN_ID:-$(/usr/bin/date '+%Y%m%d_%H%M%S')}"

/usr/bin/mkdir -p "${LOG_DIR}"
exec 9>"${LOCK_FILE}"
if ! /usr/bin/flock -n 9; then
  printf '[%s ] execution_status=failed; existing the same step task holds run lock: %s\n' \
    "$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')" "${LOCK_FILE}" >&2
  exit 75
fi
: > "${MASTER_LOG}"
exec >> "${MASTER_LOG}" 2>&1

printf '[%s ] prepare run directory: %s/%s\n' \
  "$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')" "${OUTPUT_ROOT}" "${RUN_ID}"
printf '[%s ] single task, single thread; canonical main log: %s\n' \
  "$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')" "${MASTER_LOG}"

SUCCESS_FILE="${LOG_DIR}/1.1_compare_dog5_allelome_eqtl_finemap.${RUN_ID}.success.tsv"
FAILED_FILE="${LOG_DIR}/1.1_compare_dog5_allelome_eqtl_finemap.${RUN_ID}.failed.tsv"
printf 'task\texit_code\n' > "${SUCCESS_FILE}"
printf 'task\texit_code\n' > "${FAILED_FILE}"
set +e
"${PYTHON_BIN}" "${PYTHON_SCRIPT}" \
  --allelome-links "${ALLELOME_LINKS}" \
  --allelome-summary "${ALLELOME_SUMMARY}" \
  --candidate-blocks "${CANDIDATE_BLOCKS}" \
  --valid-blocks "${VALID_BLOCKS}" \
  --credible-variants "${CREDIBLE_VARIANTS}" \
  --final-qc-summary "${FINAL_QC_SUMMARY}" \
  --new-gene-gtf "${NEW_GENE_GTF}" \
  --new-gene-transcript-map "${NEW_GENE_TRANSCRIPT_MAP}" \
  --new-gene-build-summary "${NEW_GENE_BUILD_SUMMARY}" \
  --output-root "${OUTPUT_ROOT}" \
  --run-id "${RUN_ID}" \
  --priority-ppa 0.1
exit_code=$?
set -e

if (( exit_code != 0 )); then
  printf '[%s ] execution_status=failed; Python main program exit code=%d\n' \
    "$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')" "${exit_code}"
  printf 'summary\t%s\n' "${exit_code}" >> "${FAILED_FILE}"
  exit "${exit_code}"
fi

printf 'summary\t0\n' >> "${SUCCESS_FILE}"
