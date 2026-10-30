#!/bin/bash
set -euo pipefail


PYTHON_BIN="/usr/bin/python3"
SCRIPT="/path/to/repository/scripts/03_candidate_catalog/2.2_candidate_TSS_promoter_support/2.2_candidate_TSS_promoter_support.py"
PROJECT_ROOT="/path/to/project/analysis_20260709/5-dogs/9.new-gene"
STEP21_ROOT="${PROJECT_ROOT}/2.1_transcript_first_preASE_candidate_assessment"
OUTPUT_ROOT="${PROJECT_ROOT}/2.2_candidate_TSS_promoter_support"
LOG_DIR="/path/to/repository/scripts/03_candidate_catalog/2.2_candidate_TSS_promoter_support/log"
RUN_ID="$(/usr/bin/date '+%Y%m%d_%H%M%S')"
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

RETAINED_TSV="${STEP21_ROOT}/2.1_retained_transcript_candidates.preASE.tsv"
SUMMARY_TSV="${STEP21_ROOT}/2.1_candidate_transcript_summary.preASE.tsv"
EXON_BED="/path/to/project/legacy/five-dogs/9.new-gene/1.6_gffcompare_stringtie3_merged_transcriptome-v6/transcript_exon.bed"
CHROM_SIZES="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna.fai"
TSS_ROOT="/path/to/project/legacy/five-dogs/9.new-gene/reference/TSS"
DOGA_CTSS_REGION_RESTRICTED="${TSS_ROOT}/3.DoGA-CanFam4-to-CanFam3-liftover/results/doga_CTSS_brain_CNS.CanFam3.bed"
DOGA_PROMOTER="${TSS_ROOT}/3.DoGA-CanFam4-to-CanFam3-liftover/results/doga_promoters_all.CanFam3.bed"
BARKBASE_ATAC="${TSS_ROOT}/4.BarkBase-ATAC-seq/results/barkbase_ATAC_all_selected.merged.bed"
EPICDOG_PROMOTER="${TSS_ROOT}/5.EpiC-Dog-chromatin-states/results/epicdog_brain_CNS_promoter_like_states.bed"
EPICDOG_ACTIVE="${TSS_ROOT}/5.EpiC-Dog-chromatin-states/results/epicdog_brain_CNS_active_H3K4me3_H3K27ac_states.bed"
EPICDOG_SUPER="${TSS_ROOT}/5.EpiC-Dog-chromatin-states/results/epicdog_brain_CNS_super_enhancers.bed"


/usr/bin/mkdir -p "$LOG_DIR"
RUN_LOG="${LOG_DIR}/${RUN_ID}.master.log"
SUCCESS_FILE="${LOG_DIR}/${RUN_ID}.success.tsv"
FAILED_FILE="${LOG_DIR}/${RUN_ID}.failed.tsv"
set -C
exec > "$RUN_LOG" 2>&1
printf 'task\trun_id\tstatus\texit_code\toutput_dir\n' > "$SUCCESS_FILE"
printf 'task\trun_id\tstatus\texit_code\toutput_dir\n' > "$FAILED_FILE"
set +C
CHILD_PID=""
finish() {
  local code=$?
  trap - EXIT
  if (( code != 0 )); then
    printf 'Step2.2\t%s\tfailed\t%s\t%s\n' "$RUN_ID" "$code" "$RUN_DIR" >> "$FAILED_FILE"
    printf 'execution_status=failed total=1 success=0 failed=1 skipped=0 exit_code=%s\n' "$code"
  fi
  exit "$code"
}
stop_worker() {
  local code="$1"
  if [[ -n "$CHILD_PID" ]]; then
    kill -TERM "$CHILD_PID" 2>/dev/null || true
    wait "$CHILD_PID" || true
  fi
  exit "$code"
}
trap finish EXIT
trap 'stop_worker 130' INT
trap 'stop_worker 143' TERM
printf '[%s ] START concurrent_tasks=1 threads_per_task=1(retain serial configuration)\n' "$(/usr/bin/date '+%F %T')"
printf 'script=%s\noutput=%s\n' "$SCRIPT" "$RUN_DIR"
printf 'proximity_metric=has_cage_proximity_support_500bp\nproximity_rule=same_strand_robust_promoter_distance_le_500bp\nstrict_position_support=reported_separately\nsensitivity_window_bp=100\n'

printf 'restricted_ctss=%s\nall_cage_sensitivity_analysis=not_performed\n' "$DOGA_CTSS_REGION_RESTRICTED"
printf 'all_robust_promoter=%s\npromoter_reference_scope=all_DoGA_robust_promoters_successfully_lifted\ndefault_reference_combination=region_restricted_CTSS_plus_all_robust_promoters\n' "$DOGA_PROMOTER"
for required_file in "$PYTHON_BIN" "$SCRIPT" "$RETAINED_TSV" "$SUMMARY_TSV" "$EXON_BED" "$CHROM_SIZES" "$DOGA_CTSS_REGION_RESTRICTED" "$DOGA_PROMOTER" "$BARKBASE_ATAC" "$EPICDOG_PROMOTER" "$EPICDOG_ACTIVE" "$EPICDOG_SUPER"; do
  if [[ ! -s "$required_file" ]]; then
    echo "[ERROR] required input missing or empty: $required_file" >&2
    exit 1
  fi
done
[[ -x "$PYTHON_BIN" ]] || { echo '[ERROR] Python not executable'; exit 1; }
[[ ! -e "$RUN_DIR" ]] || { echo '[ERROR] timestamp directory already exists'; exit 1; }
"$PYTHON_BIN" "$SCRIPT" \
  --retained-transcript-tsv "$RETAINED_TSV" \
  --transcript-summary-tsv "$SUMMARY_TSV" \
  --transcript-exon-bed "$EXON_BED" \
  --chrom-sizes "$CHROM_SIZES" \
  --doga-ctss-bed "$DOGA_CTSS_REGION_RESTRICTED" \
  --doga-promoter-bed "$DOGA_PROMOTER" \
  --barkbase-atac-bed "$BARKBASE_ATAC" \
  --epicdog-promoter-like-bed "$EPICDOG_PROMOTER" \
  --epicdog-active-bed "$EPICDOG_ACTIVE" \
  --epicdog-super-enhancer-bed "$EPICDOG_SUPER" \
  --outdir "$OUTPUT_ROOT" \
  --run-id "$RUN_ID" &
CHILD_PID=$!
if wait "$CHILD_PID"; then
  CHILD_PID=""
else
  WORKER_CODE=$?
  CHILD_PID=""
  exit "$WORKER_CODE"
fi

printf 'execution_status=success\n' > "$RUN_DIR/execution_completion.marker"
printf 'Step2.2\t%s\tsuccess\t0\t%s\n' "$RUN_ID" "$RUN_DIR" >> "$SUCCESS_FILE"
printf '[%s ] task completed / done execution_status=success total=1 success=1 failed=0 skipped=0\n' "$(/usr/bin/date '+%F %T')"
printf 'please check separately output completeness, distance distribution, denominator and biological interpretation; successful execution does not establish output quality.\n'
