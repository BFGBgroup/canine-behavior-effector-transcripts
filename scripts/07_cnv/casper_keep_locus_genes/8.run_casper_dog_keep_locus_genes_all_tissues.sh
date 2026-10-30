#!/usr/bin/env bash
set -euo pipefail


INPUT_PROJECT_DIR="/path/to/project/legacy/five-dogs/5.cnv"
SCRIPT_ROOT="/path/to/repository/scripts/07_cnv/casper_keep_locus_genes"
SAMPLE_META="${INPUT_PROJECT_DIR}/5.casper_input/bulk_rna.samples.keep_locus_genes.all_genes.for_casper.tsv"
SCRIPT="${SCRIPT_ROOT}/8.run_casper_dog_keep_locus_genes.sh"

CNV_SCALE="${CNV_SCALE:-3}"
LOH_SCALE="${LOH_SCALE:-3}"
BAF_MODE="${BAF_MODE:-hybrid}"
MAX_PARALLEL="${MAX_PARALLEL:-20}"

if [[ ! -s "${SAMPLE_META}" ]]; then
  echo "[ERROR] sample metadata file not found or empty: ${SAMPLE_META}" >&2
  exit 1
fi

if [[ ! -s "${SCRIPT}" ]]; then
  echo "[ERROR] Step 8 script not found or empty: ${SCRIPT}" >&2
  exit 1
fi

mapfile -t tissues < <(awk -F '\t' 'NR > 1 && $3 != "" { print $3 }' "${SAMPLE_META}" | sort -u)

if [[ "${#tissues[@]}" -eq 0 ]]; then
  echo "[ERROR] no tissues parsed from: ${SAMPLE_META}" >&2
  exit 1
fi

echo "[batch] $(date '+%F %T') start"
echo "[batch] sample_meta=${SAMPLE_META}"
echo "[batch] script=${SCRIPT}"
echo "[batch] cnv_scale=${CNV_SCALE}"
echo "[batch] loh_scale=${LOH_SCALE}"
echo "[batch] baf_mode=${BAF_MODE}"
echo "[batch] max_parallel=${MAX_PARALLEL}"
echo "[batch] check_only=${CHECK_ONLY:-0}"
echo "[batch] tissues=${tissues[*]}"

running_jobs=0
pids=()
failed=0

for tissue in "${tissues[@]}"; do
  echo "[batch] $(date '+%F %T') launch tissue=${tissue}"
  bash "${SCRIPT}" "${tissue}" "${CNV_SCALE}" "${LOH_SCALE}" "${BAF_MODE}" &
  pids+=($!)
  running_jobs=$((running_jobs + 1))

  if [[ "${running_jobs}" -ge "${MAX_PARALLEL}" ]]; then
    for pid in "${pids[@]}"; do
      if ! wait "${pid}"; then
        failed=1
      fi
    done
    pids=()
    running_jobs=0
    echo "[batch] $(date '+%F %T') completed one parallel batch"
  fi
done

if [[ "${#pids[@]}" -gt 0 ]]; then
  for pid in "${pids[@]}"; do
    if ! wait "${pid}"; then
      failed=1
    fi
  done
  echo "[batch] $(date '+%F %T') completed final parallel batch"
fi

if [[ "${failed}" -ne 0 ]]; then
  echo "[batch] $(date '+%F %T') finished with failures"
  exit 1
fi

echo "[batch] $(date '+%F %T') all done"
