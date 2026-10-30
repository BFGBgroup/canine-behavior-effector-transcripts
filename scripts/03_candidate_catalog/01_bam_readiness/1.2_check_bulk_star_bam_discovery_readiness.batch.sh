#!/bin/bash
set -euo pipefail


SCRIPT="/path/to/repository/scripts/03_candidate_catalog/01_bam_readiness/1.2_check_bulk_star_bam_discovery_readiness.py"
PYTHON_BIN="/usr/bin/python3"
SAMTOOLS="/path/to/software/miniconda3/bin/samtools"

SAMPLE_LIST="/path/to/project/legacy/five-dogs/data/statistics/five_dogs_fastq_pair_stats.tsv"
STAR_ROOT="/path/to/project/legacy/five-dogs/5.cnv/1.star-v5"
CHROM_SIZES="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna.fai"
OUTPUT_ROOT="/path/to/project/legacy/five-dogs/9.new-gene/1.2_check_bulk_star_bam_discovery_readiness-v3"
LOG_DIR="/path/to/repository/scripts/03_candidate_catalog/01_bam_readiness/logs-v3"

WORKERS="19"
MAX_ALIGNMENTS="100000"
MIN_BAM_SIZE_BYTES="10000000"
COMMAND_TIMEOUT_SEC="300"
OVERWRITE=0

if [[ "$#" -ne 0 ]]; then
  echo "[ERROR] This batch does not accept positional arguments." >&2
  echo "Usage: bash $0" >&2
  exit 2
fi

if [[ ! -f "$SCRIPT" ]]; then
  echo "[ERROR] script not found: $SCRIPT" >&2
  exit 1
fi

if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "[ERROR] python not executable: $PYTHON_BIN" >&2
  exit 1
fi

if [[ ! -x "$SAMTOOLS" ]]; then
  echo "[ERROR] samtools not executable: $SAMTOOLS" >&2
  exit 1
fi

if [[ ! -f "$SAMPLE_LIST" ]]; then
  echo "[ERROR] sample list not found: $SAMPLE_LIST" >&2
  exit 1
fi

if [[ ! -d "$STAR_ROOT" ]]; then
  echo "[ERROR] STAR root not found: $STAR_ROOT" >&2
  exit 1
fi

if [[ ! -f "$CHROM_SIZES" ]]; then
  echo "[ERROR] chrom sizes not found: $CHROM_SIZES" >&2
  exit 1
fi

mkdir -p "$LOG_DIR"

MASTER_LOG="${LOG_DIR}/1.2_check_bulk_star_bam_discovery_readiness.batch.master.log"
RUN_LOG="${LOG_DIR}/1.2_check_bulk_star_bam_discovery_readiness.run.log"
SUCCESS_LIST="${LOG_DIR}/1.2_check_bulk_star_bam_discovery_readiness.success.tsv"
FAILED_LIST="${LOG_DIR}/1.2_check_bulk_star_bam_discovery_readiness.failed.tsv"

printf "task\toutdir\trc\trun_log\n" > "$SUCCESS_LIST"
printf "task\toutdir\trc\trun_log\n" > "$FAILED_LIST"

CMD=(
  "$PYTHON_BIN"
  "$SCRIPT"
  --sample-list "$SAMPLE_LIST"
  --star-root "$STAR_ROOT"
  --outdir "$OUTPUT_ROOT"
  --samtools "$SAMTOOLS"
  --chrom-sizes "$CHROM_SIZES"
  --workers "$WORKERS"
  --max-alignments "$MAX_ALIGNMENTS"
  --min-bam-size-bytes "$MIN_BAM_SIZE_BYTES"
  --command-timeout-sec "$COMMAND_TIMEOUT_SEC"
)

if [[ "$OVERWRITE" == "1" ]]; then
  CMD+=(--overwrite)
elif [[ "$OVERWRITE" != "0" ]]; then
  echo "[ERROR] OVERWRITE must be 0 or 1, got: $OVERWRITE" >&2
  exit 2
fi

{
  echo "[batch] $(date '+%F %T') START"
  echo "[batch] SCRIPT=$SCRIPT"
  echo "[batch] PYTHON_BIN=$PYTHON_BIN"
  echo "[batch] SAMTOOLS=$SAMTOOLS"
  echo "[batch] SAMPLE_LIST=$SAMPLE_LIST"
  echo "[batch] STAR_ROOT=$STAR_ROOT"
  echo "[batch] CHROM_SIZES=$CHROM_SIZES"
  echo "[batch] OUTPUT_ROOT=$OUTPUT_ROOT"
  echo "[batch] LOG_DIR=$LOG_DIR"
  echo "[batch] WORKERS=$WORKERS"
  echo "[batch] MAX_ALIGNMENTS=$MAX_ALIGNMENTS"
  echo "[batch] MIN_BAM_SIZE_BYTES=$MIN_BAM_SIZE_BYTES"
  echo "[batch] COMMAND_TIMEOUT_SEC=$COMMAND_TIMEOUT_SEC"
  echo "[batch] OVERWRITE=$OVERWRITE"
  printf "[batch] CMD="
  printf " %q" "${CMD[@]}"
  printf "\n"
  echo "[batch] RUN_LOG=$RUN_LOG"
} | tee -a "$MASTER_LOG"

set +e
"${CMD[@]}" > "$RUN_LOG" 2>&1
RC=$?
set -e

if [[ "$RC" -eq 0 ]]; then
  printf "check_bulk_star_bam_discovery_readiness\t%s\t%s\t%s\n" "$OUTPUT_ROOT" "$RC" "$RUN_LOG" >> "$SUCCESS_LIST"
  echo "[batch] $(date '+%F %T') FINISH success rc=$RC" | tee -a "$MASTER_LOG"
else
  printf "check_bulk_star_bam_discovery_readiness\t%s\t%s\t%s\n" "$OUTPUT_ROOT" "$RC" "$RUN_LOG" >> "$FAILED_LIST"
  echo "[batch] $(date '+%F %T') FINISH failed rc=$RC" | tee -a "$MASTER_LOG"
fi

echo "[batch] RUN_LOG=$RUN_LOG" | tee -a "$MASTER_LOG"
echo "[batch] SUCCESS_LIST=$SUCCESS_LIST" | tee -a "$MASTER_LOG"
echo "[batch] FAILED_LIST=$FAILED_LIST" | tee -a "$MASTER_LOG"

exit "$RC"
