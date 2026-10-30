#!/bin/bash
set -euo pipefail


SCRIPT="/path/to/repository/scripts/03_candidate_catalog/1.5_stringtie3_bulk_transcript_discovery/1.5_stringtie3_bulk_transcript_discovery.py"
PYTHON_BIN="/usr/bin/python3"
STRINGTIE="/path/to/software/legacy_envs/stringtie3/bin/stringtie"

READINESS_TSV="/path/to/project/legacy/five-dogs/9.new-gene/1.2_check_bulk_star_bam_discovery_readiness-v3/bam_discovery_readiness.tsv"
GUIDE_GTF="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.gtf"

OUTPUT_ROOT="/path/to/project/legacy/five-dogs/9.new-gene/1.5_stringtie3_bulk_transcript_discovery-v6"
LOG_DIR="/path/to/repository/scripts/03_candidate_catalog/1.5_stringtie3_bulk_transcript_discovery/logs-v6"

WORKERS="5"
THREADS_PER_SAMPLE="8"
MERGE_THREADS="8"
EXPECTED_SAMPLE_COUNT="19"


STRAND_OPTION="--rf"


NASCENT_OPTION="-N"

RESUME=0
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

if [[ ! -x "$STRINGTIE" ]]; then
  echo "[ERROR] stringtie not executable: $STRINGTIE" >&2
  echo "[HINT] please set batch top STRINGTIE change to actual StringTie3/stringtie executable file absolute path." >&2
  exit 1
fi

for required_file in "$READINESS_TSV" "$GUIDE_GTF"; do
  if [[ ! -f "$required_file" ]]; then
    echo "[ERROR] required input not found: $required_file" >&2
    exit 1
  fi
done

if [[ "$RESUME" == "1" && "$OVERWRITE" == "1" ]]; then
  echo "[ERROR] RESUME and OVERWRITE cannot both be 1." >&2
  exit 2
fi

mkdir -p "$LOG_DIR"

MASTER_LOG="${LOG_DIR}/1.5_stringtie3_bulk_transcript_discovery.batch.master.log"
RUN_LOG="${LOG_DIR}/1.5_stringtie3_bulk_transcript_discovery.run.log"
SUCCESS_LIST="${LOG_DIR}/1.5_stringtie3_bulk_transcript_discovery.success.tsv"
FAILED_LIST="${LOG_DIR}/1.5_stringtie3_bulk_transcript_discovery.failed.tsv"

printf "task\toutdir\trc\trun_log\n" > "$SUCCESS_LIST"
printf "task\toutdir\trc\trun_log\n" > "$FAILED_LIST"

CMD=(
  "$PYTHON_BIN"
  "$SCRIPT"
  --readiness-tsv "$READINESS_TSV"
  --guide-gtf "$GUIDE_GTF"
  --stringtie "$STRINGTIE"
  --outdir "$OUTPUT_ROOT"
  --workers "$WORKERS"
  --threads-per-sample "$THREADS_PER_SAMPLE"
  --merge-threads "$MERGE_THREADS"
  --expected-sample-count "$EXPECTED_SAMPLE_COUNT"
)

if [[ -n "$STRAND_OPTION" ]]; then
  if [[ "$STRAND_OPTION" != "--rf" && "$STRAND_OPTION" != "--fr" ]]; then
    echo "[ERROR] STRAND_OPTION must be empty, --rf, or --fr; got: $STRAND_OPTION" >&2
    exit 2
  fi
  CMD+=("--strand-option=${STRAND_OPTION}")
fi

if [[ -n "$NASCENT_OPTION" ]]; then
  if [[ "$NASCENT_OPTION" != "-N" && "$NASCENT_OPTION" != "--nasc" ]]; then
    echo "[ERROR] NASCENT_OPTION must be empty, -N, or --nasc; got: $NASCENT_OPTION" >&2
    exit 2
  fi
  CMD+=("--nascent-option=${NASCENT_OPTION}")
fi

if [[ "$RESUME" == "1" ]]; then
  CMD+=(--resume)
elif [[ "$RESUME" != "0" ]]; then
  echo "[ERROR] RESUME must be 0 or 1, got: $RESUME" >&2
  exit 2
fi

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
  echo "[batch] STRINGTIE=$STRINGTIE"
  echo "[batch] READINESS_TSV=$READINESS_TSV"
  echo "[batch] GUIDE_GTF=$GUIDE_GTF"
  echo "[batch] OUTPUT_ROOT=$OUTPUT_ROOT"
  echo "[batch] LOG_DIR=$LOG_DIR"
  echo "[batch] WORKERS=$WORKERS"
  echo "[batch] THREADS_PER_SAMPLE=$THREADS_PER_SAMPLE"
  echo "[batch] MERGE_THREADS=$MERGE_THREADS"
  echo "[batch] EXPECTED_SAMPLE_COUNT=$EXPECTED_SAMPLE_COUNT"
  echo "[batch] STRAND_OPTION=${STRAND_OPTION:-NA}"
  echo "[batch] NASCENT_OPTION=${NASCENT_OPTION:-NA}"
  echo "[batch] RESUME=$RESUME"
  echo "[batch] OVERWRITE=$OVERWRITE"
  printf "[batch] CMD="
  printf " %q" "${CMD[@]}"
  printf "\n"
  echo "[batch] RUN_LOG=$RUN_LOG"
} > "$MASTER_LOG"

set +e
"${CMD[@]}" > "$RUN_LOG" 2>&1
RC=$?
set -e

if [[ "$RC" -eq 0 ]]; then
  printf "stringtie3_bulk_transcript_discovery\t%s\t%s\t%s\n" "$OUTPUT_ROOT" "$RC" "$RUN_LOG" >> "$SUCCESS_LIST"
  echo "[batch] $(date '+%F %T') FINISH success rc=$RC" | tee -a "$MASTER_LOG"
else
  printf "stringtie3_bulk_transcript_discovery\t%s\t%s\t%s\n" "$OUTPUT_ROOT" "$RC" "$RUN_LOG" >> "$FAILED_LIST"
  echo "[batch] $(date '+%F %T') FINISH failed rc=$RC" | tee -a "$MASTER_LOG"
fi

echo "[batch] RUN_LOG=$RUN_LOG" | tee -a "$MASTER_LOG"
echo "[batch] SUCCESS_LIST=$SUCCESS_LIST" | tee -a "$MASTER_LOG"
echo "[batch] FAILED_LIST=$FAILED_LIST" | tee -a "$MASTER_LOG"

exit "$RC"
