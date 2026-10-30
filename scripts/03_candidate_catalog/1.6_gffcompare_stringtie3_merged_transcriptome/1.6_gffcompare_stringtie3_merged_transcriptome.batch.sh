#!/bin/bash
set -euo pipefail


SUMMARY_SCRIPT="/path/to/repository/scripts/03_candidate_catalog/1.6_gffcompare_stringtie3_merged_transcriptome/1.6_summarize_gffcompare_stringtie3.py"
EXON_BED_SCRIPT="/path/to/repository/scripts/03_candidate_catalog/1.6_gffcompare_stringtie3_merged_transcriptome/1.6_extract_gffcompare_transcript_exon_bed.py"
PYTHON_BIN="/usr/bin/python3"
GFFCOMPARE="/path/to/software/legacy_envs/stringtie3/bin/gffcompare"

MERGED_GTF="/path/to/project/legacy/five-dogs/9.new-gene/1.5_stringtie3_bulk_transcript_discovery-v6/merged.stringtie3.gtf"
REFERENCE_GTF="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.gtf"
STEP15_QC="/path/to/project/legacy/five-dogs/9.new-gene/1.5_stringtie3_bulk_transcript_discovery-v6/stringtie3_discovery_qc_summary.tsv"

OUTPUT_ROOT="/path/to/project/legacy/five-dogs/9.new-gene/1.6_gffcompare_stringtie3_merged_transcriptome-v6"
LOG_DIR="/path/to/repository/scripts/03_candidate_catalog/1.6_gffcompare_stringtie3_merged_transcriptome/logs-v6"
OUT_PREFIX="${OUTPUT_ROOT}/gffcmp.stringtie3"
QUERY_GTF="${OUTPUT_ROOT}/merged.stringtie3.gtf"
OVERWRITE=0

if [[ "$#" -ne 0 ]]; then
  echo "[ERROR] This batch does not accept positional arguments." >&2
  echo "Usage: bash $0" >&2
  exit 2
fi

for required_script in "$SUMMARY_SCRIPT" "$EXON_BED_SCRIPT"; do
  if [[ ! -f "$required_script" ]]; then
    echo "[ERROR] required script not found: $required_script" >&2
    exit 1
  fi
done

if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "[ERROR] python not executable: $PYTHON_BIN" >&2
  exit 1
fi

if [[ ! -x "$GFFCOMPARE" ]]; then
  echo "[ERROR] gffcompare not executable: $GFFCOMPARE" >&2
  echo "[HINT] please set batch top GFFCOMPARE change to actual gffcompare executable file absolute path." >&2
  exit 1
fi

for required_file in "$MERGED_GTF" "$REFERENCE_GTF" "$STEP15_QC"; do
  if [[ ! -s "$required_file" ]]; then
    echo "[ERROR] required input missing or empty: $required_file" >&2
    exit 1
  fi
done

if [[ "$OVERWRITE" != "0" && "$OVERWRITE" != "1" ]]; then
  echo "[ERROR] OVERWRITE must be 0 or 1, got: $OVERWRITE" >&2
  exit 2
fi

if [[ -d "$OUTPUT_ROOT" && "$OVERWRITE" == "0" ]]; then
  shopt -s nullglob
  existing_outputs=("${OUT_PREFIX}".* "${OUTPUT_ROOT}/1.6_gffcompare_"*.tsv "${OUTPUT_ROOT}/1.6_gffcompare_query."* "${OUTPUT_ROOT}/transcript_exon"* )
  shopt -u nullglob
  for protected_path in "$QUERY_GTF" "$OUT_PREFIX"; do
    if [[ -e "$protected_path" ]]; then
      existing_outputs+=("$protected_path")
    fi
  done
  if [[ "${#existing_outputs[@]}" -gt 0 ]]; then
    echo "[ERROR] output files already exist under: $OUTPUT_ROOT" >&2
    echo "[HINT] such as requires overwrite the same step output, please first confirm after set batch top OVERWRITE change to 1." >&2
    printf '[ERROR] existing: %s\n' "${existing_outputs[@]}" >&2
    exit 1
  fi
fi

mkdir -p "$OUTPUT_ROOT" "$LOG_DIR"


cp -f "$MERGED_GTF" "$QUERY_GTF"

MASTER_LOG="${LOG_DIR}/1.6_gffcompare_stringtie3.batch.master.log"
RUN_LOG="${LOG_DIR}/1.6_gffcompare_stringtie3.run.log"
SUMMARY_LOG="${LOG_DIR}/1.6_gffcompare_stringtie3.summary.log"
EXON_BED_LOG="${LOG_DIR}/1.6_gffcompare_transcript_exon_bed.log"
SUCCESS_LIST="${LOG_DIR}/1.6_gffcompare_stringtie3.success.tsv"
FAILED_LIST="${LOG_DIR}/1.6_gffcompare_stringtie3.failed.tsv"

printf "task\toutdir\trc\trun_log\n" > "$SUCCESS_LIST"
printf "task\toutdir\trc\trun_log\n" > "$FAILED_LIST"

CMD=(
  "$GFFCOMPARE"
  -r "$REFERENCE_GTF"
  --chr-stats
  -o "$OUT_PREFIX"
  "$QUERY_GTF"
)

{
  echo "[batch] $(date '+%F %T') START"
  echo "[batch] SUMMARY_SCRIPT=$SUMMARY_SCRIPT"
  echo "[batch] EXON_BED_SCRIPT=$EXON_BED_SCRIPT"
  echo "[batch] PYTHON_BIN=$PYTHON_BIN"
  echo "[batch] GFFCOMPARE=$GFFCOMPARE"
  echo "[batch] MERGED_GTF=$MERGED_GTF"
  echo "[batch] QUERY_GTF=$QUERY_GTF"
  echo "[batch] REFERENCE_GTF=$REFERENCE_GTF"
  echo "[batch] STEP15_QC=$STEP15_QC"
  echo "[batch] OUTPUT_ROOT=$OUTPUT_ROOT"
  echo "[batch] LOG_DIR=$LOG_DIR"
  echo "[batch] OUT_PREFIX=$OUT_PREFIX"
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

if [[ "$RC" -ne 0 ]]; then
  printf "gffcompare_stringtie3\t%s\t%s\t%s\n" "$OUTPUT_ROOT" "$RC" "$RUN_LOG" >> "$FAILED_LIST"
  echo "[batch] $(date '+%F %T') FINISH failed rc=$RC" | tee -a "$MASTER_LOG"
  echo "[batch] RUN_LOG=$RUN_LOG" | tee -a "$MASTER_LOG"
  exit "$RC"
fi

shopt -s nullglob
TMAP_CANDIDATES=()
for candidate in \
  "${OUT_PREFIX}".*.tmap \
  "${OUT_PREFIX}".tmap; do
  if [[ -s "$candidate" ]]; then
    TMAP_CANDIDATES+=("$candidate")
  fi
done
REFMAP_CANDIDATES=()
for candidate in \
  "${OUT_PREFIX}".*.refmap \
  "${OUT_PREFIX}".refmap; do
  if [[ -s "$candidate" ]]; then
    REFMAP_CANDIDATES+=("$candidate")
  fi
done
shopt -u nullglob

if [[ "${#TMAP_CANDIDATES[@]}" -eq 1 ]]; then
  TMAP_PATH="${TMAP_CANDIDATES[0]}"
elif [[ "${#TMAP_CANDIDATES[@]}" -eq 0 ]]; then
  TMAP_PATH="NA"
else
  echo "[ERROR] expected zero or one gffcompare .tmap, found ${#TMAP_CANDIDATES[@]}" >&2
  printf '[ERROR] tmap candidate: %s\n' "${TMAP_CANDIDATES[@]}"
  printf "gffcompare_stringtie3\t%s\t%s\t%s\n" "$OUTPUT_ROOT" "3" "$RUN_LOG" >> "$FAILED_LIST"
  exit 3
fi

if [[ "${#REFMAP_CANDIDATES[@]}" -eq 1 ]]; then
  REFMAP_PATH="${REFMAP_CANDIDATES[0]}"
else
  REFMAP_PATH="NA"
fi

if [[ -s "${OUT_PREFIX}.stats" ]]; then
  STATS_PATH="${OUT_PREFIX}.stats"
elif [[ -s "${OUT_PREFIX}" ]]; then
  STATS_PATH="${OUT_PREFIX}"
else
  echo "[ERROR] expected gffcompare stats output but found neither ${OUT_PREFIX}.stats nor ${OUT_PREFIX}" >&2
  printf "gffcompare_stringtie3\t%s\t%s\t%s\n" "$OUTPUT_ROOT" "5" "$RUN_LOG" >> "$FAILED_LIST"
  exit 5
fi

if [[ -s "${OUT_PREFIX}.annotated.gtf" ]]; then
  ANNOTATED_GTF="${OUT_PREFIX}.annotated.gtf"
elif [[ -s "${OUT_PREFIX}.combined.gtf" ]]; then
  ANNOTATED_GTF="${OUT_PREFIX}.combined.gtf"
else
  echo "[ERROR] expected gffcompare annotated/combined GTF but found neither ${OUT_PREFIX}.annotated.gtf nor ${OUT_PREFIX}.combined.gtf" >&2
  printf "gffcompare_stringtie3\t%s\t%s\t%s\n" "$OUTPUT_ROOT" "4" "$RUN_LOG" >> "$FAILED_LIST"
  exit 4
fi

if [[ "$TMAP_PATH" != "NA" && -e "$TMAP_PATH" ]]; then
  ln -sfn "$TMAP_PATH" "${OUTPUT_ROOT}/1.6_gffcompare_query.tmap"
fi

if [[ "$REFMAP_PATH" != "NA" && -e "$REFMAP_PATH" ]]; then
  ln -sfn "$REFMAP_PATH" "${OUTPUT_ROOT}/1.6_gffcompare_query.refmap"
fi

SUMMARY_CMD=(
  "$PYTHON_BIN"
  "$SUMMARY_SCRIPT"
  --tmap "$TMAP_PATH"
  --annotated-gtf "$ANNOTATED_GTF"
  --tracking "${OUT_PREFIX}.tracking"
  --stats "$STATS_PATH"
  --loci "${OUT_PREFIX}.loci"
  --refmap "$REFMAP_PATH"
  --merged-gtf "$QUERY_GTF"
  --reference-gtf "$REFERENCE_GTF"
  --step15-qc "$STEP15_QC"
  --outdir "$OUTPUT_ROOT"
  --gffcompare-prefix "$OUT_PREFIX"
  --gffcompare-command "${CMD[*]}"
)

{
  echo "[batch] $(date '+%F %T') gffcompare success rc=$RC"
  echo "[batch] TMAP_PATH=$TMAP_PATH"
  echo "[batch] REFMAP_PATH=$REFMAP_PATH"
  echo "[batch] STATS_PATH=$STATS_PATH"
  echo "[batch] ANNOTATED_GTF=$ANNOTATED_GTF"
  printf "[batch] SUMMARY_CMD="
  printf " %q" "${SUMMARY_CMD[@]}"
  printf "\n"
  echo "[batch] SUMMARY_LOG=$SUMMARY_LOG"
} | tee -a "$MASTER_LOG"

set +e
"${SUMMARY_CMD[@]}" > "$SUMMARY_LOG" 2>&1
SUMMARY_RC=$?
set -e

if [[ "$SUMMARY_RC" -ne 0 ]]; then
  printf "gffcompare_stringtie3_summary\t%s\t%s\t%s\n" "$OUTPUT_ROOT" "$SUMMARY_RC" "$SUMMARY_LOG" >> "$FAILED_LIST"
  echo "[batch] $(date '+%F %T') FINISH failed rc=$RC summary_rc=$SUMMARY_RC" | tee -a "$MASTER_LOG"
  echo "[batch] SUMMARY_LOG=$SUMMARY_LOG" | tee -a "$MASTER_LOG"
  exit "$SUMMARY_RC"
fi

CLASSIFICATION_TSV="${OUTPUT_ROOT}/1.6_gffcompare_transcript_classification.tsv"
EXON_BED="${OUTPUT_ROOT}/transcript_exon.bed"
EXON_BED_QC="${OUTPUT_ROOT}/transcript_exon_bed_qc_summary.tsv"
EXON_BED_MANIFEST="${OUTPUT_ROOT}/transcript_exon_bed_manifest.tsv"
EXON_BED_COLUMNS="${OUTPUT_ROOT}/transcript_exon_bed.columns.tsv"

EXON_BED_CMD=(
  "$PYTHON_BIN"
  "$EXON_BED_SCRIPT"
  --annotated-gtf "$ANNOTATED_GTF"
  --classification-tsv "$CLASSIFICATION_TSV"
  --out-bed "$EXON_BED"
  --qc-summary "$EXON_BED_QC"
  --manifest "$EXON_BED_MANIFEST"
  --columns-tsv "$EXON_BED_COLUMNS"
)

{
  echo "[batch] $(date '+%F %T') summary success summary_rc=$SUMMARY_RC"
  echo "[batch] CLASSIFICATION_TSV=$CLASSIFICATION_TSV"
  echo "[batch] EXON_BED=$EXON_BED"
  printf "[batch] EXON_BED_CMD="
  printf " %q" "${EXON_BED_CMD[@]}"
  printf "\n"
  echo "[batch] EXON_BED_LOG=$EXON_BED_LOG"
} | tee -a "$MASTER_LOG"

set +e
"${EXON_BED_CMD[@]}" > "$EXON_BED_LOG" 2>&1
EXON_BED_RC=$?
set -e

if [[ "$EXON_BED_RC" -eq 0 ]]; then
  printf "gffcompare_stringtie3\t%s\t%s\t%s\n" "$OUTPUT_ROOT" "$RC" "$RUN_LOG" >> "$SUCCESS_LIST"
  echo "[batch] $(date '+%F %T') FINISH success rc=$RC summary_rc=$SUMMARY_RC exon_bed_rc=$EXON_BED_RC" | tee -a "$MASTER_LOG"
else
  printf "gffcompare_transcript_exon_bed\t%s\t%s\t%s\n" "$OUTPUT_ROOT" "$EXON_BED_RC" "$EXON_BED_LOG" >> "$FAILED_LIST"
  echo "[batch] $(date '+%F %T') FINISH failed rc=$RC summary_rc=$SUMMARY_RC exon_bed_rc=$EXON_BED_RC" | tee -a "$MASTER_LOG"
fi

echo "[batch] RUN_LOG=$RUN_LOG" | tee -a "$MASTER_LOG"
echo "[batch] SUMMARY_LOG=$SUMMARY_LOG" | tee -a "$MASTER_LOG"
echo "[batch] EXON_BED_LOG=$EXON_BED_LOG" | tee -a "$MASTER_LOG"
echo "[batch] SUCCESS_LIST=$SUCCESS_LIST" | tee -a "$MASTER_LOG"
echo "[batch] FAILED_LIST=$FAILED_LIST" | tee -a "$MASTER_LOG"

exit "$EXON_BED_RC"
