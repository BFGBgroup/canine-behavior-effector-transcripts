#!/bin/bash
set -euo pipefail


THREADS=5
MAPQ=30

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <animal_id> <term>" >&2
  exit 1
fi

ANIMAL_ID="$1"
TERM="$2"

IN_BASE="/path/to/project/legacy/five-dogs/4.ASE-process/batch-v7/4.3STAR+WASP"
OUT_BASE="/path/to/project/legacy/five-dogs/4.ASE-process/batch-v7/4.4filiter"

IN_DIR="${IN_BASE}/${ANIMAL_ID}/${TERM}/align_wasp"
OUT_DIR="${OUT_BASE}/${ANIMAL_ID}/${TERM}/align_wasp"
mkdir -p "${OUT_DIR}"

IN_BAM="${IN_DIR}/${ANIMAL_ID}_${TERM}.Aligned.sortedByCoord.out.bam"

WASP_PASS_BAM="${OUT_DIR}/${ANIMAL_ID}_${TERM}.WASP_pass.bam"
WASP_PASS_BAI="${WASP_PASS_BAM}.bai"

WASP_PASS_Q30_BAM="${OUT_DIR}/${ANIMAL_ID}_${TERM}.WASP_pass.q30.bam"
WASP_PASS_Q30_BAI="${WASP_PASS_Q30_BAM}.bai"

STATS_TXT="${OUT_DIR}/${ANIMAL_ID}_${TERM}.WASP_MAPQ.filter.stats.txt"
LOG_TXT="${OUT_DIR}/${ANIMAL_ID}_${TERM}.WASP_MAPQ.filter.log"

if ! command -v samtools >/dev/null 2>&1; then
  echo "[ERROR] samtools not found in PATH" >&2
  exit 1
fi

if [[ ! -f "$IN_BAM" ]]; then
  echo "[ERROR] input BAM not found: $IN_BAM" >&2
  exit 1
fi

if ! samtools quickcheck "$IN_BAM"; then
  echo "[ERROR] input BAM failed quickcheck: $IN_BAM" >&2
  exit 1
fi


if [[ -s "$WASP_PASS_Q30_BAM" && -s "$WASP_PASS_Q30_BAI" && -s "$STATS_TXT" ]]; then
  echo "[skip] ${ANIMAL_ID} ${TERM} already finished."
  exit 0
fi

{
  echo "[start] $(date '+%F %T') animal_id=${ANIMAL_ID} term=${TERM}"
  echo "[info] IN_BAM=${IN_BAM}"
  echo "[info] OUT_DIR=${OUT_DIR}"
  echo "[info] MAPQ=${MAPQ}"
} | tee "$LOG_TXT"


echo "[step] filter WASP-pass reads (vW:i:1)" | tee -a "$LOG_TXT"

samtools view -h "$IN_BAM" \
| awk 'BEGIN{OFS="\t"} ($1 ~ /^@/) || ($0 ~ /(^|\t)vW:i:1(\t|$)/)' \
| samtools view -@ "$THREADS" -b -o "$WASP_PASS_BAM" -

samtools index -@ "$THREADS" "$WASP_PASS_BAM"


echo "[step] filter MAPQ >= ${MAPQ}" | tee -a "$LOG_TXT"

samtools view -@ "$THREADS" -b -q "$MAPQ" \
  "$WASP_PASS_BAM" \
  -o "$WASP_PASS_Q30_BAM"

samtools index -@ "$THREADS" "$WASP_PASS_Q30_BAM"


echo "[step] summarize counts" | tee -a "$LOG_TXT"

TOTAL_ALIGNMENTS=$(samtools view -c "$IN_BAM")
WASP_PASS_ALIGNMENTS=$(samtools view -c "$WASP_PASS_BAM")
WASP_PASS_Q30_ALIGNMENTS=$(samtools view -c "$WASP_PASS_Q30_BAM")

{
  echo -e "animal_id\t${ANIMAL_ID}"
  echo -e "term\t${TERM}"
  echo -e "input_bam\t${IN_BAM}"
  echo -e "total_alignments\t${TOTAL_ALIGNMENTS}"
  echo -e "wasp_pass_alignments\t${WASP_PASS_ALIGNMENTS}"
  echo -e "wasp_pass_q30_alignments\t${WASP_PASS_Q30_ALIGNMENTS}"
  echo -e "wasp_pass_bam\t${WASP_PASS_BAM}"
  echo -e "wasp_pass_q30_bam\t${WASP_PASS_Q30_BAM}"
} > "$STATS_TXT"

{
  echo "[done] $(date '+%F %T') animal_id=${ANIMAL_ID} term=${TERM}"
  echo "[done] WASP_PASS_BAM=${WASP_PASS_BAM}"
  echo "[done] WASP_PASS_Q30_BAM=${WASP_PASS_Q30_BAM}"
  echo "[done] STATS_TXT=${STATS_TXT}"
} | tee -a "$LOG_TXT"
