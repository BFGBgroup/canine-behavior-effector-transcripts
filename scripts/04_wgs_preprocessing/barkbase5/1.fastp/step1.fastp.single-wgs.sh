#!/bin/bash
set -euo pipefail


THREADS=50

FASTQ_DIR="/path/to/project/legacy/five-dogs/data/WGS"
OUT_ROOT="/path/to/project/legacy/five-dogs/data/WGS/WGS_fastp"

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <animal_id>" >&2
  exit 1
fi

ANIMAL_ID="$1"

R1="${FASTQ_DIR}/${ANIMAL_ID}_1.fastq.gz"
R2="${FASTQ_DIR}/${ANIMAL_ID}_2.fastq.gz"

SAMPLE_OUT_DIR="${OUT_ROOT}/${ANIMAL_ID}"
mkdir -p "$SAMPLE_OUT_DIR"

O1="${SAMPLE_OUT_DIR}/${ANIMAL_ID}_1.fastq.gz"
O2="${SAMPLE_OUT_DIR}/${ANIMAL_ID}_2.fastq.gz"
HTML="${SAMPLE_OUT_DIR}/${ANIMAL_ID}.fastp.html"
JSON="${SAMPLE_OUT_DIR}/${ANIMAL_ID}.fastp.json"

if ! command -v fastp >/dev/null 2>&1; then
  echo "[ERROR] fastp not found in PATH" >&2
  exit 1
fi

if [[ ! -f "$R1" ]]; then
  echo "[ERROR] R1 not found: $R1" >&2
  exit 1
fi

if [[ ! -f "$R2" ]]; then
  echo "[ERROR] R2 not found: $R2" >&2
  exit 1
fi

if [[ -s "$O1" && -s "$O2" && -s "$HTML" && -s "$JSON" ]]; then
  echo "[skip] ${ANIMAL_ID} already finished."
  echo "[skip] O1=$O1"
  echo "[skip] O2=$O2"
  exit 0
fi

echo "[start] $(date '+%F %T') animal_id=$ANIMAL_ID"
echo "[info] THREADS=$THREADS"
echo "[info] R1=$R1"
echo "[info] R2=$R2"
echo "[info] O1=$O1"
echo "[info] O2=$O2"
echo "[info] HTML=$HTML"
echo "[info] JSON=$JSON"

fastp \
  -i "$R1" \
  -I "$R2" \
  -o "$O1" \
  -O "$O2" \
  -h "$HTML" \
  -j "$JSON" \
  -w "$THREADS" \
  --detect_adapter_for_pe

echo "[done] $(date '+%F %T') animal_id=$ANIMAL_ID"
echo "[done] O1=$O1"
echo "[done] O2=$O2"
