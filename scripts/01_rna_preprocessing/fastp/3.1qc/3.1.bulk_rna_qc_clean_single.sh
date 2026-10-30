#!/usr/bin/env bash
set -euo pipefail


if [[ $# -ne 2 ]]; then
  echo "Usage: bash $0 <animal_id> <term>" >&2
  exit 1
fi

ANIMAL="$1"
TERM="$2"

RAW_DIR="/path/to/project/legacy/five-dogs/data/raw-rna-data/${TERM}"
OUT_ROOT="/path/to/project/legacy/five-dogs/3.bulk-rna-qc/3.1qc/${ANIMAL}/${TERM}"

R1="${RAW_DIR}/${ANIMAL}_${TERM}_1.fastq.gz"
R2="${RAW_DIR}/${ANIMAL}_${TERM}_2.fastq.gz"

RAW_QC_DIR="${OUT_ROOT}/00.raw_fastqc"
CLEAN_DIR="${OUT_ROOT}/01.clean_fastp"
CLEAN_QC_DIR="${OUT_ROOT}/02.clean_fastqc"
MQC_DIR="${OUT_ROOT}/03.multiqc"
LOG_DIR="${OUT_ROOT}/logs"

FASTQC_THREADS=10
FASTP_THREADS=10

mkdir -p "$RAW_QC_DIR" "$CLEAN_DIR" "$CLEAN_QC_DIR" "$MQC_DIR" "$LOG_DIR"

for cmd in fastqc fastp multiqc gzip; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "[ERROR] command not found: $cmd" >&2
    exit 127
  fi
done

if [[ ! -f "$R1" ]]; then
  echo "[ERROR] R1 not found: $R1" >&2
  exit 1
fi

if [[ ! -f "$R2" ]]; then
  echo "[ERROR] R2 not found: $R2" >&2
  exit 1
fi

echo "[INFO] ANIMAL=$ANIMAL"
echo "[INFO] TERM=$TERM"
echo "[INFO] R1=$R1"
echo "[INFO] R2=$R2"
echo "[INFO] OUT_ROOT=$OUT_ROOT"

echo "[STEP] gzip check"
gzip -t "$R1"
gzip -t "$R2"

echo "[STEP] raw FastQC"
fastqc \
  -t "$FASTQC_THREADS" \
  -o "$RAW_QC_DIR" \
  "$R1" "$R2" \
  > "${LOG_DIR}/${ANIMAL}_${TERM}.raw_fastqc.log" 2>&1

CLEAN_R1="${CLEAN_DIR}/${ANIMAL}_${TERM}_1.clean.fastq.gz"
CLEAN_R2="${CLEAN_DIR}/${ANIMAL}_${TERM}_2.clean.fastq.gz"
FASTP_HTML="${CLEAN_DIR}/${ANIMAL}_${TERM}.fastp.html"
FASTP_JSON="${CLEAN_DIR}/${ANIMAL}_${TERM}.fastp.json"

echo "[STEP] fastp"
fastp \
  -w "$FASTP_THREADS" \
  -i "$R1" \
  -I "$R2" \
  -o "$CLEAN_R1" \
  -O "$CLEAN_R2" \
  --detect_adapter_for_pe \
  -q 20 \
  -u 30 \
  -n 5 \
  -l 30 \
  -h "$FASTP_HTML" \
  -j "$FASTP_JSON" \
  > "${LOG_DIR}/${ANIMAL}_${TERM}.fastp.log" 2>&1

echo "[STEP] clean FastQC"
fastqc \
  -t "$FASTQC_THREADS" \
  -o "$CLEAN_QC_DIR" \
  "$CLEAN_R1" "$CLEAN_R2" \
  > "${LOG_DIR}/${ANIMAL}_${TERM}.clean_fastqc.log" 2>&1

mkdir -p "${MQC_DIR}/raw_fastqc" "${MQC_DIR}/fastp" "${MQC_DIR}/clean_fastqc"

echo "[STEP] MultiQC"
multiqc "$RAW_QC_DIR" -o "${MQC_DIR}/raw_fastqc" \
  > "${LOG_DIR}/${ANIMAL}_${TERM}.multiqc_raw.log" 2>&1

multiqc "$CLEAN_DIR" -o "${MQC_DIR}/fastp" \
  > "${LOG_DIR}/${ANIMAL}_${TERM}.multiqc_fastp.log" 2>&1

multiqc "$CLEAN_QC_DIR" -o "${MQC_DIR}/clean_fastqc" \
  > "${LOG_DIR}/${ANIMAL}_${TERM}.multiqc_clean.log" 2>&1

echo "[DONE] ${ANIMAL} ${TERM}"
