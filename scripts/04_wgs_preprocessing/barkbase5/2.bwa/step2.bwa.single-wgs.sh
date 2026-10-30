#!/bin/bash
set -euo pipefail


THREADS=12

REF="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna"

FASTQ_ROOT="/path/to/project/legacy/five-dogs/data/WGS/WGS_fastp"
OUT_ROOT="/path/to/project/legacy/five-dogs/data/WGS/WGS_bwa"

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <animal_id>" >&2
  exit 1
fi

ANIMAL_ID="$1"

R1="${FASTQ_ROOT}/${ANIMAL_ID}/${ANIMAL_ID}_1.fastq.gz"
R2="${FASTQ_ROOT}/${ANIMAL_ID}/${ANIMAL_ID}_2.fastq.gz"

SAMPLE_OUT_DIR="${OUT_ROOT}/${ANIMAL_ID}"
mkdir -p "$SAMPLE_OUT_DIR"

BAM="${SAMPLE_OUT_DIR}/${ANIMAL_ID}.sorted.bam"
BAI="${BAM}.bai"
FLAGSTAT="${SAMPLE_OUT_DIR}/${ANIMAL_ID}.flagstat.txt"
IDXSTATS="${SAMPLE_OUT_DIR}/${ANIMAL_ID}.idxstats.txt"

RG="@RG\tID:${ANIMAL_ID}\tSM:${ANIMAL_ID}\tPL:ILLUMINA\tLB:${ANIMAL_ID}_lib1\tPU:${ANIMAL_ID}_unit1"

if ! command -v bwa >/dev/null 2>&1; then
  echo "[ERROR] bwa not found in PATH" >&2
  exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
  echo "[ERROR] samtools not found in PATH" >&2
  exit 1
fi

if [[ ! -f "$REF" ]]; then
  echo "[ERROR] REF not found: $REF" >&2
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

if [[ ! -f "${REF}.bwt" ]]; then
  echo "[ERROR] BWA index not found: ${REF}.bwt" >&2
  echo "[ERROR] Please run: bwa index $REF" >&2
  exit 1
fi

if [[ ! -f "${REF}.fai" ]]; then
  echo "[ERROR] FASTA index not found: ${REF}.fai" >&2
  echo "[ERROR] Please run: samtools faidx $REF" >&2
  exit 1
fi

if [[ -s "$BAM" && -s "$BAI" && -s "$FLAGSTAT" && -s "$IDXSTATS" ]]; then
  echo "[skip] ${ANIMAL_ID} already finished."
  echo "[skip] BAM=$BAM"
  exit 0
fi

echo "[start] $(date '+%F %T') animal_id=$ANIMAL_ID"
echo "[info] THREADS=$THREADS"
echo "[info] REF=$REF"
echo "[info] R1=$R1"
echo "[info] R2=$R2"
echo "[info] BAM=$BAM"

bwa mem -t "$THREADS" -R "$RG" "$REF" "$R1" "$R2" \
  | samtools sort -@ "$THREADS" -o "$BAM" -

samtools index -@ "$THREADS" "$BAM"
samtools flagstat "$BAM" > "$FLAGSTAT"
samtools idxstats "$BAM" > "$IDXSTATS"

echo "[done] $(date '+%F %T') animal_id=$ANIMAL_ID"
echo "[done] BAM=$BAM"
echo "[done] BAI=$BAI"
echo "[done] FLAGSTAT=$FLAGSTAT"
echo "[done] IDXSTATS=$IDXSTATS"
