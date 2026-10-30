#!/bin/bash
set -euo pipefail


THREADS=100
JAVA_MEM="60g"

REF="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna"
DICT="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.dict"

BAM_ROOT="/path/to/project/legacy/five-dogs/data/WGS/WGS_bwa"
OUT_ROOT="/path/to/project/legacy/five-dogs/data/WGS/WGS_gatk"

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <animal_id>" >&2
  exit 1
fi

ANIMAL_ID="$1"

IN_BAM="${BAM_ROOT}/${ANIMAL_ID}/${ANIMAL_ID}.sorted.bam"
IN_BAI="${IN_BAM}.bai"

SAMPLE_OUT_DIR="${OUT_ROOT}/${ANIMAL_ID}"
TMP_DIR="${SAMPLE_OUT_DIR}/tmp"
mkdir -p "$SAMPLE_OUT_DIR" "$TMP_DIR"

MARKDUP_BAM="${SAMPLE_OUT_DIR}/${ANIMAL_ID}.sorted.markdup.bam"
MARKDUP_BAI="${MARKDUP_BAM}.bai"
METRICS="${SAMPLE_OUT_DIR}/${ANIMAL_ID}.markdup.metrics.txt"

RAW_VCF="${SAMPLE_OUT_DIR}/${ANIMAL_ID}.raw.vcf.gz"
RAW_VCF_TBI="${RAW_VCF}.tbi"

RAW_SNP_VCF="${SAMPLE_OUT_DIR}/${ANIMAL_ID}.raw.snp.vcf.gz"
RAW_SNP_VCF_TBI="${RAW_SNP_VCF}.tbi"

FILTERED_SNP_VCF="${SAMPLE_OUT_DIR}/${ANIMAL_ID}.filtered.snp.vcf.gz"
FILTERED_SNP_VCF_TBI="${FILTERED_SNP_VCF}.tbi"

if ! command -v gatk >/dev/null 2>&1; then
  echo "[ERROR] gatk not found in PATH" >&2
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

if [[ ! -f "${REF}.fai" ]]; then
  echo "[ERROR] FASTA index not found: ${REF}.fai" >&2
  exit 1
fi

if [[ ! -f "$DICT" ]]; then
  echo "[ERROR] sequence dictionary not found: $DICT" >&2
  echo "[ERROR] Please run: gatk CreateSequenceDictionary -R $REF" >&2
  exit 1
fi

if [[ ! -f "$IN_BAM" ]]; then
  echo "[ERROR] input BAM not found: $IN_BAM" >&2
  exit 1
fi

if [[ ! -f "$IN_BAI" ]]; then
  echo "[ERROR] input BAI not found: $IN_BAI" >&2
  exit 1
fi

if [[ -s "$MARKDUP_BAM" && -s "$MARKDUP_BAI" && -s "$METRICS" \
   && -s "$RAW_VCF" && -s "$RAW_VCF_TBI" \
   && -s "$RAW_SNP_VCF" && -s "$RAW_SNP_VCF_TBI" \
   && -s "$FILTERED_SNP_VCF" && -s "$FILTERED_SNP_VCF_TBI" ]]; then
  echo "[skip] ${ANIMAL_ID} already finished."
  exit 0
fi

echo "[start] $(date '+%F %T') animal_id=$ANIMAL_ID"
echo "[info] THREADS=$THREADS"
echo "[info] JAVA_MEM=$JAVA_MEM"
echo "[info] REF=$REF"
echo "[info] IN_BAM=$IN_BAM"
echo "[info] OUTDIR=$SAMPLE_OUT_DIR"

echo "[step] MarkDuplicates"
gatk --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${TMP_DIR}" MarkDuplicates \
  -I "$IN_BAM" \
  -O "$MARKDUP_BAM" \
  -M "$METRICS" \
  --CREATE_INDEX true \
  --TMP_DIR "$TMP_DIR"

echo "[step] HaplotypeCaller"
gatk --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${TMP_DIR}" HaplotypeCaller \
  -R "$REF" \
  -I "$MARKDUP_BAM" \
  -O "$RAW_VCF" \
  --native-pair-hmm-threads "$THREADS"

echo "[step] SelectVariants SNP"
gatk --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${TMP_DIR}" SelectVariants \
  -R "$REF" \
  -V "$RAW_VCF" \
  --select-type-to-include SNP \
  -O "$RAW_SNP_VCF"

echo "[step] VariantFiltration SNP"
gatk --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${TMP_DIR}" VariantFiltration \
  -R "$REF" \
  -V "$RAW_SNP_VCF" \
  -O "$FILTERED_SNP_VCF" \
  --filter-name "QD_lt_2" --filter-expression "QD < 2.0" \
  --filter-name "FS_gt_60" --filter-expression "FS > 60.0" \
  --filter-name "MQ_lt_40" --filter-expression "MQ < 40.0" \
  --filter-name "SOR_gt_3" --filter-expression "SOR > 3.0" \
  --filter-name "MQRankSum_lt_-12.5" --filter-expression "MQRankSum < -12.5" \
  --filter-name "ReadPosRankSum_lt_-8" --filter-expression "ReadPosRankSum < -8.0"

echo "[done] $(date '+%F %T') animal_id=$ANIMAL_ID"
echo "[done] MARKDUP_BAM=$MARKDUP_BAM"
echo "[done] METRICS=$METRICS"
echo "[done] RAW_VCF=$RAW_VCF"
echo "[done] RAW_SNP_VCF=$RAW_SNP_VCF"
echo "[done] FILTERED_SNP_VCF=$FILTERED_SNP_VCF"
