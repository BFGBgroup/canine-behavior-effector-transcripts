#!/bin/bash
set -euo pipefail


MIN_BQ=20
MIN_MQ=30

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <animal_id> <term>" >&2
  exit 1
fi

ANIMAL_ID="$1"
TERM="$2"

REF="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna"

VCF_BASE="/path/to/project/legacy/five-dogs/4.ASE-process/batch-v7/4.5filiter-snp-vcf"
BAM_BASE="/path/to/project/legacy/five-dogs/4.ASE-process/batch-v7/4.4filiter"
OUT_BASE="/path/to/project/legacy/five-dogs/4.ASE-process/batch-v7/4.6bulk-ase-count"

VCF="${VCF_BASE}/${ANIMAL_ID}/${ANIMAL_ID}.minimal.het.snp.unique_pos_only.vcf.gz"
BAM="${BAM_BASE}/${ANIMAL_ID}/${TERM}/align_wasp/${ANIMAL_ID}_${TERM}.WASP_pass.q30.bam"

OUT_DIR="${OUT_BASE}/${ANIMAL_ID}/${TERM}"
mkdir -p "$OUT_DIR"

OUT_TABLE="${OUT_DIR}/${ANIMAL_ID}_${TERM}.ASEReadCounter.table"
OUT_LOG="${OUT_DIR}/${ANIMAL_ID}_${TERM}.ASEReadCounter.log"

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

if [[ ! -f "$VCF" ]]; then
  echo "[ERROR] VCF not found: $VCF" >&2
  exit 1
fi

if [[ ! -f "${VCF}.tbi" ]]; then
  echo "[ERROR] VCF index not found: ${VCF}.tbi" >&2
  exit 1
fi

if [[ ! -f "$BAM" ]]; then
  echo "[ERROR] BAM not found: $BAM" >&2
  exit 1
fi

if [[ ! -f "${BAM}.bai" ]]; then
  echo "[ERROR] BAM index not found: ${BAM}.bai" >&2
  exit 1
fi

if [[ -s "$OUT_TABLE" ]]; then
  echo "[skip] ${ANIMAL_ID} ${TERM} already finished."
  exit 0
fi

if ! samtools quickcheck "$BAM"; then
  echo "[ERROR] BAM failed quickcheck: $BAM" >&2
  exit 1
fi

RG_COUNT="$(samtools view -H "$BAM" | awk '/^@RG/{c++} END{print c+0}')"
if [[ "$RG_COUNT" -eq 0 ]]; then
  echo "[ERROR] BAM has no @RG header lines, but GATK ASEReadCounter requires read groups by default: $BAM" >&2
  exit 1
fi

{
  echo "[start] $(date '+%F %T') animal_id=${ANIMAL_ID} term=${TERM}"
  echo "[info] REF=${REF}"
  echo "[info] VCF=${VCF}"
  echo "[info] BAM=${BAM}"
  echo "[info] BAM_RG_COUNT=${RG_COUNT}"
  echo "[info] MIN_BQ=${MIN_BQ}"
  echo "[info] MIN_MQ=${MIN_MQ}"
  echo "[info] OUT_TABLE=${OUT_TABLE}"
} | tee "$OUT_LOG"

gatk --java-options "-Xmx32g" ASEReadCounter \
  -R "$REF" \
  -I "$BAM" \
  -V "$VCF" \
  -O "$OUT_TABLE" \
  --min-base-quality "$MIN_BQ" \
  --min-mapping-quality "$MIN_MQ" \
  > >(tee -a "$OUT_LOG") 2> >(tee -a "$OUT_LOG" >&2)

echo "[done] $(date '+%F %T') animal_id=${ANIMAL_ID} term=${TERM}" | tee -a "$OUT_LOG"
echo "[done] OUT_TABLE=${OUT_TABLE}" | tee -a "$OUT_LOG"
