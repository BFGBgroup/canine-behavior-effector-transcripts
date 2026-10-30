#!/bin/bash
set -euo pipefail


THREADS=10

RAW_BASE="/path/to/project/legacy/five-dogs/3.bulk-rna-qc/3.1qc"
STAR_INDEX_BASE="/path/to/project/legacy/five-dogs/4.ASE-process/batch-v7/4.1star-index"
VCF_BASE="/path/to/project/legacy/five-dogs/4.ASE-process/batch-v7/4.2WGS_ti_snp"
OUT_BASE="/path/to/project/legacy/five-dogs/4.ASE-process/batch-v7/4.3STAR+WASP"

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <animal_id> <term>" >&2
  exit 1
fi

ANIMAL_ID="$1"
TERM="$2"

R1="${RAW_BASE}/${ANIMAL_ID}/${TERM}/01.clean_fastp/${ANIMAL_ID}_${TERM}_1.clean.fastq.gz"
R2="${RAW_BASE}/${ANIMAL_ID}/${TERM}/01.clean_fastp/${ANIMAL_ID}_${TERM}_2.clean.fastq.gz"

STAR_INDEX_DIR="${STAR_INDEX_BASE}/${ANIMAL_ID}/${TERM}/STARindex"


VAR_VCF="${VCF_BASE}/${ANIMAL_ID}/${ANIMAL_ID}.minimal.het.snp.vcf"

SAMPLE_OUT_DIR="${OUT_BASE}/${ANIMAL_ID}/${TERM}/align_wasp"
mkdir -p "${SAMPLE_OUT_DIR}"

OUT_PREFIX="${SAMPLE_OUT_DIR}/${ANIMAL_ID}_${TERM}."
OUT_BAM="${OUT_PREFIX}Aligned.sortedByCoord.out.bam"
OUT_COUNTS="${OUT_PREFIX}ReadsPerGene.out.tab"
STAR_STD_LOG="${SAMPLE_OUT_DIR}/${ANIMAL_ID}_${TERM}.STAR.stdout.log"


RG_ID="${ANIMAL_ID}_${TERM}"
RG_SM="${ANIMAL_ID}_${TERM}"
RG_LB="${ANIMAL_ID}_${TERM}"
RG_PL="ILLUMINA"
RG_PU="${ANIMAL_ID}_${TERM}"

if ! command -v STAR >/dev/null 2>&1; then
  echo "[ERROR] STAR not found in PATH" >&2
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

if [[ ! -d "$STAR_INDEX_DIR" ]]; then
  echo "[ERROR] STAR index dir not found: $STAR_INDEX_DIR" >&2
  exit 1
fi

if [[ ! -f "$VAR_VCF" ]]; then
  echo "[ERROR] VCF not found: $VAR_VCF" >&2
  exit 1
fi

if [[ -s "$OUT_BAM" && -s "$OUT_COUNTS" ]]; then
  echo "[skip] ${ANIMAL_ID} ${TERM} already finished."
  exit 0
fi

echo "[start] $(date '+%F %T') animal_id=$ANIMAL_ID term=$TERM"
echo "[info] THREADS=$THREADS"
echo "[info] R1=$R1"
echo "[info] R2=$R2"
echo "[info] STAR_INDEX_DIR=$STAR_INDEX_DIR"
echo "[info] VAR_VCF=$VAR_VCF"
echo "[info] OUT_PREFIX=$OUT_PREFIX"
echo "[info] RG_ID=$RG_ID"
echo "[info] RG_SM=$RG_SM"
echo "[info] RG_LB=$RG_LB"
echo "[info] RG_PL=$RG_PL"
echo "[info] RG_PU=$RG_PU"

ulimit -n 65535

STAR \
  --runThreadN "${THREADS}" \
  --genomeDir "${STAR_INDEX_DIR}" \
  --readFilesIn "${R1}" "${R2}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUT_PREFIX}" \
  --outSAMtype BAM SortedByCoordinate \
  --limitBAMsortRAM 120000000000 \
  --outSAMunmapped Within \
  --outSAMattributes NH HI AS nM vA vG vW \
  --outSAMattrRGline ID:${RG_ID} SM:${RG_SM} LB:${RG_LB} PL:${RG_PL} PU:${RG_PU} \
  --waspOutputMode SAMtag \
  --varVCFfile "${VAR_VCF}" \
  --quantMode GeneCounts \
  > "${STAR_STD_LOG}" 2>&1

echo "[done] $(date '+%F %T') animal_id=$ANIMAL_ID term=$TERM"
echo "[done] BAM=${OUT_BAM}"
echo "[done] Counts=${OUT_COUNTS}"
echo "[done] STAR stdout/stderr log=${STAR_STD_LOG}"
echo "[done] STAR main log=${OUT_PREFIX}Log.out"
echo "[done] STAR final log=${OUT_PREFIX}Log.final.out"
