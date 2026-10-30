#!/usr/bin/bash
set -euo pipefail


if [[ $# -ne 3 ]]; then
  echo "usage: /usr/bin/bash $0 <animal_id> <term> <alignment_run_dir>" >&2
  exit 1
fi

ANIMAL_ID="$1"
TERM="$2"
RUN_DIR="$3"

THREADS="${THREADS:-20}"
SAMTOOLS_EXTRA_THREADS="${SAMTOOLS_EXTRA_THREADS:-19}"
LIMIT_BAM_SORT_RAM="${LIMIT_BAM_SORT_RAM:-120000000000}"

MANIFEST="/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation/1.1star-index/20260820_005233/five_dogs_star_index_manifest.tsv"
INDEX_ROOT="/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation/1.1star-index/20260820_005233"
EXPECTED_INDEX_RUN_ID="20260820_005233"
EXPECTED_FASTA="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna"
EXPECTED_GTF="/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/GCF_000002285.3_CanFam3.1_genomic_chr.plus_five_only_new_gene.gtf"
OUTPUT_ROOT="/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation/1.2.1star"

STAR_BIN="/path/to/software/miniconda3/envs/dogs/bin/STAR"
SAMTOOLS_BIN="/path/to/software/miniconda3/envs/dogs/bin/samtools"
GZIP_BIN="/usr/bin/gzip"
ZCAT_BIN="/usr/bin/zcat"
AWK_BIN="/usr/bin/awk"
HEAD_BIN="/usr/bin/head"
TR_BIN="/usr/bin/tr"
MKDIR_BIN="/usr/bin/mkdir"
DATE_BIN="/usr/bin/date"

timestamp() {
  "${DATE_BIN}" '+%F %T'
}

log() {
  printf '[%s ] %s\n' "$(timestamp)" "$*"
}

[[ "${ANIMAL_ID}" != "" && "${TERM}" != "" ]] || { log "error: animal_id and term must not be empty"; exit 1; }
[[ "${THREADS}" =~ ^[0-9]+$ && "${THREADS}" -eq 20 ]] || { log "error: THREADS must be 20 as specified by this workflow"; exit 1; }
[[ "${SAMTOOLS_EXTRA_THREADS}" =~ ^[0-9]+$ && "${SAMTOOLS_EXTRA_THREADS}" -eq 19 ]] || { log "error: SAMTOOLS_EXTRA_THREADS must be 19, so the total samtools thread count is approximately 20"; exit 1; }
[[ "${LIMIT_BAM_SORT_RAM}" =~ ^[0-9]+$ && "${LIMIT_BAM_SORT_RAM}" -gt 0 ]] || { log "error: LIMIT_BAM_SORT_RAM must be a positive integer"; exit 1; }

run_parent="${RUN_DIR%/*}"
run_id="${RUN_DIR##*/}"
[[ "${RUN_DIR}" == /* && "${run_parent}" == "${OUTPUT_ROOT}" && "${run_id}" =~ ^[0-9]{8}_[0-9]{6}$ ]] || {
  log "error: alignment_run_dir must be a timestamped directory under OUTPUT_ROOT ${RUN_DIR}"
  exit 1
}

for required_file in "${MANIFEST}" "${EXPECTED_FASTA}" "${EXPECTED_GTF}"; do
  [[ -s "${required_file}" ]] || { log "error: required input is missing or empty ${required_file}"; exit 1; }
done
[[ -d "${RUN_DIR}" ]] || { log "error: batch run directory does not exist ${RUN_DIR}"; exit 1; }
[[ -x "${STAR_BIN}" ]] || { log "error: STAR not executable ${STAR_BIN}"; exit 1; }
[[ -x "${SAMTOOLS_BIN}" ]] || { log "error: samtools not executable ${SAMTOOLS_BIN}"; exit 1; }
[[ -x "${GZIP_BIN}" ]] || { log "error: gzip not executable ${GZIP_BIN}"; exit 1; }
[[ -x "${ZCAT_BIN}" ]] || { log "error: zcat not executable ${ZCAT_BIN}"; exit 1; }

EXPECTED_MANIFEST_HEADER=$'run_id\tanimal_id\tterm\tmax_read_length\tsjdbOverhang\tfq1_path\tfq2_path\tstar_index_dir\tgenome_fasta\tgtf_file\tannotation_scope\treadlen_source\thistorical_readlen_gtf'
ACTUAL_MANIFEST_HEADER="$(${HEAD_BIN} -n 1 "${MANIFEST}" | ${TR_BIN} -d '\r')"
[[ "${ACTUAL_MANIFEST_HEADER}" == "${EXPECTED_MANIFEST_HEADER}" ]] || {
  log "error: shared index manifest unexpected header: ${ACTUAL_MANIFEST_HEADER}"
  exit 1
}

MATCH_COUNT="$(${AWK_BIN} -F'\t' -v animal="${ANIMAL_ID}" -v term="${TERM}" '
  NR > 1 && $2 == animal && $3 == term {count++}
  END {print count + 0}
' "${MANIFEST}")"
[[ "${MATCH_COUNT}" -eq 1 ]] || {
  log "error: manifest in animal_id=${ANIMAL_ID} term=${TERM} of matching count is ${MATCH_COUNT}, expected 1"
  exit 1
}

IFS=$'\t' read -r MANIFEST_RUN_ID MAX_READ_LENGTH SJDB_OVERHANG R1 R2 STAR_INDEX_DIR GENOME_FASTA GTF_FILE ANNOTATION_SCOPE < <(
  "${AWK_BIN}" -F'\t' -v animal="${ANIMAL_ID}" -v term="${TERM}" '
    NR > 1 && $2 == animal && $3 == term {
      gsub(/\r/, "", $0)
      print $1 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11
      exit
    }
  ' "${MANIFEST}"
)

[[ "${MANIFEST_RUN_ID}" == "${EXPECTED_INDEX_RUN_ID}" ]] || { log "error: index manifest run ID does not satisfy expected ${MANIFEST_RUN_ID}"; exit 1; }
[[ "${ANNOTATION_SCOPE}" == "five-dogs-only" ]] || { log "error: annotation_scope is not five-dogs-only ${ANNOTATION_SCOPE}"; exit 1; }
[[ "${GENOME_FASTA}" == "${EXPECTED_FASTA}" ]] || { log "error: manifest reference FASTA does not satisfy five-only workflow contract ${GENOME_FASTA}"; exit 1; }
[[ "${GTF_FILE}" == "${EXPECTED_GTF}" ]] || { log "error: manifest GTF does not satisfy five-only workflow contract ${GTF_FILE}"; exit 1; }
[[ "${MAX_READ_LENGTH}" =~ ^[0-9]+$ && "${MAX_READ_LENGTH}" -gt 1 ]] || { log "error: manifest in max_read_length invalid"; exit 1; }
[[ "${SJDB_OVERHANG}" =~ ^[0-9]+$ && "${SJDB_OVERHANG}" -gt 0 ]] || { log "error: manifest in sjdbOverhang invalid"; exit 1; }
[[ "$((MAX_READ_LENGTH - 1))" -eq "${SJDB_OVERHANG}" ]] || { log "error: sjdbOverhang does not equal max_read_length-1"; exit 1; }

EXPECTED_STAR_INDEX_DIR="${INDEX_ROOT}/shared/sjdbOverhang_${SJDB_OVERHANG}/STARindex"
[[ "${STAR_INDEX_DIR}" == "${EXPECTED_STAR_INDEX_DIR}" ]] || {
  log "error: manifest index path does not match sjdbOverhang: ${STAR_INDEX_DIR}"
  exit 1
}

[[ -s "${R1}" ]] || { log "error: R1 missing or is empty ${R1}"; exit 1; }
[[ -s "${R2}" ]] || { log "error: R2 missing or is empty ${R2}"; exit 1; }
[[ -s "${STAR_INDEX_DIR}/Genome" ]] || { log "error: STAR index missing nonempty Genome ${STAR_INDEX_DIR}"; exit 1; }
[[ -s "${STAR_INDEX_DIR}/SA" ]] || { log "error: STAR index missing nonempty SA ${STAR_INDEX_DIR}"; exit 1; }
[[ -s "${STAR_INDEX_DIR}/SAindex" ]] || { log "error: STAR index missing nonempty SAindex ${STAR_INDEX_DIR}"; exit 1; }
[[ -s "${STAR_INDEX_DIR}/genomeParameters.txt" ]] || { log "error: STAR index missing nonempty genomeParameters.txt ${STAR_INDEX_DIR}"; exit 1; }

SAMPLE_OUT_DIR="${RUN_DIR}/${ANIMAL_ID}/${TERM}/align_star"
OUT_PREFIX="${SAMPLE_OUT_DIR}/${ANIMAL_ID}_${TERM}."
OUT_BAM="${OUT_PREFIX}Aligned.sortedByCoord.out.bam"
STAR_STD_LOG="${SAMPLE_OUT_DIR}/${ANIMAL_ID}_${TERM}.STAR.stdout.log"

if [[ -e "${SAMPLE_OUT_DIR}" ]]; then
  log "error: sample output directory already exists, refuse to overwrite or mixed run ${SAMPLE_OUT_DIR}"
  exit 1
fi
"${MKDIR_BIN}" -p "${SAMPLE_OUT_DIR}"

RG_ID="${ANIMAL_ID}_${TERM}"
RG_SM="${ANIMAL_ID}_${TERM}"
RG_LB="${ANIMAL_ID}_${TERM}"
RG_PL="ILLUMINA"
RG_PU="${ANIMAL_ID}_${TERM}"

log "start: five dogs five-only paired-end STAR alignment"
log "ALIGNMENT_RUN_ID=${run_id}"
log "INDEX_RUN_ID=${MANIFEST_RUN_ID}"
log "ANIMAL_ID=${ANIMAL_ID}"
log "TERM=${TERM}"
log "THREADS=${THREADS}"
log "SAMTOOLS_ADDITIONAL_THREADS=${SAMTOOLS_EXTRA_THREADS}"
log "LIMIT_BAM_SORT_RAM=${LIMIT_BAM_SORT_RAM}"
log "MAX_READ_LENGTH=${MAX_READ_LENGTH}"
log "SJDB_OVERHANG=${SJDB_OVERHANG}"
log "R1=${R1}"
log "R2=${R2}"
log "STAR_INDEX_DIR=${STAR_INDEX_DIR}"
log "GENOME_FASTA=${GENOME_FASTA}"
log "GTF_FILE=${GTF_FILE}"
log "ANNOTATION_SCOPE=${ANNOTATION_SCOPE}"
log "OUT_PREFIX=${OUT_PREFIX}"

log "stage 1/3: check paired-end FASTQ of gzip completeness"
"${GZIP_BIN}" -t "${R1}"
"${GZIP_BIN}" -t "${R2}"

ulimit -n 65535

log "stage 2/3: run paired-end STAR alignment and GeneCounts; STAR standard output write ${STAR_STD_LOG}"
"${STAR_BIN}" \
  --runThreadN "${THREADS}" \
  --genomeDir "${STAR_INDEX_DIR}" \
  --readFilesIn "${R1}" "${R2}" \
  --readFilesCommand "${ZCAT_BIN}" \
  --outFileNamePrefix "${OUT_PREFIX}" \
  --outSAMtype BAM SortedByCoordinate \
  --limitBAMsortRAM "${LIMIT_BAM_SORT_RAM}" \
  --outSAMattributes NH HI AS nM \
  --outSAMattrRGline ID:"${RG_ID}" SM:"${RG_SM}" LB:"${RG_LB}" PL:"${RG_PL}" PU:"${RG_PU}" \
  --quantMode GeneCounts \
  > "${STAR_STD_LOG}" 2>&1

log "stage 3/3: is coordinates sort BAM create BAI index"
"${SAMTOOLS_BIN}" index -@ "${SAMTOOLS_EXTRA_THREADS}" "${OUT_BAM}"

log "BAM=${OUT_BAM}"
log "BAM_INDEX=${OUT_BAM}.bai"
log "COUNTS=${OUT_PREFIX}ReadsPerGene.out.tab"
log "STAR_FINAL_LOG=${OUT_PREFIX}Log.final.out"
log "task completed / done"
log "execution_status=success"
log "reminder: execution completed does not imply BAM, BAI, GeneCounts, alignment rate or biological quality already passed; please in run after separately check result."
