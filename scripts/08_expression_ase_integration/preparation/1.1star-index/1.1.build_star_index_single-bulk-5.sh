#!/usr/bin/bash
set -euo pipefail


if [[ $# -ne 5 ]]; then
  echo "usage: /usr/bin/bash $0 <max_read_length> <sjdbOverhang> <run_dir> <genome_fasta> <gtf_file>" >&2
  exit 1
fi

MAX_READ_LENGTH="$1"
SJDB_OVERHANG="$2"
RUN_DIR="$3"
GENOME_FASTA="$4"
GTF_FILE="$5"

THREADS="${THREADS:-30}"
OUTPUT_ROOT="/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation/1.1star-index"
EXPECTED_GTF="/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/GCF_000002285.3_CanFam3.1_genomic_chr.plus_five_only_new_gene.gtf"
EXPECTED_FASTA="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna"
STAR_BIN="/path/to/software/miniconda3/envs/dogs/bin/STAR"

DATE_BIN="/usr/bin/date"
MKDIR_BIN="/usr/bin/mkdir"

timestamp() {
  "${DATE_BIN}" '+%F %T'
}

log() {
  printf '[%s ] %s\n' "$(timestamp)" "$*"
}

[[ "${MAX_READ_LENGTH}" =~ ^[0-9]+$ && "${MAX_READ_LENGTH}" -gt 1 ]] || { log "error: max_read_length must be an integer greater than 1"; exit 1; }
[[ "${SJDB_OVERHANG}" =~ ^[0-9]+$ && "${SJDB_OVERHANG}" -gt 0 ]] || { log "error: sjdbOverhang must be a positive integer"; exit 1; }
[[ "$((MAX_READ_LENGTH - 1))" -eq "${SJDB_OVERHANG}" ]] || { log "error: sjdbOverhang must be equal to max_read_length-1"; exit 1; }
[[ "${THREADS}" =~ ^[0-9]+$ && "${THREADS}" -gt 0 ]] || { log "error: THREADS must be a positive integer"; exit 1; }

run_parent="${RUN_DIR%/*}"
run_id="${RUN_DIR##*/}"
[[ "${RUN_DIR}" == /* && "${run_parent}" == "${OUTPUT_ROOT}" && "${run_id}" =~ ^[0-9]{8}_[0-9]{6}$ ]] || {
  log "error: run_dir must be a timestamped directory under OUTPUT_ROOT ${RUN_DIR}"
  exit 1
}
[[ "${GENOME_FASTA}" == "${EXPECTED_FASTA}" ]] || { log "error: reference FASTA does not satisfy five-only workflow contract ${GENOME_FASTA}"; exit 1; }
[[ "${GTF_FILE}" == "${EXPECTED_GTF}" ]] || { log "error: GTF does not satisfy five-only workflow contract ${GTF_FILE}"; exit 1; }
[[ -d "${RUN_DIR}" ]] || { log "error: batch run directory does not exist ${RUN_DIR}"; exit 1; }
[[ -s "${GENOME_FASTA}" ]] || { log "error: reference FASTA missing or is empty ${GENOME_FASTA}"; exit 1; }
[[ -s "${GTF_FILE}" ]] || { log "error: five-only GTF missing or is empty ${GTF_FILE}"; exit 1; }
[[ -x "${STAR_BIN}" ]] || { log "error: STAR not executable ${STAR_BIN}"; exit 1; }

GROUP_DIR="${RUN_DIR}/shared/sjdbOverhang_${SJDB_OVERHANG}"
STAR_INDEX_DIR="${GROUP_DIR}/STARindex"
INPUT_RECORD="${GROUP_DIR}/index_build_inputs.tsv"
STAR_LOG="${GROUP_DIR}/STAR_genomeGenerate.log"

if [[ -e "${GROUP_DIR}" ]]; then
  log "error: index group directory already exists, refuse to overwrite or mixed run ${GROUP_DIR}"
  exit 1
fi
"${MKDIR_BIN}" -p "${STAR_INDEX_DIR}"

printf 'run_id\tmax_read_length\tsjdbOverhang\tthreads\tstar_bin\tgenome_fasta\tgtf_file\tstar_index_dir\n' > "${INPUT_RECORD}"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "${run_id}" "${MAX_READ_LENGTH}" "${SJDB_OVERHANG}" "${THREADS}" \
  "${STAR_BIN}" "${GENOME_FASTA}" "${GTF_FILE}" "${STAR_INDEX_DIR}" >> "${INPUT_RECORD}"

log "start: build five dogs five-only shared STAR index"
log "RUN_ID=${run_id}"
log "MAX_READ_LENGTH=${MAX_READ_LENGTH}"
log "SJDB_OVERHANG=${SJDB_OVERHANG}"
log "THREADS=${THREADS}"
log "STAR_BIN=${STAR_BIN}"
log "GENOME_FASTA=${GENOME_FASTA}"
log "GTF_FILE=${GTF_FILE}"
log "STAR_INDEX_DIR=${STAR_INDEX_DIR}"
log "STAR_LOG=${STAR_LOG}"

log "run STAR genomeGenerate"
if ! (
  cd "${GROUP_DIR}"
  "${STAR_BIN}" \
    --runThreadN "${THREADS}" \
    --runMode genomeGenerate \
    --genomeDir "${STAR_INDEX_DIR}" \
    --genomeFastaFiles "${GENOME_FASTA}" \
    --sjdbGTFfile "${GTF_FILE}" \
    --sjdbOverhang "${SJDB_OVERHANG}" \
    > "${STAR_LOG}" 2>&1
); then
  log "error: STAR genomeGenerate failed, see ${STAR_LOG}"
  exit 1
fi

log "task completed / done"
log "execution_status=success"
log "reminder: execution completed does not imply index completeness, format, count or biological quality already passed; please in run after separately check result."
