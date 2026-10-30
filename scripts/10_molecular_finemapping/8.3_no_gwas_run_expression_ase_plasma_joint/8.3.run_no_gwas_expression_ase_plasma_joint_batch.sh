#!/usr/bin/env bash
set -euo pipefail


SCRIPT_DIR="/path/to/repository/scripts/10_molecular_finemapping/8.3_no_gwas_run_expression_ase_plasma_joint"
PLASMA_PYTHON="/path/to/software/environments/plasma_py27/bin/python2.7"
PLASMA_SCRIPT="/path/to/software/plasma/run_plasma.py"

PREPARE_RUN_DIR="/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.2_no_gwas_prepare_expression_ase_plasma_inputs/20260821_181104"
MANIFEST="${PREPARE_RUN_DIR}/no_gwas_plasma_candidate_blocks.tsv"

OUTPUT_ROOT="/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.3_no_gwas_run_expression_ase_plasma_joint"
RUN_ID="${RUN_ID:-$(/usr/bin/date '+%Y%m%d_%H%M%S')}"
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"
RESULT_DIR="${RUN_DIR}/plasma_results"
TASKS_FILE="${RUN_DIR}/8.3.no_gwas_expression_ase_plasma_joint.tasks.tsv"
RUN_SUMMARY_FILE="${RUN_DIR}/8.3.no_gwas_expression_ase_plasma_joint.run_summary.tsv"
COMPLETION_MARKER="${RUN_DIR}/execution_completion.marker"
INTERRUPTED_MARKER="${RUN_DIR}/execution_interrupted.marker"
LAUNCHER_STATE_FILE="${RUN_DIR}/launcher_state.tsv"
STATUS_LOCK="${RUN_DIR}/execution_status.lock"

LOG_DIR="${SCRIPT_DIR}/log"
MASTER_LOG="${LOG_DIR}/8.3.no_gwas_expression_ase_plasma_joint.${RUN_ID}.log"
SUCCESS_FILE="${LOG_DIR}/8.3.no_gwas_expression_ase_plasma_joint.${RUN_ID}.success.tsv"
FAILED_FILE="${LOG_DIR}/8.3.no_gwas_expression_ase_plasma_joint.${RUN_ID}.failed.tsv"
SKIPPED_FILE="${LOG_DIR}/8.3.no_gwas_expression_ase_plasma_joint.${RUN_ID}.skipped.tsv"
LAUNCHER_LOCK="${LOG_DIR}/8.3.no_gwas_expression_ase_plasma_joint.launcher.lock"

PARALLEL_JOBS=17
THREADS_PER_TASK=1
HEARTBEAT_SECONDS=300
MAX_CAUSAL="${MAX_CAUSAL:-1}"
CONFIDENCE="${CONFIDENCE:-0.95}"
SEARCH_MODE="${SEARCH_MODE:-exhaustive}"
OVERDISPERSION_GLOBAL="${OVERDISPERSION_GLOBAL:-0}"


export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

declare -A TASK_BLOCKS=()
TERMINATING=0
NORMAL_COMPLETION=0

log() {
  printf '[%s ] %s\n' "$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')" "$*"
}

die() {
  log "[ERROR] $*" >&2
  exit 1
}

require_file() {
  local path="$1"
  local label="$2"
  [[ -s "${path}" ]] || die "${label} does not exist or is empty: ${path}"
}

require_dir() {
  local path="$1"
  local label="$2"
  [[ -d "${path}" ]] || die "${label} directory does not exist: ${path}"
}

validate_integer() {
  local value="$1"
  local label="$2"
  [[ "${value}" =~ ^[0-9]+$ ]] || die "${label} must be nonnegative integer, current value: ${value}"
}

validate_probability() {
  local value="$1"
  local label="$2"
  /usr/bin/awk -v value="${value}" -v label="${label}" '
    BEGIN {
      if (value !~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ || value <= 0 || value >= 1) {
        printf("[ERROR] %s must be a numeric value strictly between 0 and 1; received: %s\n", label, value) > "/dev/stderr"
        exit 1
      }
    }
  ' || exit 1
}

validate_overdispersion() {
  local value="$1"
  /usr/bin/awk -v value="${value}" '
    BEGIN {
      if (value !~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ || value < 0 || value >= 1) {
        printf("[ERROR] OVERDISPERSION_GLOBAL must be a numeric value in [0,1); received: %s\n", value) > "/dev/stderr"
        exit 1
      }
    }
  ' || exit 1
}

prepare_run() {
  /usr/bin/mkdir -p "${LOG_DIR}" "${OUTPUT_ROOT}"

  [[ "${RUN_ID}" =~ ^[0-9]{8}_[0-9]{6}$ ]] || die "RUN_ID format must be YYYYMMDD_HHMMSS: ${RUN_ID}"
  exec 9>"${LAUNCHER_LOCK}"
  /usr/bin/flock -n 9 || die "Another instance of this step holds the lock: ${LAUNCHER_LOCK}"

  for path in "${RUN_DIR}" "${MASTER_LOG}" "${SUCCESS_FILE}" "${FAILED_FILE}" "${SKIPPED_FILE}"; do
    [[ ! -e "${path}" ]] || die "Refusing to overwrite an existing run; target already exists: ${path}"
  done

  /usr/bin/mkdir "${RUN_DIR}"
  /usr/bin/mkdir "${RESULT_DIR}"
  : > "${MASTER_LOG}"
  exec >> "${MASTER_LOG}" 2>&1
}

status_header() {
  printf 'run_id\tblock_id\tselected_sample_count\tmarker_count\ttotal_exp_unique_count\tanalysis_mode\tstatus\texit_code\tstarted_at\tfinished_at\tinput_dir\toutput_path\ttask_log\treason\n'
}

write_status_header() {
  local status_file="$1"
  status_header > "${status_file}"
}

validate_config() {
  require_file "${PLASMA_PYTHON}" "PLASMA Python interpreter"
  [[ -x "${PLASMA_PYTHON}" ]] || die "PLASMA Python interpreter not executable: ${PLASMA_PYTHON}"
  require_file "${PLASMA_SCRIPT}" "PLASMA main program"
  require_file "${MANIFEST}" "8.2 candidate block manifest"
  require_file "${PREPARE_RUN_DIR}/execution_completion.marker" "8.2 execution completion marker"
  require_dir "${PREPARE_RUN_DIR}/plasma_inputs" "8.2 PLASMA input"

  validate_integer "${MAX_CAUSAL}" "MAX_CAUSAL"
  [[ "${MAX_CAUSAL}" -gt 0 ]] || die "MAX_CAUSAL must be greater than or equal to 1"
  validate_probability "${CONFIDENCE}" "CONFIDENCE"
  validate_overdispersion "${OVERDISPERSION_GLOBAL}"
  [[ "${SEARCH_MODE}" == "exhaustive" || "${SEARCH_MODE}" == "shotgun" ]] || \
    die "SEARCH_MODE must be exhaustive or shotgun"
}

write_tasks() {
  /usr/bin/awk -F'\t' '
    BEGIN { OFS = "\t" }
    NR == 1 {
      for (i = 1; i <= NF; i++) h[$i] = i
      required[1] = "block_id"
      required[2] = "plasma_input_allowed"
      required[3] = "input_dir"
      required[4] = "selected_sample_count"
      required[5] = "common_marker_count"
      required[6] = "total_exp_values"
      for (i = 1; i <= 6; i++) {
        if (!(required[i] in h)) {
          printf("[ERROR] manifest missing required columns: %s\n", required[i]) > "/dev/stderr"
          exit 2
        }
      }
      print "block_id", "input_dir", "selected_sample_count", "common_marker_count", "total_exp_unique_count"
      next
    }
    $h["plasma_input_allowed"] == "1" && $h["input_dir"] != "" {
      delete seen
      unique_count = 0
      value_count = split($h["total_exp_values"], values, ",")
      for (j = 1; j <= value_count; j++) {
        if (!(values[j] in seen)) {
          seen[values[j]] = 1
          unique_count++
        }
      }
      print $h["block_id"], $h["input_dir"], $h["selected_sample_count"], $h["common_marker_count"], unique_count
    }
  ' "${MANIFEST}" > "${TASKS_FILE}"
}

reset_status_files() {
  write_status_header "${SUCCESS_FILE}"
  write_status_header "${FAILED_FILE}"
  write_status_header "${SKIPPED_FILE}"
}

load_task_blocks() {
  local block_id _rest

  while IFS=$'\t' read -r block_id _rest; do
    [[ "${block_id}" == "block_id" ]] && continue
    [[ -z "${block_id}" ]] && continue
    [[ -z "${TASK_BLOCKS[${block_id}]+x}" ]] || die "task table contains duplicate block_id: ${block_id}"
    TASK_BLOCKS["${block_id}"]=1
  done < "${TASKS_FILE}"
}

write_launcher_state() {
  local execution_status="$1"
  local signal_name="${2:-}"
  {
    printf 'key\tvalue\n'
    printf 'run_id\t%s\n' "${RUN_ID}"
    printf 'launcher_pid\t%s\n' "$$"
    printf 'execution_status\t%s\n' "${execution_status}"
    printf 'signal\t%s\n' "${signal_name}"
    printf 'updated_at\t%s\n' "$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')"
  } > "${LAUNCHER_STATE_FILE}"
}

handle_signal() {
  local signal_name="$1"
  local exit_code="$2"
  local -a active_pids

  [[ "${TERMINATING}" -eq 0 ]] || return
  TERMINATING=1
  log "[interrupt] received ${signal_name}, stop current launcher and its active task"
  write_launcher_state "interrupted" "${signal_name}"
  {
    printf 'run_id=%s\n' "${RUN_ID}"
    printf 'execution_status=interrupted\n'
    printf 'signal=%s\n' "${signal_name}"
    printf 'interrupted_at=%s\n' "$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')"
  } > "${INTERRUPTED_MARKER}"

  mapfile -t active_pids < <(jobs -pr)
  if [[ "${#active_pids[@]}" -gt 0 ]]; then
    kill -TERM "${active_pids[@]}" 2>/dev/null || true
  fi
  wait || true
  log "execution_status=interrupted"
  exit "${exit_code}"
}

handle_exit() {
  local exit_code="$1"
  if [[ "${exit_code}" -ne 0 && "${TERMINATING}" -eq 0 && "${NORMAL_COMPLETION}" -eq 0 ]]; then
    write_launcher_state "failed" ""
    log "execution_status=failed; exit_code=${exit_code}"
  fi
}

write_command_file() {
  local command_file="$1"
  shift
  printf '%q' "$1" > "${command_file}"
  shift
  local arg
  for arg in "$@"; do
    printf ' %q' "${arg}" >> "${command_file}"
  done
  printf '\n' >> "${command_file}"
}

append_status() {
  local status_file="$1"
  shift
  {
    /usr/bin/flock -x 200
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$@" >> "${status_file}"
  } 200>"${STATUS_LOCK}"
}

run_one_block() {
  local block_id="$1"
  local input_dir="$2"
  local selected_sample_count="$3"
  local marker_count="$4"
  local total_exp_unique_count="$5"
  local result_dir="${RESULT_DIR}/${block_id}"
  local stdout_log="${result_dir}/run.stdout.log"
  local stderr_log="${result_dir}/run.stderr.log"
  local started_at finished_at rc required reason plasma_pid

  started_at="$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')"
  for required in hap_A.txt counts_A.txt hap_B.txt counts_B.txt total_exp.txt markers.tsv samples.tsv; do
    if [[ ! -s "${input_dir}/${required}" ]]; then
      finished_at="$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')"
      reason="missing_input_${required}"
      append_status "${FAILED_FILE}" "${RUN_ID}" "${block_id}" "${selected_sample_count}" "${marker_count}" \
        "${total_exp_unique_count}" "not_started" "failed" "1" "${started_at}" "${finished_at}" \
        "${input_dir}" "" "" "${reason}"
      return 0
    fi
  done

  if [[ "${total_exp_unique_count}" -lt 2 ]]; then
    finished_at="$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')"
    append_status "${SKIPPED_FILE}" "${RUN_ID}" "${block_id}" "${selected_sample_count}" "${marker_count}" \
      "${total_exp_unique_count}" "not_run" "skipped" "0" "${started_at}" "${finished_at}" \
      "${input_dir}" "" "" "constant_total_expression_not_joint_estimable"
    return 0
  fi

  if [[ -e "${result_dir}" ]]; then
    finished_at="$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')"
    append_status "${FAILED_FILE}" "${RUN_ID}" "${block_id}" "${selected_sample_count}" "${marker_count}" \
      "${total_exp_unique_count}" "not_started" "failed" "1" "${started_at}" "${finished_at}" \
      "${input_dir}" "" "" "result_dir_already_exists"
    return 0
  fi

  if ! /usr/bin/mkdir "${result_dir}"; then
    finished_at="$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')"
    append_status "${FAILED_FILE}" "${RUN_ID}" "${block_id}" "${selected_sample_count}" "${marker_count}" \
      "${total_exp_unique_count}" "not_started" "failed" "1" "${started_at}" "${finished_at}" \
      "${input_dir}" "" "" "cannot_create_result_dir"
    return 0
  fi

  local -a cmd=(
    /usr/bin/env
    "OPENBLAS_NUM_THREADS=${THREADS_PER_TASK}"
    "OMP_NUM_THREADS=${THREADS_PER_TASK}"
    "MKL_NUM_THREADS=${THREADS_PER_TASK}"
    "NUMEXPR_NUM_THREADS=${THREADS_PER_TASK}"
    "VECLIB_MAXIMUM_THREADS=${THREADS_PER_TASK}"
    "${PLASMA_PYTHON}"
    "${PLASMA_SCRIPT}"
    "${input_dir}/hap_A.txt"
    "${input_dir}/counts_A.txt"
    "${input_dir}/hap_B.txt"
    "${input_dir}/counts_B.txt"
    "${result_dir}"
    "--total_exp_path" "${input_dir}/total_exp.txt"
    "--search_mode" "${SEARCH_MODE}"
    "--max_causal" "${MAX_CAUSAL}"
    "--confidence" "${CONFIDENCE}"
    "--overdispersion_global" "${OVERDISPERSION_GLOBAL}"
  )

  write_command_file "${result_dir}/command.txt" "${cmd[@]}"

  set +e
  "${cmd[@]}" > "${stdout_log}" 2> "${stderr_log}" &
  plasma_pid=$!
  trap 'kill -TERM "${plasma_pid}" 2>/dev/null || true' TERM INT
  wait "${plasma_pid}"
  rc=$?
  trap - TERM INT
  set -e
  finished_at="$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')"

  if [[ "${rc}" -eq 0 ]]; then
    append_status "${SUCCESS_FILE}" "${RUN_ID}" "${block_id}" "${selected_sample_count}" "${marker_count}" \
      "${total_exp_unique_count}" "joint_expression_ase" "success" "0" "${started_at}" "${finished_at}" \
      "${input_dir}" "${result_dir}" "${stderr_log}" ""
  else
    append_status "${FAILED_FILE}" "${RUN_ID}" "${block_id}" "${selected_sample_count}" "${marker_count}" \
      "${total_exp_unique_count}" "joint_expression_ase" "failed" "${rc}" "${started_at}" "${finished_at}" \
      "${input_dir}" "${result_dir}" "${stderr_log}" "plasma_exit_nonzero"
  fi

  return 0
}

run_tasks() {
  local total_tasks pending_tasks active_jobs launched
  local last_heartbeat_epoch now_epoch success_count failed_count skipped_count
  local -a active_pids
  total_tasks=${#TASK_BLOCKS[@]}
  [[ "${total_tasks}" -gt 0 ]] || die "manifest contains no schedulable task"
  pending_tasks=${total_tasks}

  log "[stage:tasks] tasks=${total_tasks}; tasks_file=${TASKS_FILE}"
  log "[stage:run] parallel_jobs=${PARALLEL_JOBS}; threads_per_task=${THREADS_PER_TASK}"
  log "[stage:run] PLASMA mode=expression+ASE joint; total_exp=enabled; as_only=disabled"
  log "[stage:run] Blocks with constant total expression will be skipped; this step does not fall back to ASE-only analysis."

  launched=0
  last_heartbeat_epoch="$(/usr/bin/date '+%s')"
  while IFS=$'\t' read -r block_id input_dir selected_sample_count marker_count total_exp_unique_count; do
    [[ "${block_id}" == "block_id" ]] && continue
    run_one_block "${block_id}" "${input_dir}" "${selected_sample_count}" "${marker_count}" \
      "${total_exp_unique_count}" &
    launched=$((launched + 1))

    if [[ "${launched}" -eq 1 || $((launched % 500)) -eq 0 || "${launched}" -eq "${pending_tasks}" ]]; then
      log "[progress] dispatched=${launched}/${pending_tasks}"
    fi

    while true; do
      mapfile -t active_pids < <(jobs -rp)
      active_jobs=${#active_pids[@]}
      [[ "${active_jobs}" -lt "${PARALLEL_JOBS}" ]] && break

      now_epoch="$(/usr/bin/date '+%s')"
      if [[ $((now_epoch - last_heartbeat_epoch)) -ge "${HEARTBEAT_SECONDS}" ]]; then
        success_count=$(( $(/usr/bin/wc -l < "${SUCCESS_FILE}") - 1 ))
        failed_count=$(( $(/usr/bin/wc -l < "${FAILED_FILE}") - 1 ))
        skipped_count=$(( $(/usr/bin/wc -l < "${SKIPPED_FILE}") - 1 ))
        log "[heartbeat] dispatched=${launched}/${pending_tasks}; active=${active_jobs}; success=${success_count}; failed=${failed_count}; skipped=${skipped_count}"
        last_heartbeat_epoch="${now_epoch}"
        write_launcher_state "running" ""
      fi
      /usr/bin/sleep 2
    done
  done < "${TASKS_FILE}"

  [[ "${launched}" -eq "${pending_tasks}" ]] || \
    die "dispatched task count differs from the expected pending task count: launched=${launched}, pending=${pending_tasks}"

  while true; do
    mapfile -t active_pids < <(jobs -rp)
    active_jobs=${#active_pids[@]}
    [[ "${active_jobs}" -eq 0 ]] && break

    now_epoch="$(/usr/bin/date '+%s')"
    if [[ $((now_epoch - last_heartbeat_epoch)) -ge "${HEARTBEAT_SECONDS}" ]]; then
      success_count=$(( $(/usr/bin/wc -l < "${SUCCESS_FILE}") - 1 ))
      failed_count=$(( $(/usr/bin/wc -l < "${FAILED_FILE}") - 1 ))
      skipped_count=$(( $(/usr/bin/wc -l < "${SKIPPED_FILE}") - 1 ))
      log "[heartbeat] dispatched=${launched}/${pending_tasks}; active=${active_jobs}; success=${success_count}; failed=${failed_count}; skipped=${skipped_count}"
      last_heartbeat_epoch="${now_epoch}"
      write_launcher_state "running" ""
    fi
    /usr/bin/sleep 2
  done

  wait
  log "[stage:run] All dispatched tasks have returned."
}

write_run_summary() {
  local tasks success failed skipped terminal_total
  tasks=$(( $(/usr/bin/wc -l < "${TASKS_FILE}") - 1 ))
  success=$(( $(/usr/bin/wc -l < "${SUCCESS_FILE}") - 1 ))
  failed=$(( $(/usr/bin/wc -l < "${FAILED_FILE}") - 1 ))
  skipped=$(( $(/usr/bin/wc -l < "${SKIPPED_FILE}") - 1 ))
  terminal_total=$((success + failed + skipped))

  {
    printf 'key\tvalue\n'
    printf 'run_id\t%s\n' "${RUN_ID}"
    printf 'run_mode\tfresh_full_run\n'
    printf 'tasks\t%s\n' "${tasks}"
    printf 'success\t%s\n' "${success}"
    printf 'failed\t%s\n' "${failed}"
    printf 'skipped\t%s\n' "${skipped}"
    printf 'terminal_total\t%s\n' "${terminal_total}"
    printf 'parallel_jobs\t%s\n' "${PARALLEL_JOBS}"
    printf 'threads_per_task\t%s\n' "${THREADS_PER_TASK}"
    printf 'heartbeat_seconds\t%s\n' "${HEARTBEAT_SECONDS}"
    printf 'search_mode\t%s\n' "${SEARCH_MODE}"
    printf 'max_causal\t%s\n' "${MAX_CAUSAL}"
    printf 'confidence\t%s\n' "${CONFIDENCE}"
    printf 'overdispersion_global\t%s\n' "${OVERDISPERSION_GLOBAL}"
    printf 'plasma_python\t%s\n' "${PLASMA_PYTHON}"
    printf 'plasma_script\t%s\n' "${PLASMA_SCRIPT}"
    printf 'prepare_manifest\t%s\n' "${MANIFEST}"
    printf 'result_dir\t%s\n' "${RESULT_DIR}"
  } > "${RUN_SUMMARY_FILE}"

  log "[summary] tasks=${tasks}; success=${success}; failed=${failed}; skipped=${skipped}"
  [[ "${terminal_total}" -eq "${tasks}" ]] || die "task terminal state count mismatch: tasks=${tasks}, terminal_total=${terminal_total}"
  [[ "${failed}" -eq 0 ]] || die "${failed} task(s) failed; no execution completion marker will be written."
}

write_completion_marker() {
  {
    printf 'run_id=%s\n' "${RUN_ID}"
    printf 'execution_status=success\n'
    printf 'finished_at=%s\n' "$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')"
  } > "${COMPLETION_MARKER}"
}

main() {
  prepare_run
  log "[start] start run no-GWAS PLASMA expression + ASE joint"
  log "[config] run_id=${RUN_ID}"
  log "[config] run_mode=fresh_full_run"
  log "[config] PLASMA_PYTHON=${PLASMA_PYTHON}"
  log "[config] PLASMA_SCRIPT=${PLASMA_SCRIPT}"
  log "[config] MANIFEST=${MANIFEST}"
  log "[config] RUN_DIR=${RUN_DIR}"
  log "[config] parallel_jobs=${PARALLEL_JOBS}; threads_per_task=${THREADS_PER_TASK}"
  log "[config] search_mode=${SEARCH_MODE}; max_causal=${MAX_CAUSAL}; confidence=${CONFIDENCE}; overdispersion_global=${OVERDISPERSION_GLOBAL}"

  write_launcher_state "running" ""
  trap 'handle_signal TERM 143' TERM
  trap 'handle_signal INT 130' INT
  trap 'handle_exit $?' EXIT

  validate_config
  log "[stage:prepare] configuration and upstream execution completion marker passed check"
  write_tasks
  reset_status_files
  load_task_blocks

  run_tasks
  write_run_summary

  log "task completed / done"
  log "execution_status=success"
  log "[reminder] Successful execution does not imply that output completeness, format, counts, or biological quality meet requirements; review the outputs separately."
  write_completion_marker
  write_launcher_state "success" ""
  NORMAL_COMPLETION=1
  trap - EXIT TERM INT
}

main "$@"
