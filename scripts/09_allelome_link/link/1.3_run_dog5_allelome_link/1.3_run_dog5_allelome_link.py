#!/usr/bin/python3

import csv
import os
import signal
import subprocess
import traceback
from collections import Counter
from contextlib import ExitStack, redirect_stdout, redirect_stderr
from datetime import datetime
from pathlib import Path

ROOT = Path('/path/to/project/analysis_20260817/5-dogs/7.Allelome.LINK')
STEP = '1.3_run_dog5_allelome_link'
SCRIPT_DIR = ROOT / 'scripts' / STEP
INPUT = ROOT / '1.2_build_dog5_allelome_locus_tables/20260907_202637'
RSCRIPT = '/path/to/software/miniconda3/envs/dogs/bin/Rscript'
LINK = '/path/to/software/Allelome.LINK/Allelome.LINK.R'
FIELDS = ['task_id', 'execution_status', 'exit_code', 'started_at', 'finished_at',
          'input_file', 'locus_metadata', 'output_path', 'task_log', 'reason']
child = None


def now():
    return datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def log(message):
    print(f'[{now()} ] {message}', flush=True)


def stop(signum, frame):

    raise InterruptedError(f'received_signal={signum}')


def cleanup_child():

    global child
    if child is None:
        return None

    previous = {s: signal.signal(s, signal.SIG_IGN)
                for s in (signal.SIGTERM, signal.SIGINT)}
    try:
        if child.poll() is None:
            try:
                os.killpg(child.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
        try:
            code = child.wait(timeout=10)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            code = child.wait()
        child = None
        return code
    finally:
        for s, handler in previous.items():
            signal.signal(s, handler)


def run_tasks(run_dir, writers, handles):
    global child
    marker = dict(line.split('=', 1) for line in
                  (INPUT / 'execution_completion.marker').read_text().splitlines() if '=' in line)
    if marker.get('execution_status') != 'success':
        raise ValueError('step two did not finish normally')
    with (INPUT / 'dog5_allelome_locus_table_manifest.tsv').open(newline='') as handle:
        tasks = list(csv.DictReader(handle, delimiter='\t'))
    log(f'[1/2] input={INPUT}; tasks={len(tasks)}; concurrent_tasks=1 threads_per_task=1(new script default)')
    log(f'Rscript={RSCRIPT}; LINK={LINK}; reads=10 bias=0.6 window_kb=500 duplicates=TRUE')
    env = dict(os.environ)
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        env[key] = '1'
    counts = Counter()
    interrupted = False
    for index, task in enumerate(tasks, 1):
        name = task['task_id']
        parent = run_dir / task['animal_id'] / task['tissue']
        output = parent / name
        task_log = parent / f'{name}.Rscript.log'
        row = dict(task_id=name, started_at=now(), input_file=task['locus_table'],
                   locus_metadata=task['locus_metadata'], output_path=str(output), task_log=str(task_log),
                   execution_status='', exit_code='', finished_at='', reason='')
        if interrupted or int(task['parameter_pass_rows']) < 2:
            row.update(execution_status='skipped', reason='batch_interrupted' if interrupted else 'fewer_than_two_biased_loci')
        else:
            log(f'[2/2] [{index}/{len(tasks)}] START {name}')
            try:
                parent.mkdir(parents=True, exist_ok=True)
                if output.exists():
                    raise FileExistsError(str(output))

                command = [RSCRIPT, '--vanilla', LINK, '--input', task['locus_table'],
                           '--name', name, '--output', str(parent), '--window-size', '500',
                           '--total-reads', '10', '--allelic-bias', '0.6', '--duplicates', 'TRUE']
                with task_log.open('x', encoding='utf-8') as capture:
                    child = subprocess.Popen(command, stdout=capture, stderr=subprocess.STDOUT,
                                             stdin=subprocess.DEVNULL, cwd=run_dir, env=env,
                                             start_new_session=True)
                    code = child.wait()
                child = None
                row.update(execution_status='success' if code == 0 else 'failed', exit_code=code,
                           reason='Rscript_exit_zero' if code == 0 else 'Rscript_nonzero_exit')
            except (InterruptedError, KeyboardInterrupt) as error:
                interrupted = True
                row.update(execution_status='failed', reason=str(error))
                code = cleanup_child()
                if code is not None:
                    row['exit_code'] = code
            except Exception as error:
                row.update(execution_status='failed', reason=str(error))
                code = cleanup_child()
                if code is not None:
                    row['exit_code'] = code
                traceback.print_exc()
        row['finished_at'] = now()
        status = row['execution_status']
        writers[status].writerow(row)
        handles[status].flush()
        counts[status] += 1
        log(f'[{index}/{len(tasks)}] {name}: {status}; reason={row["reason"]}')
    log(f'execution_tasks={len(tasks)} success={counts["success"]} failed={counts["failed"]} skipped={counts["skipped"]}')
    if interrupted or counts['failed']:
        raise RuntimeError('batch has failed or is interrupted; please check failed.tsv and task log')


def main():
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_dir = SCRIPT_DIR / 'log'
    log_dir.mkdir(parents=True, exist_ok=True)
    run_dir = ROOT / STEP / stamp
    with ExitStack() as stack:
        master = stack.enter_context((log_dir / f'{STEP}.{stamp}.log').open('x', encoding='utf-8', buffering=1))
        stack.enter_context(redirect_stdout(master))
        stack.enter_context(redirect_stderr(master))
        handles, writers = {}, {}
        for status in ('success', 'failed', 'skipped'):
            handle = stack.enter_context((log_dir / f'{STEP}.{stamp}.{status}.tsv').open('x', encoding='utf-8', newline=''))
            handles[status] = handle
            writers[status] = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t', lineterminator='\n')
            writers[status].writeheader()
            handle.flush()
        signal.signal(signal.SIGTERM, stop)
        signal.signal(signal.SIGINT, stop)
        try:
            run_dir.mkdir(parents=True, exist_ok=False)
            log(f'run_id={stamp}; output_dir={run_dir}')
            run_tasks(run_dir, writers, handles)
        except BaseException:
            cleanup_child()
            traceback.print_exc()
            log('execution_status=failed; no generate batch completion marker.')
            raise
        with (run_dir / 'execution_completion.marker').open('x', encoding='utf-8') as handle:
            handle.write(f'run_id={stamp}\nexecution_status=success\ninput_run={INPUT}\n')
        log('task completed / done; execution_status=success')
        log('please independent review output completeness, format, link count and biological quality; execution success does not equal result passed QC.')


if __name__ == '__main__':
    main()
