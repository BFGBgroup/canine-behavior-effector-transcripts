#!/usr/bin/python3

import csv
import math
import signal
import traceback
from collections import Counter, defaultdict
from contextlib import redirect_stdout, redirect_stderr
from datetime import datetime
from functools import lru_cache
from pathlib import Path

ROOT = Path('/path/to/project/analysis_20260817/5-dogs/7.Allelome.LINK')
STEP = '1.2_build_dog5_allelome_locus_tables'
SCRIPT_DIR = ROOT / 'scripts' / STEP
AUDIT = ROOT / '1.1_validate_dog5_allelome_inputs/20260907_195800'
GENE_ROOT = Path('/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation/1.5ase-cnv-clean-gene-id-haplotype-table/20260821_161805')
LINK_FIELDS = ['chr', 'start', 'end', 'name', 'A1_reads', 'A2_reads',
               'total_reads', 'allelic_score', 'allelic_ratio']
EXTRA_FIELDS = ['input_file', 'input_line', 'assembly', 'start0', 'end0',
                'coordinate_source', 'allelic_score', 'allelic_ratio_raw',
                'allelic_ratio', 'passes_link_parameters']


def log_message(message):

    print(f'[{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} ] {message}', flush=True)


def read_rows(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def write_rows(path, fields, rows):
    with path.open('x', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


@lru_cache(maxsize=131072)
def score(h1, h2):

    total = h1 + h2
    if total < 10:
        return 0.0
    term = -total * math.log(2.0)
    log_sum = term
    for k in range(1, min(h1, h2) + 1):
        term += math.log(total - k + 1) - math.log(k)
        maximum = max(log_sum, term)
        log_sum = maximum + math.log(math.exp(log_sum - maximum) + math.exp(term - maximum))
    value = log_sum / math.log(10.0)
    return round(-value if h1 > h2 else value, 3)


def build(run_dir):
    log_message('[1/3] Read the completed input audit; concurrent_tasks=1 threads_per_task=1(serial settings retained)')
    log_message(f'input_audit_run={AUDIT}; gene_haplotype_run={GENE_ROOT}')
    marker = dict(line.split('=', 1) for line in
                  (AUDIT / 'execution_completion.marker').read_text().splitlines() if '=' in line)
    summary = {r['metric']: r['value'] for r in read_rows(AUDIT / 'dog5_allelome_validation_summary.tsv')}
    if marker.get('execution_status') != 'success' or summary.get('audit_status') != 'PASS':
        raise ValueError('The specified step-one run did not pass; review the inputs before proceeding')
    samples = list(read_rows(AUDIT / 'dog5_allelome_sample_manifest.tsv'))
    tasks = list(read_rows(AUDIT / 'dog5_allelome_chromosome_task_manifest.tsv'))
    reference = {}
    for row in read_rows(AUDIT / 'dog5_allelome_locus_reference.tsv'):
        key = row['gene_id']
        if key in reference:
            raise ValueError(f'The interval reference contains a duplicate GeneID: {key}')
        reference[key] = row
    task_map = defaultdict(list)
    for task in tasks:
        task_map[task['sample_key']].append(task)
    manifests = []
    totals = Counter()
    log_message(f'[2/3] Build nine-column tables for {len(samples)} samples and {len(tasks)} chromosome tasks')
    for index, sample in enumerate(samples, 1):
        if sample['audit_status'] != 'PASS':
            raise ValueError(f"sample failed input review: {sample['sample_key']}")
        source = Path(sample['gene_haplotype_tsv'])
        if not source.is_relative_to(GENE_ROOT):
            raise ValueError(f'Input does not belong to the specified 1.5 run: {source}')
        log_message(f"[{index}/{len(samples)}] Start reading and processing: {source}")
        grouped = defaultdict(list)
        with source.open(encoding='utf-8-sig', newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            source_fields = list(reader.fieldnames)
            for line, row in enumerate(reader, 2):
                ref = reference[row['gene_id']]
                if (row['chromosome'] != ref['chromosome'] or
                        row['record_type'] != ref['record_type'] or
                        row['gene_strand'] != ref['gene_strand'] or
                        row['sample'] != sample['animal_id'] or row['tissue'] != sample['tissue']):
                    raise ValueError(f'Input does not match the interval/sample list: {source}:{line}')
                h1, h2, total = (int(row[f]) for f in ('hap1_count_sum', 'hap2_count_sum', 'total_count_sum'))
                if min(h1, h2) < 0 or total != h1 + h2:
                    raise ValueError(f'haplotype count mismatch: {source}:{line}')
                start, end = int(ref['start0']), int(ref['end0'])
                ratio_raw = h1 / total if total else 0.0
                ratio = round(ratio_raw, 2)
                value = score(h1, h2)
                passed = total >= 10 and (ratio >= 0.6 or ratio <= 0.4)
                link = dict(zip(LINK_FIELDS, [row['chromosome'], start, end, row['group_id'],
                                              h1, h2, total, value, ratio]))
                metadata = dict(row)
                metadata.update(input_file=str(source), input_line=line, assembly=ref['assembly'],
                                start0=start, end0=end, coordinate_source=ref['coordinate_source'],
                                allelic_score=value, allelic_ratio_raw=ratio_raw, allelic_ratio=ratio,
                                passes_link_parameters='yes' if passed else 'no')
                grouped[row['chromosome']].append((link, metadata))
        sample_tasks = task_map[sample['sample_key']]
        if set(grouped) - {t['chromosome'] for t in sample_tasks}:
            raise ValueError(f"chromosome is absent from the input task list: {sample['sample_key']}")
        directory = run_dir / sample['animal_id'] / sample['tissue']
        directory.mkdir(parents=True, exist_ok=False)
        for task in sample_tasks:
            if task['audit_status'] != 'PASS' or task['gene_haplotype_tsv'] != str(source):
                raise ValueError(f"chromosome input task mismatch: {task['task_id']}")
            records = sorted(grouped[task['chromosome']], key=lambda x: (x[0]['start'], x[0]['end'], x[0]['name']))
            prefix = f"{sample['sample_key']}.{task['chromosome']}"
            table = directory / f'{prefix}.locus_table.tsv'
            metadata_path = directory / f'{prefix}.locus_metadata.tsv'
            write_rows(table, LINK_FIELDS, (r[0] for r in records))
            write_rows(metadata_path, source_fields + EXTRA_FIELDS, (r[1] for r in records))
            counts = Counter(r[1]['record_type'] for r in records)
            biased = sum(r[1]['passes_link_parameters'] == 'yes' for r in records)
            multi = sum(int(r[1]['multi_gene_variant_count']) > 0 for r in records)

            manifests.append(dict(task_id=task['task_id'], sample_key=sample['sample_key'],
                                  animal_id=sample['animal_id'], tissue=sample['tissue'],
                                  chromosome=task['chromosome'], input_file=str(source),
                                  locus_table=str(table), locus_metadata=str(metadata_path),
                                  written_rows=len(records), known_gene_rows=counts['known_gene'],
                                  manual_new_gene_rows=counts['manual_new_gene'],
                                  unannotated_variant_rows=counts['unannotated_variant'],
                                  parameter_pass_rows=biased, rows_with_multi_gene_variants=multi,
                                  link_candidate_status='ENOUGH_BIASED_LOCI' if biased >= 2 else 'FEWER_THAN_TWO_BIASED_LOCI'))
            totals.update(counts)
            totals['written_rows'] += len(records)
            totals['parameter_pass_rows'] += biased
            totals['rows_with_multi_gene_variants'] += multi
        log_message(f"[{index}/{len(samples)}] {sample['sample_key']} write {len(sample_tasks)} input sets")
    log_message('[3/3] write task manifest and parameter description(not a post-run QC)')
    write_rows(run_dir / 'dog5_allelome_locus_table_manifest.tsv', list(manifests[0]), manifests)
    totals.update(sample_count=len(samples), chromosome_task_count=len(manifests))
    write_rows(run_dir / 'dog5_allelome_locus_table_summary.tsv', ['metric', 'value'],
               ({'metric': k, 'value': v} for k, v in totals.items()))
    parameters = dict(input_audit_run=str(AUDIT), gene_haplotype_run=str(GENE_ROOT),
                      coordinate_convention='0-based_half-open', count_unit='ASEReadCounter_site_allele_support_sum',
                      allele_orientation='A1=sample_HAP1;A2=sample_HAP2;not_parental_origin',
                      min_total_reads=10, ratio_digits=2, bias_lower=0.4, bias_upper=0.6,
                      filtering='no_rows_removed;thresholds_flagged_only',
                      score='signed_log10_binomial_lower_tail;log_space;rounded_3_digits',
                      concurrent_tasks=1, threads_per_task=1)
    write_rows(run_dir / 'parameter_contract.tsv', ['parameter', 'value'],
               ({'parameter': k, 'value': v} for k, v in parameters.items()))


def interrupted(signum, frame):
    raise InterruptedError(f'received termination signal {signum}')


def main():
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_dir = SCRIPT_DIR / 'log'
    log_dir.mkdir(parents=True, exist_ok=True)
    run_dir = ROOT / STEP / stamp
    fields = ['task_id', 'execution_status', 'output_dir', 'message']
    with (log_dir / f'{STEP}.{stamp}.log').open('x', encoding='utf-8', buffering=1) as log:
        with redirect_stdout(log), redirect_stderr(log):
            success_path = log_dir / f'{STEP}.{stamp}.success.tsv'
            failed_path = log_dir / f'{STEP}.{stamp}.failed.tsv'
            write_rows(success_path, fields, [])
            write_rows(failed_path, fields, [])
            signal.signal(signal.SIGTERM, interrupted)
            signal.signal(signal.SIGINT, interrupted)
            marker = run_dir / 'execution_completion.marker'
            own_run = False
            try:
                run_dir.mkdir(parents=True, exist_ok=False)
                own_run = True
                log_message(f'run_id={stamp}; output_dir={run_dir}')
                build(run_dir)
                with success_path.open('a', encoding='utf-8', newline='') as handle:
                    csv.writer(handle, delimiter='\t', lineterminator='\n').writerow([STEP, 'success', str(run_dir), 'normal_end'])
                with marker.open('x', encoding='utf-8') as handle:
                    handle.write(f'run_id={stamp}\nexecution_status=success\n')
                log_message('task completed / done; execution_status=success; execution_tasks=1 success=1 failed=0')
                log_message('please separately check output completeness, format, counts and biological quality.')
            except BaseException as error:
                if own_run and marker.exists():
                    marker.unlink()

                with success_path.open('w', encoding='utf-8', newline='') as handle:
                    csv.writer(handle, delimiter='\t', lineterminator='\n').writerow(fields)
                with failed_path.open('a', encoding='utf-8', newline='') as handle:
                    csv.writer(handle, delimiter='\t', lineterminator='\n').writerow([STEP, 'failed', str(run_dir), str(error)])
                traceback.print_exc()
                log_message('execution_status=failed; execution_tasks=1 success=0 failed=1; not completed output unavailable in under one step.')
                raise


if __name__ == '__main__':
    main()
