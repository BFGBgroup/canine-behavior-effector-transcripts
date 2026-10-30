#!/usr/bin/python3

import csv
import signal
import traceback
from collections import Counter, defaultdict
from contextlib import ExitStack, redirect_stdout, redirect_stderr
from datetime import datetime
from pathlib import Path

ROOT = Path('/path/to/project/analysis_20260817/5-dogs/7.Allelome.LINK')
STEP = '1.4_annotate_and_merge_dog5_allelome_links'
SCRIPT_DIR = ROOT / 'scripts' / STEP
LINK_RUN = ROOT / '1.3_run_dog5_allelome_link/20260907_205139'
STATUS_ROOT = ROOT / 'scripts/1.3_run_dog5_allelome_link/log'
SUCCESS = STATUS_ROOT / '1.3_run_dog5_allelome_link.20260907_205139.success.tsv'
TRANSCRIPTS = Path('/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.transcript_map.tsv')
META = ['gene_id', 'record_type', 'gene_symbol', 'gene_biotype', 'gene_strand',
        'annotation_source', 'source_gtf_gene_ids', 'variant_count', 'unambiguous_variant_count',
        'multi_gene_variant_count', 'unambiguous_hap1_count_sum', 'unambiguous_hap2_count_sum',
        'multi_gene_hap1_count_sum', 'multi_gene_hap2_count_sum', 'input_file', 'input_line']
NOVEL = ['class_code', 'novelty_category', 'reference_gene_id', 'FEELnc_biotype', 'gene_locus_filter_reason']
RAW = [f'{f}_{side}' for side in ('base', 'target') for f in
       ('chr', 'start', 'end', 'name', 'A1_reads', 'A2_reads', 'total_reads', 'allelic_score', 'allelic_ratio')]
RAW += ['mechanism', 'linkage_score']
EXTRA = ['task_id', 'animal_id', 'tissue', 'source_link_file', 'source_link_line', 'source_metadata_file',
         'endpoint1_gene_id', 'endpoint2_gene_id', 'pair_type', 'overlap_bp', 'interval_gap_bp',
         'strand_relation', 'either_endpoint_multigene', 'both_endpoints_unambiguous_count_positive',
         'identical_haplotype_counts', 'evidence_stratum', 'shared_variant_identity_checked']
FIELDS = RAW + EXTRA + [f'{side}_{f}' for side in ('base', 'target') for f in META + ['endpoint_type'] + NOVEL]


def log(message):
    print(f'[{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} ] {message}', flush=True)


def rows(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def writer(stack, path, fields):
    handle = stack.enter_context(path.open('x', encoding='utf-8', newline=''))
    out = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
    out.writeheader()
    return out


def endpoint_type(meta):

    if meta['record_type'] == 'manual_new_gene':
        return 'candidate_new_gene'
    if meta['record_type'] == 'unannotated_variant':
        return 'unannotated_variant'
    if meta['gene_biotype'] == 'protein_coding':
        return 'known_protein_coding'
    if meta['gene_biotype'] == 'lncRNA':
        return 'known_lncRNA'
    return 'known_other_biotype'


def build(run_dir):
    log(f'[1/3] LINK_RUN={LINK_RUN}; transcript_map={TRANSCRIPTS}; serial=1 threads=1')
    marker = dict(x.split('=', 1) for x in (LINK_RUN / 'execution_completion.marker').read_text().splitlines() if '=' in x)
    if marker.get('execution_status') != 'success':
        raise ValueError('the specified step-three run is incomplete')
    novelty = defaultdict(lambda: defaultdict(set))
    for r in rows(TRANSCRIPTS):
        for f in NOVEL:
            novelty[r['final_gene_id']][f].add(r[f])
    tasks = list(rows(SUCCESS))
    pairs = {}
    summary = Counter()
    with ExitStack() as stack:
        all_rows = writer(stack, run_dir / 'dog5_allelome_links.annotated.tsv', FIELDS)
        separated = {s: writer(stack, run_dir / f'dog5_allelome_links.{s}.tsv', FIELDS) for s in
                     ('overlap_review', 'nonoverlap_multigene_review', 'nonoverlap_no_multigene_flag')}
        task_writer = writer(stack, run_dir / 'dog5_allelome_annotation.task_manifest.tsv',
                             ['task_id', 'source_links', 'source_metadata', 'annotated_rows'])
        for i, task in enumerate(tasks, 1):
            if task['execution_status'] != 'success':
                raise ValueError(f"An unsuccessful task appears in the success list: {task['task_id']}")
            link_file = Path(task['output_path']) / (task['task_id'] + '_links_full_table.txt')
            metadata_file = Path(task['locus_metadata'])

            meta = {}
            for m in rows(metadata_file):
                if m['group_id'] in meta:
                    raise ValueError(f"duplicate metadata group_id: {m['group_id']}")
                meta[m['group_id']] = m
            count = 0
            for line, original in enumerate(rows(link_file), 2):
                a, b = meta[original['name_base']], meta[original['name_target']]
                for side, m in (('base', a), ('target', b)):
                    if (original[f'chr_{side}'] != m['chromosome'] or
                            int(original[f'start_{side}']) != int(m['start0']) or
                            int(original[f'end_{side}']) != int(m['end0'])):
                        raise ValueError(f'Link intervals do not match the metadata: {link_file}:{line}')
                if (a['sample'], a['tissue'], a['chromosome']) != (b['sample'], b['tissue'], b['chromosome']):
                    raise ValueError(f'across sample or chromosome link: {link_file}:{line}')
                overlap = max(0, min(int(a['end0']), int(b['end0'])) - max(int(a['start0']), int(b['start0'])))
                gap = max(0, max(int(a['start0']), int(b['start0'])) - min(int(a['end0']), int(b['end0'])))
                multi = int(a['multi_gene_variant_count']) > 0 or int(b['multi_gene_variant_count']) > 0

                stratum = 'overlap_review' if overlap else ('nonoverlap_multigene_review' if multi else 'nonoverlap_no_multigene_flag')
                first, second = sorted((a['gene_id'], b['gene_id']))
                pair_type = '--'.join(sorted((endpoint_type(a), endpoint_type(b))))
                strand = ('unknown' if a['gene_strand'] not in ('+', '-') or b['gene_strand'] not in ('+', '-')
                          else ('same' if a['gene_strand'] == b['gene_strand'] else 'opposite'))
                out = dict(original)
                out.update(task_id=task['task_id'], animal_id=a['sample'], tissue=a['tissue'],
                           source_link_file=str(link_file), source_link_line=line, source_metadata_file=str(metadata_file),
                           endpoint1_gene_id=first, endpoint2_gene_id=second, pair_type=pair_type,
                           overlap_bp=overlap, interval_gap_bp=gap, strand_relation=strand,
                           either_endpoint_multigene='yes' if multi else 'no',
                           both_endpoints_unambiguous_count_positive='yes' if all(
                               int(m['unambiguous_hap1_count_sum']) + int(m['unambiguous_hap2_count_sum']) > 0 for m in (a, b)) else 'no',
                           identical_haplotype_counts='yes' if (a['hap1_count_sum'], a['hap2_count_sum']) == (b['hap1_count_sum'], b['hap2_count_sum']) else 'no',
                           evidence_stratum=stratum, shared_variant_identity_checked='no')
                for side, m in (('base', a), ('target', b)):
                    out.update({f'{side}_{f}': m[f] for f in META})
                    out[f'{side}_endpoint_type'] = endpoint_type(m)
                    if m['record_type'] == 'manual_new_gene' and m['gene_id'] not in novelty:
                        raise ValueError(f"new gene missing upstream classification: {m['gene_id']}")
                    out.update({f'{side}_{f}': ';'.join(sorted(novelty.get(m['gene_id'], {}).get(f, ()))) or 'NA' for f in NOVEL})
                all_rows.writerow(out)
                separated[stratum].writerow(out)
                key = (a['chromosome'], first, second)
                if key not in pairs:
                    pairs[key] = dict(rows=0, samples=set(), animals=set(), tissues=set(), enhancing=0, repressing=0,
                                      strata=Counter(), max_score=float('-inf'), pair_type=pair_type,
                                      mechanism_samples=defaultdict(set))
                p = pairs[key]
                p['rows'] += 1
                p['samples'].add((a['sample'], a['tissue']))
                p['animals'].add(a['sample']); p['tissues'].add(a['tissue'])
                p[original['mechanism']] += 1
                p['mechanism_samples'][original['mechanism']].add((a['sample'], a['tissue']))
                p['strata'][stratum] += 1
                p['max_score'] = max(p['max_score'], float(original['linkage_score']))
                summary['annotated_rows'] += 1; summary[stratum] += 1
                summary['mechanism_' + original['mechanism']] += 1
                summary['pair_type_' + pair_type] += 1
                count += 1
            task_writer.writerow(dict(task_id=task['task_id'], source_links=str(link_file),
                                      source_metadata=str(metadata_file), annotated_rows=count))
            if i == 1 or i % 20 == 0 or i == len(tasks):
                log(f'[2/3] {i}/{len(tasks)} task already annotation; cumulative {summary["annotated_rows"]} row')
        log('[3/3] summarize undirected endpoint pairs and animals—tissue support, do not merge haplotype count')
        pair_fields = ['chromosome', 'endpoint1_gene_id', 'endpoint2_gene_id', 'pair_type', 'expanded_rows',
                       'sample_count', 'animal_count', 'tissue_count', 'sample_list', 'animal_list', 'tissue_list',
                       'enhancing_rows', 'repressing_rows', 'enhancing_sample_count', 'repressing_sample_count',
                       'mechanism_pattern', 'max_linkage_score', 'overlap_review_rows',
                       'nonoverlap_multigene_review_rows', 'nonoverlap_no_multigene_flag_rows']
        pair_writer = writer(stack, run_dir / 'dog5_allelome_pairs.summary.tsv', pair_fields)
        for (chrom, a, b), p in sorted(pairs.items()):
            pair_writer.writerow(dict(chromosome=chrom, endpoint1_gene_id=a, endpoint2_gene_id=b,
                pair_type=p['pair_type'], expanded_rows=p['rows'], sample_count=len(p['samples']),
                animal_count=len(p['animals']), tissue_count=len(p['tissues']),
                sample_list=';'.join(f'{animal}|{tissue}' for animal, tissue in sorted(p['samples'])),
                animal_list=';'.join(sorted(p['animals'])), tissue_list=';'.join(sorted(p['tissues'])),
                enhancing_rows=p['enhancing'], repressing_rows=p['repressing'],
                enhancing_sample_count=len(p['mechanism_samples']['enhancing']), repressing_sample_count=len(p['mechanism_samples']['repressing']),
                mechanism_pattern='mixed' if p['enhancing'] and p['repressing'] else ('enhancing_only' if p['enhancing'] else 'repressing_only'),
                max_linkage_score=p['max_score'], **{s + '_rows': p['strata'][s] for s in separated}))
        summary['unique_endpoint_pairs'] = len(pairs)
        summary['input_tasks'] = len(tasks)
        sw = writer(stack, run_dir / 'dog5_allelome_annotation.summary.tsv', ['metric', 'value'])
        sw.writerows({'metric': k, 'value': v} for k, v in sorted(summary.items()))


def stop(signum, frame):
    raise InterruptedError(f'received_signal={signum}')


def main():
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_dir = SCRIPT_DIR / 'log'
    log_dir.mkdir(parents=True, exist_ok=True)
    run_dir = ROOT / STEP / stamp
    with ExitStack() as stack:
        master = stack.enter_context((log_dir / f'{STEP}.{stamp}.log').open('x', encoding='utf-8', buffering=1))
        stack.enter_context(redirect_stdout(master)); stack.enter_context(redirect_stderr(master))
        status_fields = ['task_id', 'execution_status', 'output_path', 'reason']
        success_handle = stack.enter_context((log_dir / f'{STEP}.{stamp}.success.tsv').open(
            'x', encoding='utf-8', newline='', buffering=1))
        success = csv.DictWriter(success_handle, fieldnames=status_fields, delimiter='\t', lineterminator='\n')
        success.writeheader()
        success_header_end = success_handle.tell()
        failed = writer(stack, log_dir / f'{STEP}.{stamp}.failed.tsv', status_fields)
        signal.signal(signal.SIGINT, stop); signal.signal(signal.SIGTERM, stop)
        completion = run_dir / 'execution_completion.marker'
        completion_owned = False
        try:
            run_dir.mkdir(parents=True, exist_ok=False)
            log(f'run_id={stamp}; output_dir={run_dir}')
            build(run_dir)

            with completion.open('x', encoding='utf-8') as handle:
                completion_owned = True
                handle.write(f'run_id={stamp}\nexecution_status=success\ninput_run={LINK_RUN}\n')
            success.writerow(dict(task_id=STEP, execution_status='success', output_path=str(run_dir), reason='normal_end'))
            success_handle.flush()
            log('task completed / done; execution_status=success')
            log('please independent check output completeness, annotation and counting scope; candidate link is not causal regulation validate.')
        except BaseException as error:

            if completion_owned:
                completion.unlink(missing_ok=True)
            success_handle.seek(success_header_end)
            success_handle.truncate()
            failed.writerow(dict(task_id=STEP, execution_status='failed', output_path=str(run_dir), reason=str(error)))
            traceback.print_exc(); log('execution_status=failed; not completed output unavailable in downstream analysis.')
            raise


if __name__ == '__main__':
    main()
