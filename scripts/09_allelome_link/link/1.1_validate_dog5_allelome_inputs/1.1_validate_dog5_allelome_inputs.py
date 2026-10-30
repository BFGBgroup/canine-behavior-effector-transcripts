#!/usr/bin/python3


import csv
import re
import signal
import traceback
from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime
from pathlib import Path


ROOT = Path('/path/to/project/analysis_20260817/5-dogs/7.Allelome.LINK')
STEP = '1.1_validate_dog5_allelome_inputs'
SCRIPT_DIR = ROOT / 'scripts' / STEP
OUTPUT_ROOT = ROOT / STEP
TASK_LIST = Path('/path/to/project/legacy/five-dogs/data/statistics/five_dogs_fastq_pair_stats.tsv')
GENE_ROOT = Path('/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation/1.5ase-cnv-clean-gene-id-haplotype-table/20260821_161805')
VARIANT_ROOT = Path('/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation/1.4ase-cnv-qc-filter/20260821_154914')
KNOWN_REFERENCE = Path('/path/to/project/analysis_20260709/5-dogs/3.ASE-process/3.7_build_canfam31_gene_span_reference/CanFam3.1.autosomes.gene_span.annotation.tsv')
NEW_GTF = Path('/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.new_only.gtf')
CHROMS = tuple(f'chr{i}' for i in range(1, 39))
TYPES = ('known_gene', 'manual_new_gene', 'unannotated_variant')
COUNT_UNIT = 'ASEReadCounter_site_allele_support_sum'
AGGREGATION = 'sum_distinct_variant_allele_support_within_gene_id'
META_MAP = {
    'gene_id': 'GeneID', 'gene_symbol': 'GeneSymbol', 'gene_id_base': 'GeneID_base',
    'ncbi_gene_id': 'NCBI_GeneID', 'gene_feature_type': 'GeneFeatureType',
    'gene_biotype': 'GeneBiotype', 'annotation_source': 'AnnotationSource',
    'source_gtf_gene_ids': 'SourceGTFGeneIDs', 'gene_strand': 'GeneStrand',
    'chromosome': 'CHROM',
}
ACTION_COLUMNS = tuple(f'{action}_variant_count' for action in (
    'keep', 'flag_no_cnv_consensus', 'flag_soft_cnv_overlap',
    'flag_strong_cnv_overlap', 'flag_cnv_conflict', 'flag_not_evaluable',
))
COUNT_COLUMNS = (
    'variant_count', 'unambiguous_variant_count', 'multi_gene_variant_count',
    *ACTION_COLUMNS, 'hap1_count_sum', 'hap2_count_sum', 'total_count_sum',
    'unambiguous_hap1_count_sum', 'unambiguous_hap2_count_sum',
    'multi_gene_hap1_count_sum', 'multi_gene_hap2_count_sum',
)
QC_SUM_MAP = {
    'output_groups': 'row_count',
    **{f'{t}_groups': f'{t}_rows' for t in TYPES},
    'output_variant_references': 'variant_count',
    'output_unambiguous_variant_references': 'unambiguous_variant_count',
    'output_multi_gene_variant_references': 'multi_gene_variant_count',
    **{f'output_{c}': c for c in COUNT_COLUMNS if c.endswith('hap1_count_sum') or c.endswith('hap2_count_sum')},
}
QC_CONTRACT = {
    'aggregation_key': 'GeneID+CHROM+SAMPLE+TISSUE',
    'phase_status_required': 'PHASED_HET_EXACT_MATCH',
    'vcf_gt_required': '0|1_or_1|0', 'phase_set_scope': 'chromosome',
    'cnv_clean_pass_required': 'TRUE',
    'cnv_input_policy': 'known_cnv_removed:exclude_only_single_state_consensus_cnv_at_7_of_9_scales',
    'count_unit': COUNT_UNIT, 'aggregation_method': AGGREGATION,
}
STAT_FIELDS = ('row_count', *(f'{t}_rows' for t in TYPES), *COUNT_COLUMNS,
               'count_ge10_rows', 'bias_pass_raw_rows', 'bias_pass_rounded_rows',
               'rounding_changes_pass_rows', 'rows_with_multi_gene_variants')
REFERENCE_FIELDS = ('assembly', 'gene_id', 'record_type', 'gene_symbol', 'gene_id_base',
                    'ncbi_gene_id', 'gene_feature_type', 'gene_biotype', 'annotation_source',
                    'source_gtf_gene_ids', 'gene_strand', 'chromosome', 'start0', 'end0',
                    'start1', 'end1', 'span_bp', 'coordinate_source', 'observed_sample_count')
ISSUE_FIELDS = ('severity', 'sample_key', 'input_path', 'line_number', 'code', 'detail')


def log(message):
    print(f'[{datetime.now():%Y-%m-%d %H:%M:%S} ] {message}', flush=True)


def read_tsv(path, required):

    with path.open('r', encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        names = reader.fieldnames or []
        if len(names) != len(set(names)) or not set(required).issubset(names):
            raise ValueError(f'{path}: duplicate headers or missing fields {sorted(set(required) - set(names))}')
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f'{path}:{reader.line_num}: Column count does not match the header')
            yield reader.line_num, row


def write_tsv(path, fields, rows):
    with path.open('x', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def metrics(path):
    result = {}
    for line, row in read_tsv(path, ('metric', 'value')):
        if row['metric'] in result:
            raise ValueError(f'{path}:{line}: duplicate metric {row["metric"]}')
        result[row['metric']] = row['value']
    return result


def make_reference(meta, kind, start1, end1, source):
    if not (1 <= start1 <= end1):
        raise ValueError(f'{source}: invalid coordinates {meta["gene_id"]}: {start1}-{end1}')
    return dict(meta, assembly='CanFam3.1', record_type=kind, start0=start1 - 1,
                end0=end1, start1=start1, end1=end1, span_bp=end1 - start1 + 1,
                coordinate_source=str(source))


def load_references():

    references = {}
    for line, row in read_tsv(KNOWN_REFERENCE, (*META_MAP.values(), 'START1', 'END1', 'START0', 'END0')):
        meta = {key: row[value] for key, value in META_MAP.items()}
        gene = meta['gene_id']
        if gene in references or meta['chromosome'] not in CHROMS or meta['gene_strand'] not in ('+', '-'):
            raise ValueError(f'{KNOWN_REFERENCE}:{line}: duplicate ID or invalid chromosome/strand')
        rec = make_reference(meta, 'known_gene', int(row['START1']), int(row['END1']), KNOWN_REFERENCE)
        if (rec['start0'], rec['end0']) != (int(row['START0']), int(row['END0'])):
            raise ValueError(f'{KNOWN_REFERENCE}:{line}: 0-based/1-based coordinates mismatch')
        references[gene] = rec
    manual = {}
    with NEW_GTF.open('r', encoding='utf-8') as handle:
        for line, text in enumerate(handle, 1):
            if text.startswith('#') or not text.strip():
                continue
            cols = text.rstrip('\n').split('\t')
            if len(cols) != 9 or cols[2] not in ('transcript', 'exon'):
                raise ValueError(f'{NEW_GTF}:{line}: is not checked of transcript/exon GTF row')
            attrs = dict(re.findall(r'(\S+)\s+"([^"]*)"', cols[8]))
            required_attrs = ('gene_id', 'gene_name', 'transcript_id', 'source_dataset', 'original_gene_id', 'selection_stage')
            if any(not attrs.get(key) for key in required_attrs):
                raise ValueError(f'{NEW_GTF}:{line}: required GTF attribute missing')
            gene = attrs['gene_id']
            if (not gene.startswith('DOG5_NEWGENE_') or attrs['gene_name'] != gene
                    or attrs['source_dataset'] != 'five-dogs'
                    or attrs['selection_stage'] != '2.1_transcript_first_preASE_candidate_assessment'):
                raise ValueError(f'{NEW_GTF}:{line}: is not five-dog-specific candidate annotation')
            if gene in references or cols[6] not in ('+', '-'):
                raise ValueError(f'{NEW_GTF}:{line}: ID conflict or strand information unexpected')
            start, end = int(cols[3]), int(cols[4])
            if not 1 <= start <= end:
                raise ValueError(f'{NEW_GTF}:{line}: GTF coordinates unexpected')
            rec = manual.setdefault(gene, dict(chrom=cols[0], strand=cols[6], name=attrs['gene_name'],
                                               start=start, end=end, source_ids=set(), features=set()))
            if (rec['chrom'], rec['strand'], rec['name']) != (cols[0], cols[6], attrs['gene_name']):
                raise ValueError(f'{NEW_GTF}:{line}: the same candidate of chromosome, strand or name conflict')
            rec['start'], rec['end'] = min(rec['start'], start), max(rec['end'], end)
            rec['source_ids'].add(f'{attrs["source_dataset"]}:{attrs["original_gene_id"]}')
            rec['features'].add(cols[2])
    for gene, rec in manual.items():
        if rec['features'] != {'transcript', 'exon'}:
            raise ValueError(f'{NEW_GTF}: {gene} missing transcript or exon')
        if rec['chrom'] not in CHROMS:
            continue
        meta = dict(gene_id=gene, gene_symbol=rec['name'], gene_id_base=gene, ncbi_gene_id='NA',
                    gene_feature_type='manual_new_gene', gene_biotype='manual_new_gene',
                    annotation_source='ManualNewGene', source_gtf_gene_ids='|'.join(sorted(rec['source_ids'])),
                    gene_strand=rec['strand'], chromosome=rec['chrom'])
        references[gene] = make_reference(meta, 'manual_new_gene', rec['start'], rec['end'], NEW_GTF)
    return references, len(manual), sum(rec['chrom'] not in CHROMS for rec in manual.values())


def audit(run_dir, issues):
    error_counts = Counter()
    inventory = []

    def issue(sample, path, line, code, detail):
        error_counts[sample] += 1
        issues.writerow(dict(severity='ERROR', sample_key=sample, input_path=str(path),
                             line_number=line, code=code, detail=detail))

    def record_input(role, path, sample='ALL'):
        inventory.append(dict(role=role, sample_key=sample, input_path=str(path), exists=path.is_file()))
        if not path.is_file():
            issue(sample, path, 0, 'missing_input', role)

    log('[stage 1/4] check sample manifest and annotation source')
    for role, path in [('task_list', TASK_LIST), ('known_gene_reference', KNOWN_REFERENCE), ('five_only_new_gtf', NEW_GTF)]:
        record_input(role, path)
    references, manual_all, manual_nonauto = load_references()
    reference_counts = Counter(rec['record_type'] for rec in references.values())
    samples, seen_samples = [], set()
    for line, row in read_tsv(TASK_LIST, ('animal_id', 'term', 'raw_data_name')):
        animal, tissue = row['animal_id'], row['term']
        if not all(re.fullmatch(r'[A-Za-z0-9_-]+', value) for value in (animal, tissue)):
            raise ValueError(f'{TASK_LIST}:{line}: sample identifier cannot used for safe path build')
        if (animal, tissue) in seen_samples:
            raise ValueError(f'{TASK_LIST}:{line}: animal_id + term duplicate')
        seen_samples.add((animal, tissue))
        samples.append((animal, tissue))
    if not samples:
        raise ValueError(f'{TASK_LIST}: no sample')
    log(f'sample={len(samples)}; known gene={reference_counts["known_gene"]}; autosome candidate={reference_counts["manual_new_gene"]}; non-autosomal candidate={manual_nonauto}')

    summaries, task_rows, observed = [], [], Counter()
    log('[stage 2/4] per sample read 1.5 main table, check QC, count with coordinates; no repeated filtering')
    for index, (animal, tissue) in enumerate(samples, 1):
        sample = f'{animal}_{tissue}'
        stem = GENE_ROOT / animal / tissue / sample
        main = Path(str(stem) + '.cnv_clean.gene_id_haplotype.tsv')
        qc_path = Path(str(stem) + '.cnv_clean.gene_id_haplotype.qc.tsv')
        variant = VARIANT_ROOT / animal / tissue / f'{sample}.ase_haplotype.cnv_clean.tsv'
        for role, path in [('gene_haplotype', main), ('gene_haplotype_qc', qc_path), ('upstream_cnv_clean_exists_only', variant)]:
            record_input(role, path, sample)
        qc = metrics(qc_path)
        expected_qc = dict(QC_CONTRACT, animal_id=animal, tissue=tissue, input_tsv=str(variant),
                           output_tsv=str(main), qc_tsv=str(qc_path))
        for key, expected in expected_qc.items():
            if qc.get(key) != expected:
                issue(sample, qc_path, 0, 'qc_contract_mismatch', f'{key}: observed={qc.get(key)!r}; expected={expected!r}')
        counts, per_chrom, seen = Counter(), {chrom: Counter() for chrom in CHROMS}, set()
        required = (*META_MAP, 'group_id', 'record_type', 'sample', 'tissue', *COUNT_COLUMNS, 'count_unit', 'aggregation_method')
        for line, row in read_tsv(main, required):
            gene, chrom, kind = row['gene_id'], row['chromosome'], row['record_type']
            if kind not in TYPES or chrom not in CHROMS:
                raise ValueError(f'{main}:{line}: unsupported of record_type or chromosome')
            checks = {
                'sample_identity': row['sample'] == animal and row['tissue'] == tissue,
                'group_id': row['group_id'] == f'{gene}|{chrom}|{animal}|{tissue}',
                'unique_group': row['group_id'] not in seen,
                'count_unit': row['count_unit'] == COUNT_UNIT,
                'aggregation_method': row['aggregation_method'] == AGGREGATION,
            }
            if row['group_id'] not in seen:
                observed[gene] += 1
            seen.add(row['group_id'])
            vals = {key: int(row[key]) for key in COUNT_COLUMNS}
            checks['nonnegative_counts'] = all(value >= 0 for value in vals.values())
            checks['nonempty_group'] = vals['variant_count'] > 0 and vals['total_count_sum'] > 0
            checks['hap_total'] = vals['hap1_count_sum'] + vals['hap2_count_sum'] == vals['total_count_sum']
            checks['variant_partition'] = vals['unambiguous_variant_count'] + vals['multi_gene_variant_count'] == vals['variant_count']
            checks['cnv_action_partition'] = sum(vals[key] for key in ACTION_COLUMNS) == vals['variant_count']
            for hap in ('hap1', 'hap2'):
                checks[f'{hap}_partition'] = vals[f'unambiguous_{hap}_count_sum'] + vals[f'multi_gene_{hap}_count_sum'] == vals[f'{hap}_count_sum']
            if kind == 'unannotated_variant':
                match = re.fullmatch(r'\(NA\)(chr\d+):(\d+)_([ACGT])>([ACGT])', gene)
                if not match or match[1] != chrom or int(match[2]) < 1 or match[3] == match[4]:
                    raise ValueError(f'{main}:{line}: unannotated SNP ID format unexpected: {gene}')
                meta = dict(gene_id=gene, gene_symbol=gene, gene_id_base=gene, ncbi_gene_id='NA',
                            gene_feature_type='unannotated', gene_biotype='unannotated',
                            annotation_source='unannotated', source_gtf_gene_ids='NA', gene_strand='NA', chromosome=chrom)
                references.setdefault(gene, make_reference(meta, kind, int(match[2]), int(match[2]), GENE_ROOT))
                checks['unannotated_single_variant'] = vals['variant_count'] == 1
            ref = references.get(gene)
            checks['reference_found'] = ref is not None
            if ref is not None:
                for key in (*META_MAP, 'record_type'):

                    checks[f'reference_{key}'] = (set(row[key].split('|')) == set(ref[key].split('|'))
                                                 if key == 'source_gtf_gene_ids' else row[key] == ref[key])
            for code, passed in checks.items():
                if not passed:
                    detail = f'gene_id={gene}; group_id={row["group_id"]}'
                    if ref is not None and code.startswith('reference_') and code != 'reference_found':
                        field = code[len('reference_'):]
                        detail += f'; observed={row[field]!r}; expected={ref[field]!r}'
                    issue(sample, main, line, code, detail)

            total = vals['total_count_sum']
            ratio = vals['hap1_count_sum'] / total if total > 0 else 0.5
            raw_pass = total >= 10 and (ratio >= 0.6 or ratio <= 0.4)
            rounded = round(ratio, 2)
            rounded_pass = total >= 10 and (rounded >= 0.6 or rounded <= 0.4)
            delta = dict(vals, row_count=1, **{f'{kind}_rows': 1}, count_ge10_rows=int(total >= 10),
                         bias_pass_raw_rows=int(raw_pass), bias_pass_rounded_rows=int(rounded_pass),
                         rounding_changes_pass_rows=int(raw_pass != rounded_pass),
                         rows_with_multi_gene_variants=int(vals['multi_gene_variant_count'] > 0))
            counts.update(delta)
            per_chrom[chrom].update(delta)
        for key, field in QC_SUM_MAP.items():
            if qc.get(key) != str(counts[field]):
                issue(sample, qc_path, 0, 'qc_count_mismatch', f'{key}: QC={qc.get(key)!r}; table={counts[field]}')
        summaries.append(dict(sample_key=sample, animal_id=animal, tissue=tissue, gene_haplotype_tsv=str(main),
                              gene_haplotype_qc_tsv=str(qc_path), upstream_cnv_clean_tsv=str(variant),
                              **{key: counts[key] for key in STAT_FIELDS}))
        for chrom in CHROMS:
            task_rows.append(dict(task_id=f'{sample}_{chrom}', sample_key=sample, animal_id=animal,
                                  tissue=tissue, chromosome=chrom, gene_haplotype_tsv=str(main),
                                  **{key: per_chrom[chrom][key] for key in STAT_FIELDS}))
        log(f'[{index}/{len(samples)}] {sample}: rows={counts["row_count"]}; input_errors={error_counts[sample]}')

    log('[stage 3/4] output standardized coordinates reference and task list; do not start downstream tasks')
    for row in summaries:
        row['audit_status'] = 'FAIL' if error_counts[row['sample_key']] or error_counts['ALL'] else 'PASS'
        row['error_count'] = error_counts[row['sample_key']]
    for row in task_rows:
        row['audit_status'] = 'FAIL' if error_counts[row['sample_key']] or error_counts['ALL'] else 'PASS'

        row['preparation_status'] = 'BLOCKED' if row['audit_status'] == 'FAIL' else ('EMPTY' if row['row_count'] == 0 else 'READY')
    write_tsv(run_dir / 'dog5_allelome_input_audit.tsv', ('role', 'sample_key', 'input_path', 'exists'), inventory)
    reference_rows = (dict(rec, observed_sample_count=observed[gene]) for gene, rec in
                      sorted(references.items(), key=lambda item: (int(item[1]['chromosome'][3:]), item[1]['start0'], item[0])))
    write_tsv(run_dir / 'dog5_allelome_locus_reference.tsv', REFERENCE_FIELDS, reference_rows)
    write_tsv(run_dir / 'dog5_allelome_sample_manifest.tsv', (*summaries[0],), summaries)
    write_tsv(run_dir / 'dog5_allelome_chromosome_task_manifest.tsv', (*task_rows[0],), task_rows)
    status = 'FAIL' if sum(error_counts.values()) else 'PASS'
    summary = dict(audit_status=status, sample_count=len(samples), animal_count=len({a for a, t in samples}),
                   chromosome_task_count=len(task_rows), input_gene_sample_rows=sum(r['row_count'] for r in summaries),
                   known_reference_genes=reference_counts['known_gene'], candidate_genes_all_contigs=manual_all,
                   candidate_genes_autosomes=reference_counts['manual_new_gene'], candidate_genes_non_autosomes=manual_nonauto,
                   observed_unique_loci=len(observed), error_count=sum(error_counts.values()))
    write_tsv(run_dir / 'dog5_allelome_validation_summary.tsv', ('metric', 'value'),
              (dict(metric=key, value=value) for key, value in summary.items()))
    contract = dict(assembly='CanFam3.1', chromosome_scope='chr1-chr38', key='GeneID+CHROM+SAMPLE+TISSUE',
                    count_unit=COUNT_UNIT, aggregation_method=AGGREGATION, phase_set_scope='chromosome',
                    coordinate_0based='start0,end0: 0-based half-open',
                    coordinate_1based='start1,end1: 1-based inclusive',
                    unannotated_snv_interval='start0=POS-1;end0=POS;start1=end1=POS',
                    gene_interval='known: reference START1/END1; manual: min/max transcript+exon GTF',
                    haplotype_orientation='keep HAP1/HAP2; do not flip on gene strand',
                    multi_gene_policy='retain input assignments; gene pairs may share counted SNPs',
                    preview_only='no row filtering; total>=10 AND (HAP1/total<=0.4 OR HAP1/total>=0.6); report raw and rounded(2)',
                    next_step='prepare LINK inputs separately; this reference is NOT a LINK input table',
                    phase_limit='QC contract checked; no BCF revalidation or independent phasing assessment',
                    source_variant_limit='1.4 path and existence checked only; no variant-level recount',
                    execution_vs_audit='execution success means audit completed; proceed only after audit_status=PASS')
    write_tsv(run_dir / 'dog5_allelome_parameter_contract.tsv', ('parameter', 'value'),
              (dict(parameter=key, value=value) for key, value in contract.items()))
    log(f'[stage 4/4] input audit completed: audit_status={status}; errors={sum(error_counts.values())}')
    return status


def interrupted(signum, frame):
    raise InterruptedError(f'task received termination signal {signum}')


def main():
    run_id = datetime.now().strftime('%Y%m%d_%H%M%S')
    run_dir = OUTPUT_ROOT / run_id
    log_dir = SCRIPT_DIR / 'log'
    log_dir.mkdir(parents=True, exist_ok=True)
    master_log = log_dir / f'{STEP}.{run_id}.log'
    fields = ('run_id', 'task_id', 'status', 'exit_code', 'started_at', 'finished_at', 'output_path', 'task_log', 'reason')

    with master_log.open('x', encoding='utf-8', buffering=1) as handle, redirect_stdout(handle), redirect_stderr(handle):
        started = datetime.now().isoformat(timespec='seconds')
        signal.signal(signal.SIGTERM, interrupted)
        signal.signal(signal.SIGINT, interrupted)
        with (log_dir / f'{STEP}.{run_id}.success.tsv').open('x', encoding='utf-8', newline='') as success, \
             (log_dir / f'{STEP}.{run_id}.failed.tsv').open('x', encoding='utf-8', newline='') as failed:
            writers = [csv.DictWriter(h, fieldnames=fields, delimiter='\t', lineterminator='\n') for h in (success, failed)]
            for writer in writers:
                writer.writeheader()
            success.flush()
            failed.flush()
            state = dict(run_id=run_id, task_id=STEP, started_at=started, output_path=str(run_dir), task_log=str(master_log))
            marker_path = run_dir / 'execution_completion.marker'
            marker_owned = False
            try:
                log(f'start; concurrent_tasks=1; threads_per_task=1(new script default, pure standard library serial); output={run_dir}')
                log(f'input_1.5={GENE_ROOT}; known_reference={KNOWN_REFERENCE}; candidate_gtf={NEW_GTF}')
                run_dir.mkdir(parents=True, exist_ok=False)
                with (run_dir / 'dog5_allelome_validation_issues.tsv').open('x', encoding='utf-8', newline='') as issue_file:
                    issues = csv.DictWriter(issue_file, fieldnames=ISSUE_FIELDS, delimiter='\t', lineterminator='\n')
                    issues.writeheader()
                    audit_status = audit(run_dir, issues)

                writers[0].writerow(dict(state, status='success', exit_code=0, finished_at=datetime.now().isoformat(timespec='seconds'), reason=f'audit_completed;audit_status={audit_status}'))
                success.flush()
                with marker_path.open('x', encoding='utf-8') as marker:
                    marker_owned = True
                    marker.write(f'run_id={run_id}\nexecution_status=success\naudit_status={audit_status}\ncompleted_at={datetime.now().isoformat(timespec="seconds")}\n')
                log('task completed / done; execution_status=success; total_tasks=1; success=1; failed=0')
                log('please separately check output completeness, format and biological interpretation; only audit_status=PASS and review confirm after only then enter under one step.')
            except BaseException as exc:
                traceback.print_exc()


                if marker_owned:
                    marker_path.unlink()
                success.seek(0)
                success.truncate()
                writers[0].writeheader()
                success.flush()
                writers[1].writerow(dict(state, status='failed', exit_code=1, finished_at=datetime.now().isoformat(timespec='seconds'), reason=str(exc)))
                failed.flush()
                log('execution_status=failed; total_tasks=1; success=0; failed=1; this run did not finish normally, do not retain the completion marker.')
                return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
