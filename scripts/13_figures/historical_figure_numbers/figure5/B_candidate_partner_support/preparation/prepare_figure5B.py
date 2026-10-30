#!/usr/bin/python3


import csv
import gzip
import os
import re
import sys
import traceback
from collections import defaultdict
from contextlib import redirect_stdout, redirect_stderr
from datetime import datetime
from pathlib import Path

ROOT = Path('/path/to/project/analysis_20260817')
OUT_ROOT = ROOT / 'article/figure/figure5/B_candidate_partner_support/preparation'
CANDIDATE_INPUT = ROOT / 'results/5-dogs/8.new-gene-allelome-eqtl-gwas-finemap-integration/1.1_integrate_dog5_newgene_allelome_eqtl_gwas_finemap/20260910_195332/dog5_newgene_eqtl_gwas_finemap_all_candidates.tsv'
LOCUS_INPUT = Path('/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.10_run_dog_behavior_susie_rss_finemapping/20260802_205352/locus_run_manifest.tsv')
M5 = ROOT / 'results/5-dogs/5.allelome-compare-gwas-finemap/summary/20260910_182124/5dogs_gwas_finemap_allelome_main_focus.tsv'
M8 = ROOT / 'results/5-dogs/8.new-gene-allelome-eqtl-gwas-finemap-integration/1.1_integrate_dog5_newgene_allelome_eqtl_gwas_finemap/20260910_195332/dog5_integrated_newgene_pcGene_links.tsv.gz'
LINK_ROOT = ROOT / '5-dogs/7.Allelome.LINK'
MANIFEST = LINK_ROOT / '1.2_build_dog5_allelome_locus_tables/20260907_202637/dog5_allelome_locus_table_manifest.tsv'
ANNOTATED = LINK_ROOT / '1.4_annotate_and_merge_dog5_allelome_links/20260907_213429/dog5_allelome_links.annotated.tsv'
ANNOTATION_TASKS = LINK_ROOT / '1.4_annotate_and_merge_dog5_allelome_links/20260907_213429/dog5_allelome_annotation.task_manifest.tsv'
SUCCESS = LINK_ROOT / 'scripts/1.3_run_dog5_allelome_link/log/1.3_run_dog5_allelome_link.20260907_205139.success.tsv'
INPUTS = dict(candidate_input=CANDIDATE_INPUT, locus_manifest=LOCUS_INPUT,
              module5=M5, module8=M8, endpoint_manifest=MANIFEST,
              annotated_links=ANNOTATED, annotation_tasks=ANNOTATION_TASKS, success_list=SUCCESS)
DOGS = ['Ace_396226', 'CodRou_346627', 'Crak_358336', 'Dozer_386519', 'SanMcg_349891']
TISSUES = ['cerebellum', 'cortex_frontal', 'cortex_occipital', 'pituitary_gland']
MIN_SUPPORT = 10
WINDOW_BP = 500000
LEAD_P_THRESHOLD = 1e-5


def read_tsv(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def write_tsv(path, rows, fields):
    with open(path, 'x', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: 'NA' if row.get(k) is None else row.get(k, 'NA') for k in fields})


def collapse(values):
    return ';'.join(sorted({str(v) for v in values if v not in (None, '', 'NA')})) or 'NA'


def log(message):
    print(datetime.now().strftime('[%Y-%m-%d %H:%M:%S ]'), message, flush=True)


def eligibility(candidate, partner):

    if len(candidate) > 1 or len(partner) > 1:
        return None, 'ambiguous_endpoint'
    if not candidate or not partner:
        return False, 'both_endpoints_missing' if not candidate and not partner else 'one_endpoint_missing'
    c, p = candidate[0], partner[0]
    if c['record_type'] != 'manual_new_gene' or p['gene_biotype'] != 'protein_coding':
        return False, 'endpoint_type_mismatch'
    if p['record_type'] in ('manual_new_gene', 'unannotated_variant'):
        return False, 'endpoint_type_mismatch'
    if float(c['total_reads']) < MIN_SUPPORT or float(p['total_reads']) < MIN_SUPPORT:
        return False, 'support_below_10'
    if c['chr'] != p['chr'] or c['chr'] not in {f'chr{i}' for i in range(1, 39)}:
        return False, 'not_same_autosome'
    if int(p['end']) < int(c['start']) - WINDOW_BP or int(p['start']) > int(c['end']) + WINDOW_BP:
        return False, 'outside_window'
    return True, 'eligible'


def select_candidates():

    molecular = []
    tiers = ['A_BOTH_95PCT_CREDIBLE_PLUS_ALLELOME',
             'B_GWAS_NONCREDIBLE_PLUS_ALLELOME',
             'C_BOTH_95PCT_CREDIBLE_NO_ALLELOME',
             'D_GWAS_NONCREDIBLE_NO_ALLELOME']
    seen = set()
    for number, r in enumerate(read_tsv(CANDIDATE_INPUT), 1):
        if not r['candidate_id'] or r['candidate_id'] in seen:
            raise ValueError('No. 8 module candidate_id missing or duplicate')
        seen.add(r['candidate_id'])

        if (r['gwas_branch'] != 'main' or r['eqtl_is_credible_set_member'] != '1'
                or r['eqtl_priority_ppa_ge_0.1'] != '1'
                or r['eqtl_tissue'] not in TISSUES or r['evidence_tier'] not in tiers):
            raise ValueError('No. 8 module branch, tissue or molecular CS/PPA input contract change')
        lead = float(r['gwas_locus_lead_p_lrt'])
        ppa = float(r['eqtl_ppa'])
        if not (re.fullmatch(r'DOG5_NEWGENE_[0-9]+', r['new_gene_id'])
                and 0 < lead < LEAD_P_THRESHOLD and 0.1 <= ppa <= 1
                and r['gwas_overlap_snp_in_credible_set'] == 'yes'
                and r['allelome_block_match'] == 'yes'):
            continue
        molecular.append((r, number))
    molecular.sort(key=lambda z: (tiers.index(z[0]['evidence_tier']), z[0]['new_gene_id'],
        TISSUES.index(z[0]['eqtl_tissue']), z[0]['gwas_trait'], z[0]['gwas_locus_id'],
        z[0]['gwas_overlap_snp_id'], z[0]['candidate_id']))
    candidates = []
    for r, number in molecular:
        candidates.append(dict(candidate_id=r['candidate_id'], new_gene_id=r['new_gene_id'],
            eqtl_tissue=r['eqtl_tissue'], gwas_trait=r['gwas_trait'],
            gwas_locus_id=r['gwas_locus_id'], gwas_overlap_snp_id=r['gwas_overlap_snp_id'],
            gwas_branch=r['gwas_branch'], gwas_locus_lead_p_lrt=r['gwas_locus_lead_p_lrt'],
            evidence_route='Molecular', inclusion_status='Primary',
            selection_reason='regional_lead_P_lt_1e-5_and_both_CS_and_LINK',
            source_file=str(CANDIDATE_INPUT), source_record_numbers=str(number),
            regional_lead_source=str(CANDIDATE_INPUT)))
    loci = {}
    for r in read_tsv(LOCUS_INPUT):
        if r['branch'] == 'main':
            if r['locus_id'] in loci:
                raise ValueError('Primary-region manifest contains a duplicate region key')
            loci[r['locus_id']] = r
    positions = defaultdict(list)
    for number, r in enumerate(read_tsv(M5), 1):
        if (r['gwas_branch'] == 'main'
                and re.fullmatch(r'DOG5_NEWGENE_[0-9]+', r['manual_new_gene'])
                and r['gwas_credible_sets'] not in ('', 'NA')):
            snp = r['gwas_variant_display'].split(' ')[0]
            key = '|'.join((r['manual_new_gene'], r['gwas_trait'], snp, r['gwas_locus_id']))
            positions[key].append((number, r))

    for key in sorted(positions):
        rows = positions[key]
        r = rows[0][1]
        lead = loci[r['gwas_locus_id']]['gwas_lead_p_lrt']
        primary = float(lead) < LEAD_P_THRESHOLD
        candidates.append(dict(candidate_id='module5::' + key, new_gene_id=r['manual_new_gene'],
            eqtl_tissue=None, gwas_trait=r['gwas_trait'], gwas_locus_id=r['gwas_locus_id'],
            gwas_overlap_snp_id=r['gwas_variant_display'].split(' ')[0],
            gwas_branch='main', gwas_locus_lead_p_lrt=lead,
            evidence_route='Positional', inclusion_status='Primary' if primary else 'Exploratory',
            selection_reason='position_route_regional_lead_P_lt_1e-5' if primary
                else 'user_retained_original_focus_exploratory_relation',
            source_file=str(M5), source_record_numbers=';'.join(str(n) for n, _ in rows),
            regional_lead_source=str(LOCUS_INPUT)))
    for index, r in enumerate(candidates, 1):
        r['plot_row'] = index
    return candidates


def prepare(out):
    log('[1/5] from original integrate result reconstruct candidate and read link; concurrent_tasks=1 threads_per_task=1')
    for role, path in INPUTS.items():
        log(f'{role}={path}')
    candidates = select_candidates()
    candidate_by_id = {a['candidate_id']: a for a in candidates}
    molecular = {a['candidate_id']: a for a in candidates if a['evidence_route'] == 'Molecular'}
    positional = {(a['new_gene_id'], a['gwas_trait'], a['gwas_locus_id'],
                   a['gwas_overlap_snp_id']): a for a in candidates if a['evidence_route'] == 'Positional'}
    selected = []

    def add(a, partner, symbol, chrom, dog, tissue, mechanism, score, stratum, path, row_number):
        selected.append(dict(candidate_id=a['candidate_id'], new_gene_id=a['new_gene_id'],
            evidence_route=a['evidence_route'], inclusion_status=a['inclusion_status'],
            pair_id=a['new_gene_id'] + '::' + partner, partner_gene_id=partner,
            partner_symbol=symbol, chromosome=chrom, animal_id=dog, tissue=tissue,
            mechanism=mechanism, linkage_score=score, evidence_stratum=stratum,
            source_file=str(path), source_record_number=row_number))

    for number, r in enumerate(read_tsv(M8), 1):
        if r['candidate_id'] in molecular:
            a = molecular[r['candidate_id']]
            add(a, r['partner_gene_id'], r['partner_gene_symbol'], r['chrom'], r['animal_id'],
                r['link_tissue'], r['mechanism'], r['linkage_score'], r['evidence_stratum'], M8, number)
    for number, r in enumerate(read_tsv(M5), 1):
        key = (r['manual_new_gene'], r['gwas_trait'], r['gwas_locus_id'], r['gwas_variant_display'].split(' ')[0])
        a = positional.get(key)
        if a is not None and r['gwas_branch'] == a['gwas_branch']:
            dog, tissue = r['sample_tissue'].split('|')
            add(a, r['protein_coding_gene_id'], r['protein_coding_gene'], r['new_gene_span'].split(':')[0],
                dog, tissue, r['mechanism'], r['linkage_score'], r['evidence_stratum'], M5, number)
    pairs = {}
    selected_by_pair = defaultdict(list)
    selected_cells = defaultdict(set)
    for r in selected:
        pair = r['pair_id']
        selected_by_pair[pair].append(r)
        selected_cells[(pair, r['animal_id'], r['tissue'])].add(r['candidate_id'])
        pairs[pair] = {k: r[k] for k in ('pair_id', 'new_gene_id', 'partner_gene_id', 'partner_symbol', 'chromosome')}
    log(f'candidate record={len(candidates)}; existing expanded records={len(selected)}; candidate—partner pair={len(pairs)}')

    log('[2/5] read these pair of all tissue LINK background, retain original source with assignment marker')
    observations = defaultdict(list)
    raw_rows = []
    for number, r in enumerate(read_tsv(ANNOTATED), 1):
        if r['pair_type'] != 'candidate_new_gene--known_protein_coding':
            continue
        c, p = ('base', 'target') if r['base_endpoint_type'] == 'candidate_new_gene' else ('target', 'base')
        pair = r[c + '_gene_id'] + '::' + r[p + '_gene_id']
        if pair not in pairs:
            continue
        item = dict(pair_id=pair, animal_id=r['animal_id'], tissue=r['tissue'], task_id=r['task_id'],
            candidate_endpoint=r['name_' + c], partner_endpoint=r['name_' + p],
            mechanism=r['mechanism'], linkage_score=r['linkage_score'], evidence_stratum=r['evidence_stratum'],
            either_endpoint_multigene=r['either_endpoint_multigene'], overlap_bp=r['overlap_bp'],
            source_link_file=r['source_link_file'], source_link_line=r['source_link_line'],
            annotated_source_record_number=number)
        observations[(pair, r['animal_id'], r['tissue'])].append(item)
        raw_rows.append(item)

    log('[3/5] by endpoint and upstream execution record reconstruct per sample status')
    tasks = list(read_tsv(MANIFEST))
    task_by_key = {(r['animal_id'], r['tissue'], r['chromosome']): r for r in tasks}
    sample_set = {(r['animal_id'], r['tissue']) for r in tasks}
    success = {r['task_id']: r for r in read_tsv(SUCCESS)}
    annotated_tasks = {r['task_id']: r for r in read_tsv(ANNOTATION_TASKS)}
    genes = {p['new_gene_id'] for p in pairs.values()} | {p['partner_gene_id'] for p in pairs.values()}
    chromosomes = sorted({p['chromosome'] for p in pairs.values()})
    used_tasks, cells = [], []
    endpoint_fields = ['name', 'start', 'end', 'total_reads', 'allelic_ratio', 'variant_count',
                       'multi_gene_variant_count', 'count_unit', 'record_type', 'gene_biotype']
    for dog in DOGS:
        for tissue in TISSUES:
            endpoints = defaultdict(list)
            ready = {}
            for chrom in chromosomes:
                task = task_by_key.get((dog, tissue, chrom))
                if task is None:
                    ready[chrom] = False
                    continue
                used_tasks.append(task)
                s = success.get(task['task_id'], {})
                a = annotated_tasks.get(task['task_id'], {})
                ready[chrom] = (s.get('execution_status') == 'success' and s.get('exit_code') == '0'
                    and s.get('input_file') == task['locus_table']
                    and s.get('locus_metadata') == task['locus_metadata']
                    and a.get('source_metadata') == task['locus_metadata'])
                metadata = {m['group_id']: m for m in read_tsv(task['locus_metadata'])}
                for ep in read_tsv(task['locus_table']):
                    meta = metadata[ep['name']]
                    if meta['gene_id'] in genes:
                        combined = dict(meta)
                        combined.update(ep)
                        endpoints[(chrom, meta['gene_id'])].append(combined)
            for pair, catalog in pairs.items():
                chrom = catalog['chromosome']
                task_present = (dog, tissue, chrom) in task_by_key
                c = endpoints[(chrom, catalog['new_gene_id'])]
                p = endpoints[(chrom, catalog['partner_gene_id'])]
                eligible, reason = eligibility(c, p) if task_present else (None, 'task_missing')
                if (dog, tissue) not in sample_set:
                    reason = 'sample_missing'
                observed = observations[(pair, dog, tissue)]

                if observed:
                    state = 'detected'
                elif not ready.get(chrom, False) or eligible is None:
                    state = 'unknown'
                elif eligible:
                    state = 'not_detected'
                else:
                    state = 'not_evaluable'
                records = selected_by_pair[pair]
                scope = any(candidate_by_id[r['candidate_id']]['evidence_route'] == 'Positional'
                    or candidate_by_id[r['candidate_id']]['eqtl_tissue'] == tissue for r in records)
                selected_ids = selected_cells[(pair, dog, tissue)]
                identity_ok = (len(c) == len(p) == 1 and all(
                    r['candidate_endpoint'] == c[0]['name'] and r['partner_endpoint'] == p[0]['name'] for r in observed))
                cell = dict(catalog, animal_id=dog, tissue=tissue,
                    sample_present=(dog, tissue) in sample_set, task_present=task_present,
                    upstream_call_complete=ready.get(chrom, False), pair_eligible=eligible,
                    eligibility_reason=reason, pair_state=state, observed_link_rows=len(observed),
                    reported_mechanisms=collapse(r['mechanism'] for r in observed),
                    linkage_scores=collapse(r['linkage_score'] for r in observed),
                    evidence_strata=collapse(r['evidence_stratum'] for r in observed),
                    has_unflagged_link=any(r['evidence_stratum'] == 'nonoverlap_no_multigene_flag' for r in observed),
                    observed_endpoint_match=identity_ok if observed else None,
                    detected_and_eligible=(eligible is True) if observed else None,
                    within_A_tissue_scope=scope, A_selected_link=bool(selected_ids),
                    A_candidate_ids=collapse(selected_ids),
                    background_only_link=bool(observed) and not bool(selected_ids),
                    phase_comparability='not_assessed')
                for prefix, eps in [('candidate', c), ('partner', p)]:
                    for field in endpoint_fields:
                        cell[prefix + '_' + field] = eps[0][field] if len(eps) == 1 else None
                cells.append(cell)
            log(f'{dog} / {tissue}; sample exists={(dog, tissue) in sample_set}')

    log('[4/5] summary independent dog, tissue, shared partner; retain all pair, do not select main-figure rows')
    summaries = []
    for pair, catalog in pairs.items():
        pcells = [c for c in cells if c['pair_id'] == pair]
        detected = [c for c in pcells if c['pair_state'] == 'detected']
        original = selected_by_pair[pair]
        summaries.append(dict(catalog,
            A_candidate_ids=collapse(r['candidate_id'] for r in original),
            evidence_routes=collapse(r['evidence_route'] for r in original),
            inclusion_statuses=collapse(r['inclusion_status'] for r in original),
            A_dog_count=len({r['animal_id'] for r in original}),
            A_sample_count=len({(r['animal_id'], r['tissue']) for r in original}),
            A_dogs=collapse(r['animal_id'] for r in original),
            A_tissues=collapse(r['tissue'] for r in original),
            background_dog_count=len({c['animal_id'] for c in detected}),
            background_sample_count=len(detected),
            background_dogs=collapse(c['animal_id'] for c in detected),
            background_tissues=collapse(c['tissue'] for c in detected),
            evaluable_sample_count=sum(c['pair_eligible'] is True and c['upstream_call_complete'] for c in pcells),
            unknown_sample_count=sum(c['pair_state'] == 'unknown' for c in pcells),
            not_evaluable_sample_count=sum(c['pair_state'] == 'not_evaluable' for c in pcells)))
    partners = defaultdict(list)
    for p in pairs.values():
        partners[p['partner_gene_id']].append(p)
    shared = [dict(partner_gene_id=partner, partner_symbol=collapse(p['partner_symbol'] for p in ps),
        candidate_count=len({p['new_gene_id'] for p in ps}),
        candidate_ids=collapse(p['new_gene_id'] for p in ps), pair_ids=collapse(p['pair_id'] for p in ps))
        for partner, ps in sorted(partners.items())]

    log('[5/5] save data prepare result with input path, no post-run quality auditing')

    outputs = [('Figure5B_selected_candidates.tsv', candidates, ['candidate_id']),
        ('Figure5B_A_selected_link_records.tsv', selected, ['candidate_id', 'pair_id']),
        ('Figure5B_pair_catalogue.tsv', summaries, ['pair_id']),
        ('Figure5B_pair_sample_evidence.tsv', cells, ['pair_id', 'animal_id', 'tissue', 'pair_state']),
        ('Figure5B_observed_background_links.tsv', raw_rows, ['pair_id']),
        ('Figure5B_shared_partners.tsv', shared, ['partner_gene_id']),
        ('Figure5B_used_endpoint_tasks.tsv', used_tasks, ['task_id'])]
    for name, rows, fallback in outputs:
        fields = list(rows[0]) if rows else fallback
        write_tsv(out / name, rows, fields)
        log(f'{name}: {len(rows)} rows')
    write_tsv(out / 'Figure5B_input_paths.tsv',
        [dict(input_role=k, absolute_path=str(v)) for k, v in INPUTS.items()], ['input_role', 'absolute_path'])
    notes = [
        'Candidates are read directly from module 8 all_candidates and module 5 main_focus; positional-path regional leads come from the main locus manifest, without relying on historical panel A source tables.',
        'Figure5B_selected_candidates retains only identity, scope, and source fields; it does not copy panel A PPA/PIP or derived display columns. A_* fields denote the same integration-rule scope as A, not input read from panel A outputs.',
        'Module 7 only supplements panel A provenance annotations; it does not participate in candidate selection here and is not required as input.',
        'This step prepares data only; it does not generate plots, select representative pairs, or modify PPA/PIP or panel A inclusion criteria.',
        'A_selected_link_records retains expanded integration-source records; row counts must not be interpreted directly as independent support counts.',
        'pair_catalogue defines pairs by candidate and partner gene IDs; A_* counts are kept separate from cross-tissue background_* counts.',
        'Each pair_sample_evidence row represents a pair×dog×tissue; LINK detection in another tissue does not mean that tissue passes molecular–GWAS integration.',
        'detected=called in the original LINK results; not_detected=task completed and evaluable but not called; not_evaluable=task completed but ineligible; unknown=task/sample missing or identity unresolved.',
        'pair_eligible requires each specific endpoint to have >=10 site-level allelic support and to meet the same-autosome ±500 kb window condition; allelic bias is not pre-screened.',
        'Endpoint coordinates retain upstream start/end values (0-based intervals); the window follows LINK\'s inclusive-boundary comparisons and is not converted to variant coordinates.',
        'Neither unknown nor NA is treated as 0. Crak pituitary availability is determined from the manifest; negatives are not imputed to reach an expected 20 samples.',
        'detected prioritizes original calls; observed_endpoint_match and detected_and_eligible are accompanying annotations, not additional exclusion criteria.',
        'within_A_tissue_scope indicates tissue scope; A_selected_link indicates that the pair–sample record actually occurs in A\'s integration sources. These have different meanings.',
        'Partner sharing is assessed by gene ID; it does not imply shared dog–tissue evidence or establish specificity or causal regulation.',
        'enhancing/repressing are only original tool patterns; multi-gene assignment and overlap flags are retained, and phase comparability is not assessed in this step.',
        'Missing molecular information in the Positional path does not indicate a negative molecular result; issues in the original PLASMA implementation are not corrected in this step.',
        'Python 3.10 standard library: csv (TSV and quoted-field parsing), gzip (text decompression), pathlib, datetime, contextlib; one task and one thread.',
        'After execution, separately review output completeness, fields, counts, endpoint assignments, and biological interpretation.',
    ]
    with open(out / 'Figure5B_data_notes.txt', 'x', encoding='utf-8') as handle:
        handle.write('\n'.join(notes) + '\n')


def main():
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ[name] = '1'
    out = OUT_ROOT / datetime.now().strftime('%Y%m%d_%H%M%S')
    out.mkdir(exist_ok=False)
    log_path = out / 'prepare_figure5B.log'
    with open(log_path, 'x', encoding='utf-8', buffering=1) as handle:
        with redirect_stdout(handle), redirect_stderr(handle):
            log(f'output directory={out}; Python={sys.version}; executable={sys.executable}')
            try:
                prepare(out)
                with open(out / 'execution_completion.marker', 'x', encoding='utf-8') as marker:
                    marker.write(f'execution_status=success\noutput_dir={out}\nmaster_log={log_path}\n')
                log('task completed / done; execution_status=success')
                log('execution completed does not equal quality verify passed, please check separately output completeness, format, counts and biological quality.')
            except BaseException:
                log('execution_status=failed')
                traceback.print_exc()
                raise


if __name__ == '__main__':
    main()
