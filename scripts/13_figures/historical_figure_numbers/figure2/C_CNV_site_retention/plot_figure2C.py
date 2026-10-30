#!/path/to/software/miniconda3/envs/dogs/bin/python

import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[name] = '1'
import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from importlib.metadata import version
import logging
import math
import sys
from pathlib import Path
from datetime import datetime
from contextlib import redirect_stdout, redirect_stderr
ROOT = Path('/path/to/repository/scripts/13_figures/historical_figure_numbers/figure2/C_CNV_site_retention')
DOGS = ['Ace_396226', 'CodRou_346627', 'Crak_358336', 'Dozer_386519', 'SanMcg_349891']
TISSUES = ['cerebellum', 'cortex_frontal', 'cortex_occipital', 'pituitary_gland']
SHORT = dict(zip(TISSUES, ['CB', 'FC', 'OC', 'Pit']))


MANIFEST = Path('/path/to/repository/scripts/08_expression_ase_integration/preparation/1.4ase-cnv-qc-filter/log/1.4.ase_cnv_qc_filter.20260821_154914.success.tsv')
FAI = Path('/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna.fai')
PHASE = Path('/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.10_reconcile_dog11_whatsHap_shapeit5/reconciled_bcf')
STATES = {'amp', 'del', 'cnloh'}
FLAGS = {'keep', 'flag_cnv_conflict', 'flag_strong_cnv_overlap', 'flag_soft_cnv_overlap', 'flag_no_cnv_consensus'}
FIELDS = ['HAP1_COUNT', 'HAP2_COUNT', 'REF_COUNT', 'ALT_COUNT', 'TOTAL_COUNT', 'PHASE_STATUS', 'VCF_GT', 'cnv_result_status', 'cnv_evaluable', 'cnv_clean_pass', 'cnv_filter_action', 'cnv_step7_consensus_state', 'cnv_step7_consensus_support', 'cnv_step7_consensus_required_support', 'cnv_filter_reason']


WIDTH_MM = 87.0
HEIGHT_MM = 75.0
DPI = 600
FONT = 'Arial'
FONT_FILES = [Path('/path/to/fonts/Arial/Arial.TTF'),
              Path('/path/to/fonts/Arial/Arialbd.TTF')]


def configure_font():

    from matplotlib import font_manager
    for path in FONT_FILES:
        if not path.is_file():
            raise FileNotFoundError(f'missing Arial font file: {path}')
        font_manager.fontManager.addfont(str(path))
    resolved = {}
    for weight in ('normal', 'bold'):
        path = font_manager.findfont(
            font_manager.FontProperties(family=FONT, weight=weight), fallback_to_default=False)
        if font_manager.FontProperties(fname=path).get_name() != FONT:
            raise ValueError('font is not Arial, do not silently substitute')
        resolved[weight] = path
    if resolved['normal'] == resolved['bold']:
        raise ValueError('Arial regular and bold fonts must not resolve to the same file')
    return resolved


def records(path):
    with Path(path).open(encoding='utf-8', newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')

def key(row):
    chrom = row['CHROM']
    if not chrom.startswith('chr'):
        chrom = 'chr' + chrom
    return chrom, int(row['POS']), row['REF'], row['ALT']

def base(row, lengths):
    k = key(row)
    return (k[0] in lengths and row['PHASE_STATUS'] == 'PHASED_HET_EXACT_MATCH'
            and row['VCF_GT'] in {'0|1', '1|0'} and int(row['TOTAL_COUNT']) > 0)

def boolean(value):
    return {'TRUE': True, 'FALSE': False, '1': True, '0': False}.get(value)

def write_table(path, rows):
    with path.open('x', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def snp_keys(entries, lengths):

    import pysam
    targets = defaultdict(set)
    for e in entries:
        for row in records(e['full_tsv']):
            if base(row, lengths):
                k = key(row)
                targets[k[0]].add(k)
    valid, found = set(), set()
    for chrom in lengths:
        path = PHASE / f'{chrom}.dog11.ref1519.whatsHap_reconciled.phased.bcf'
        with pysam.VariantFile(str(path), drop_samples=True, threads=1) as handle:
            for rec in handle:
                name = rec.contig if rec.contig.startswith('chr') else 'chr' + rec.contig
                for alt in rec.alts or ():
                    candidate = (name, rec.pos, rec.ref, alt)
                    if candidate in targets[chrom]:
                        found.add(candidate)
                if rec.alts is None or len(rec.alts) != 1:
                    continue
                if rec.ref not in {'A', 'C', 'G', 'T'} or rec.alts[0] not in {'A', 'C', 'G', 'T'} or rec.ref == rec.alts[0]:
                    continue
                name = rec.contig if rec.contig.startswith('chr') else 'chr' + rec.contig
                k = (name, rec.pos, rec.ref, rec.alts[0])
                if k in targets[chrom]:
                    valid.add(k)
        logging.info('checking biallelic SNP: %s; target key=%s', chrom, len(targets[chrom]))
    missing = set().union(*targets.values()) - found
    if missing:
        raise ValueError(f'Strict phase keys are missing from the current BCF and cannot be classified as non-SNPs: {len(missing)}; example: {next(iter(missing))}')
    return valid

def summarize(e, lengths, valid, audit, conflicts):
    config = {r['metric']: r['value'] for r in records(e['summary_tsv'])}
    if config.get('clean_mode') != 'known_cnv_removed' or float(config.get('consensus_fraction', 'nan')) != 0.75:
        raise ValueError(f'filter mode does not satisfy current scope: {e["summary_tsv"]}')
    seen, rejected, retained, before = {}, set(), {}, set()
    counts = Counter()
    for row in records(e['full_tsv']):
        if row['SAMPLE'] != e['animal_id'] or row['TISSUE'] != e['term']:
            raise ValueError('full table sample identity mismatch')
        k = key(row)
        signature = tuple(row[f] for f in FIELDS)
        if k in seen:
            if seen[k] != signature:
                conflicts.writerow([e['animal_id'], e['term'], *k, repr(seen[k]), repr(signature)])
                raise ValueError(f'duplicate variant record conflict: {k}; see conflicts.tsv')
            continue
        seen[k] = signature
        if not base(row, lengths):
            continue
        if not 1 <= k[1] <= lengths[k[0]]:
            raise ValueError(f'coordinates exceeds FAI: {k}')
        if k not in valid:
            rejected.add(k)
            continue
        h1, h2, ref, alt, total = (int(row[f]) for f in FIELDS[:5])
        if min(h1, h2, ref, alt) < 0 or h1 + h2 != total or ref + alt != total:
            raise ValueError(f'allelic count mismatch: {k}')
        counts['N_all_informative'] += 1
        evaluable, keep = boolean(row['cnv_evaluable']), boolean(row['cnv_clean_pass'])
        action = row['cnv_filter_action']
        if row['cnv_result_status'] == 'OK' and evaluable is True:
            before.add(k)
            if keep is None or (keep and action not in FLAGS) or (not keep and action != 'exclude_cnv_overlap'):
                raise ValueError(f'evaluable set contains unknown/mismatch determine: {k}')
            if keep:
                retained[k] = signature
                decision = 'retained'
            else:
                if row['cnv_step7_consensus_state'] not in STATES or int(row['cnv_step7_consensus_required_support']) != 7 or not 7 <= int(row['cnv_step7_consensus_support']) <= 9:
                    raise ValueError(f'exclude reason does not satisfy 7/9 rule: {k}')
                counts['N_removed_cnv'] += 1
                decision = 'removed_cnv'
        elif evaluable is False and row['cnv_result_status'] in {'NO_STEP7_SEGMENT_RESULT', 'NO_STEP7_SAMPLE_RESULT', 'INCOMPLETE_STEP7_SCALE_SET'}:
            counts['N_cnv_not_evaluable'] += 1
            decision = 'not_evaluable'
        else:
            counts['N_cnv_error_or_unknown'] += 1
            decision = 'error_or_unknown'
        audit.writerow([e['animal_id'], e['term'], *k, h1, h2, total,
                        row['cnv_result_status'], row['cnv_evaluable'], row['cnv_clean_pass'], action,
                        row['cnv_filter_reason'], decision, e['full_tsv']])

    clean = set()
    for row in records(e['clean_tsv']):
        if row['SAMPLE'] != e['animal_id'] or row['TISSUE'] != e['term']:
            raise ValueError('clean table sample identity mismatch')
        if not base(row, lengths) or key(row) not in valid:
            continue
        k = key(row)
        if row['cnv_result_status'] != 'OK' or boolean(row['cnv_evaluable']) is not True:
            continue
        if k not in retained or tuple(row[f] for f in FIELDS) != retained[k]:
            conflicts.writerow([e['animal_id'], e['term'], *k, repr(retained.get(k)), repr(tuple(row[f] for f in FIELDS))])
            raise ValueError(f'clean exists additional or is modify of variant: {k}')
        clean.add(k)
    if clean != set(retained) or not clean <= before or len(before) != len(clean) + counts['N_removed_cnv']:
        missing = set(retained) - clean
        for k in sorted(missing):
            conflicts.writerow([e['animal_id'], e['term'], *k, 'expected retained', 'absent from clean'])
        raise ValueError('Before/After set comparison failed; see conflicts.tsv')
    n = len(before)
    return dict(Dog=e['animal_id'], Tissue=e['term'], scope='all_informative_biallelic_snps',
                N_phase_coverage_nonbiallelic_or_nonSNP=len(rejected),
                N_all_informative=counts['N_all_informative'], N_cnv_not_evaluable=counts['N_cnv_not_evaluable'],
                N_cnv_error_or_unknown=counts['N_cnv_error_or_unknown'], N_before=n,
                N_removed_cnv=counts['N_removed_cnv'], N_after=len(clean),
                retention_pct=100*len(clean)/n if n else 'NA',
                removed_pct=100*counts['N_removed_cnv']/n if n else 'NA', count_qc_status='pass')

def summarize_from_raw(run):
    entries = list(records(MANIFEST))

    expected = {(d, t) for d in DOGS for t in TISSUES} - {('Crak_358336', 'pituitary_gland')}
    if len(entries) != 19 or {(e['animal_id'], e['term']) for e in entries} != expected:
        raise ValueError('sample list does not satisfy 19 items actual combination')
    entries.sort(key=lambda e: (DOGS.index(e['animal_id']), TISSUES.index(e['term'])))
    with FAI.open() as handle:
        all_lengths = {v[0]: int(v[1]) for line in handle if (v := line.split())}
    lengths = {f'chr{i}': all_lengths[f'chr{i}'] for i in range(1, 39)}
    inputs = [MANIFEST, FAI, Path(__file__).resolve()]
    inputs += [Path(e[f]) for e in entries for f in ('full_tsv', 'clean_tsv', 'summary_tsv')]
    inputs += [PHASE / f'{c}.dog11.ref1519.whatsHap_reconciled.phased.bcf' for c in lengths]
    provenance = [dict(path=str(p), size=p.stat().st_size, mtime_ns=p.stat().st_mtime_ns) for p in inputs]
    (run / 'input_parameters.json').write_text(json.dumps(dict(inputs=provenance, python=sys.version,
        packages={x: version(x) for x in ('pysam', 'matplotlib')}, scope='all_informative_biallelic_snps',
        total_count_min_exclusive=0, consensus_required=7, consensus_scales=9,
        threads=1, audit_requested_by='Panel C user brief'), indent=2), encoding='utf-8')
    logging.info('start; concurrent_tasks=1; threads_per_task=1; check phased BCF SNP type')
    valid = snp_keys(entries, lengths)
    with gzip.open(run / 'panel_C_cnv_retention.audit.tsv.gz', 'wt', encoding='utf-8', newline='') as ah, (run / 'conflicts.tsv').open('x', encoding='utf-8', newline='') as ch:
        audit, conflicts = csv.writer(ah, delimiter='\t'), csv.writer(ch, delimiter='\t')
        audit.writerow(['Dog','Tissue','CHROM','POS','REF','ALT','HAP1_COUNT','HAP2_COUNT','TOTAL_COUNT','cnv_result_status','cnv_evaluable','cnv_clean_pass','cnv_filter_action','reason','classification','source'])
        conflicts.writerow(['Dog','Tissue','CHROM','POS','REF','ALT','expected','observed'])
        summaries = []
        for e in entries:
            result = summarize(e, lengths, valid, audit, conflicts)
            summaries.append(result)
            logging.info('%s/%s Before=%s After=%s Removed=%s', result['Dog'], result['Tissue'], result['N_before'], result['N_after'], result['N_removed_cnv'])
    summary = run / 'panel_C_cnv_retention.summary.tsv'
    write_table(summary, summaries)
    logging.info('statistics with set audit completed; plot only consume summary.tsv')
    draw(summary, run)


def draw(summary, run):
    logging.info('read C figure summary: %s; single task, numerical-library thread counts remain 1', summary)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    with summary.open(encoding='utf-8', newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    expected = {(d, t) for d in DOGS for t in TISSUES} - {('Crak_358336', 'pituitary_gland')}
    if len(rows) != 19 or {(r['Dog'], r['Tissue']) for r in rows} != expected:
        raise ValueError('summary table is not current 19 row design')
    rows.sort(key=lambda r: (DOGS.index(r['Dog']), TISSUES.index(r['Tissue'])))
    for r in rows:
        if r['count_qc_status'] != 'pass':
            raise ValueError('statistics audit failed, do not generate the final figure')
        r['before'], r['after'], r['removed'] = (int(r[k]) for k in ('N_before','N_after','N_removed_cnv'))
        if not 0 <= r['after'] <= r['before'] or r['before'] != r['after'] + r['removed']:
            raise ValueError('summary count does not satisfy subset relationship')
    import numpy as np
    import hashlib
    import json
    from matplotlib.colors import LinearSegmentedColormap


    matrix = np.full((len(DOGS), len(TISSUES)), np.nan)
    for r in rows:
        if r['before'] > 0:
            matrix[DOGS.index(r['Dog']), TISSUES.index(r['Tissue'])] = 100*r['removed']/r['before']
    finite = matrix[np.isfinite(matrix)]
    if not finite.size:
        raise ValueError('No samples with a positive denominator')
    vmax = max(12, 2*math.ceil(float(finite.max())/2))
    cmap = LinearSegmentedColormap.from_list('removal_blue', ['#F2F6FA', '#719FC1', '#174A73'])
    cmap.set_bad('#DADDE1')
    font_paths = configure_font()
    plt.rcParams.update({'font.family':'sans-serif',
        'font.sans-serif':[FONT],
        'font.size':7, 'pdf.fonttype':42, 'svg.fonttype':'none',
        'figure.facecolor':'white', 'savefig.facecolor':'white',
        'savefig.bbox':None, 'savefig.transparent':False})

    fig = plt.figure(figsize=(WIDTH_MM/25.4, HEIGHT_MM/25.4))

    grid = fig.add_gridspec(1, 3, left=0.27, right=0.93,
                           bottom=8/HEIGHT_MM, top=62/HEIGHT_MM,
                           width_ratios=[4.2, 1, 0.12], wspace=0.16)
    count_ax = fig.add_subplot(grid[0,0])

    legend_ax = fig.add_axes([0.19, 0.5/HEIGHT_MM, 0.79, 7/HEIGHT_MM])
    legend_ax.axis('off')
    from matplotlib.lines import Line2D

    positions, group_centres, separators = [], [], []
    cursor = 0
    group_gap = 0.5
    for dog in DOGS:
        group_rows = [r for r in rows if r['Dog'] == dog]
        ys = [cursor + j for j in range(len(group_rows))]
        positions.extend(ys)
        group_centres.append((dog.split('_')[0], (ys[0] + ys[-1])/2))
        cursor += len(group_rows) + group_gap

        separators.append(cursor - (1 + group_gap)/2)
    xmax = max(50000, 50000 * math.ceil(max(r['before'] for r in rows)/50000))
    for r, y in zip(rows, positions):
        count_ax.plot([r['after'], r['before']], [y,y], color='#B8BEC3', lw=0.7, zorder=1)
        count_ax.scatter(r['before'], y, s=25, facecolors='white', edgecolors='#707070', linewidths=1, zorder=2)
        count_ax.scatter(r['after'], y, s=11, color='#0072B2', edgecolors='none', zorder=3)
        count_ax.text(-0.035,y,SHORT[r['Tissue']],transform=count_ax.get_yaxis_transform(),
                      ha='right',va='center',fontsize=7)
    for name,y in group_centres:
        count_ax.text((2/WIDTH_MM-count_ax.get_position().x0)/count_ax.get_position().width,y,name,transform=count_ax.get_yaxis_transform(),
                      ha='left',va='center',fontsize=7,fontweight='bold')
    for y in separators[:-1]:
        count_ax.axhline(y,color='#E5E7E9',lw=0.5)
    count_ax.set_xlim(0,xmax)
    count_ax.set_ylim(positions[-1]+0.7,-0.9)
    count_ax.set_yticks([])
    ticks = list(range(0,xmax+1,50000))
    count_ax.set_xticks(ticks,['0' if v==0 else f'{v//1000}k' for v in ticks])
    count_ax.xaxis.tick_top()
    count_ax.xaxis.set_label_position('top')
    count_ax.set_xlabel('Number of informative SNPs',fontsize=7.5)

    count_ax.xaxis.set_label_coords(0.5,1.105)
    count_ax.tick_params(axis='x',labelsize=7,length=3,width=0.6)
    for side in ('left','right','bottom'):
        count_ax.spines[side].set_visible(False)
    count_ax.spines['top'].set_linewidth(0.6)
    handles = [Line2D([],[],marker='o',linestyle='none',markerfacecolor='white',
                       markeredgecolor='#707070',markersize=4,label='Before CNV filtering'),
               Line2D([],[],marker='o',linestyle='none',color='#0072B2',markersize=3,
                       label='After CNV filtering')]
    legend_ax.legend(handles=handles,loc='center',ncol=2,frameon=False,fontsize=7,
                     handlelength=1.0,handletextpad=0.35,columnspacing=1.1)
    from matplotlib.patches import Rectangle
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    ax = fig.add_subplot(grid[0,1], sharey=count_ax)
    norm = Normalize(vmin=0, vmax=vmax)
    for r,y in zip(rows,positions):
        value = matrix[DOGS.index(r['Dog']), TISSUES.index(r['Tissue'])]
        face = '#DADDE1' if not np.isfinite(value) else cmap(norm(value))
        ax.add_patch(Rectangle((0,y-0.48),1,0.96,facecolor=face,edgecolor='white',linewidth=0.5))
        label = 'NA' if not np.isfinite(value) else ('<0.01' if 0 < value < 0.005 else f'{value:.2f}')
        color = 'white' if np.isfinite(value) and value/vmax > 0.65 else '#172B3A'

        ax.text(0.5,y,label,ha='center',va='center',fontsize=6.5,color=color)
    ax.set_xlim(0,1)
    ax.set_xticks([])
    ax.tick_params(axis='y',left=False,labelleft=False)
    ax.set_title('Removed (%)',fontsize=7.5,y=1.105,pad=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    scale_grid = grid[0,2].subgridspec(3, 1, height_ratios=[1,1,1], hspace=0)
    cax = fig.add_subplot(scale_grid[1,0])

    cb = fig.colorbar(ScalarMappable(norm=norm,cmap=cmap), cax=cax,
                      orientation='vertical',ticks=[0, vmax/2, vmax])
    cb.ax.tick_params(labelsize=6.5,width=0.5,length=2)
    cb.outline.set_visible(False)
    fig.text(2/WIDTH_MM, 1-2/HEIGHT_MM, 'C', fontsize=10, fontweight='bold', va='top')

    sys.path.insert(0, '/path/to/repository/support/figure_alignment')
    from audit_panel_alignment import require_matplotlib_panel_alignment
    require_matplotlib_panel_alignment(
        fig, json_out=run/'Figure2C.alignment.json',
        exclude_axes=[legend_ax, cax], tolerance_pt=1.5, strict=True)
    fig.savefig(run/'Figure2C_CNV_site_retention.pdf')
    fig.savefig(run/'Figure2C_CNV_site_retention.svg')
    fig.savefig(run/'Figure2C_CNV_site_retention.png', dpi=DPI)

    from PIL import Image
    with Image.open(run/'Figure2C_CNV_site_retention.png') as raster:
        with raster.convert('RGB') as rgb:
            rgb.save(run/'Figure2C_CNV_site_retention.tiff', compression='tiff_lzw', dpi=(DPI,DPI))
    plt.close(fig)


    with (run/'Figure2C.removed_percentage_matrix.tsv').open('w') as h:
        writer = csv.writer(h, delimiter='\t')
        writer.writerow(['Dog']+[SHORT[t] for t in TISSUES])
        for i, dog in enumerate(DOGS):
            writer.writerow([dog]+['NA' if not np.isfinite(v) else str(v) for v in matrix[i]])
    (run/'Figure2C.render_parameters.json').write_text(json.dumps({
        'input_summary':str(summary.resolve()), 'input_sha256':hashlib.sha256(summary.read_bytes()).hexdigest(),
        'metric':'100 * N_removed_cnv / N_before', 'color_scale':'linear',
        'source_mode':'raw_1_4_tables_and_phased_bcf', 'layout':'counts_left_shared_rows_removal_column', 'count_metrics':['N_before','N_after'],
        'count_axis_scale':'linear', 'count_axis_limits':[0,xmax], 'retention_column':False,
        'cell_height_row_units':0.96, 'colorbar_height_fraction':1/3, 'colorbar_ticks':[0,vmax/2,vmax],
        'data_area_bottom_mm':8, 'data_area_top_mm':62,
        'data_area_height_mm':54, 'group_gap_row_units':group_gap,
        'layout_role':'second_row_75mm_trial',
        'panel_wspace':0.16, 'column_title_y':1.105, 'dog_label_left_mm':2,
        'vmin':0, 'vmax':vmax, 'width_mm':WIDTH_MM, 'height_mm':HEIGHT_MM,
        'font_family':FONT, 'resolved_font_files':font_paths, 'raster_dpi':DPI,
        'font_sizes_pt':{'panel':10,'titles':7.5,'body':7,'percentage':6.5,'colorbar':6.5},
        'quality_review':{'font_embedding':'NOT_CHECKED','collisions':'NOT_CHECKED',
                          'figure2_assembly':'NOT_CHECKED','figure_and_caption':'NOT_CHECKED'},
        'observed_samples':len(rows), 'unavailable_sample':'Crak pituitary',
        'row_normalization':False, 'upstream_analysis_rerun':False}, indent=2)+'\n')
    caption = (
        'C. Informative SNP counts and removal percentages following CNV filtering. '
        'Left: each row is one of 19 dog–tissue combinations. Grey open circles show '
        'N_before and smaller blue circles show N_after on a zero-based linear count axis; '
        'connecting lines indicate paired counts without positional offsets. '
        'Right: a single removal-percentage column shares the same 19 sample rows. '
        'Tissues are cerebellum (CB), frontal cortex (FC), occipital cortex (OC), '
        'and pituitary (Pit). Cell values and the common linear color scale show '
        '100 x N_removed_cnv / N_before, without row normalization. '
        'The denominator comprises unique autosomal biallelic SNPs with exact phased heterozygous '
        'matching, positive RNA allele counts and evaluable CNV status. Filtering removes sites with '
        'a single amplification, deletion or copy-neutral LOH consensus supported by at least seven '
        'of nine analysis scales. N_before = N_after + N_removed_cnv. '
        'Crak pituitary is absent and has no displayed row. A zero denominator would be shown as grey NA. '
        '0.00 denotes no removal; positive values that would round to zero are labelled <0.01. '
        'The 19 observed samples are shown separately, with no significance testing. '
        'Removal percentages do not estimate genome-wide CNV burden, and retention does not establish '
        'diploid copy number or statistically significant ASE.\n'
    )
    notes = (
        'Summary recomputed from 1.4 full/clean tables (run 20260821_154914) and phased BCF '
         'before rendering. This is the same CNV-filtering run used by panel B and upstream '
         'of the 1.5 gene-level tables used by panel A. '
        'Upstream CNV calls and filtering policy are unchanged. '
        'Each sample/CHROM/POS/REF/ALT is counted once, including intergenic sites. '
        'The known_cnv_removed policy uses cnv_clean_pass; soft/strong/conflict flags are not '
        'automatically excluded. Percentages are recomputed from integer counts. '
        'The left paired-count plot and right removal column use the same summary; full counts remain in the summary generated in this run. '
        'Pattern reference: Matplotlib official annotated heatmap example; independent implementation.\n'
        'https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html\n'
        'https://github.com/matplotlib/matplotlib/blob/main/galleries/examples/images_contours_and_fields/image_annotated_heatmap.py\n'
    )
    (run/'Figure2C_CNV_site_retention.caption.md').write_text(caption, encoding='utf-8')
    (run/'Figure2C_CNV_site_retention.methods_notes.md').write_text(notes, encoding='utf-8')
    logging.info('left and right shared sample row of combination figure, shared source table, final figure legends and technical notes have been written')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarize the original step 1.4 table and phased BCF to plot panel C.')
    parser.parse_args()
    run = ROOT/datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    run.mkdir(parents=True, exist_ok=False)
    os.environ['MPLCONFIGDIR'] = str(run/'.matplotlib')
    with (run/'figure2C.log').open('x', encoding='utf-8') as log:
        with redirect_stdout(log), redirect_stderr(log):
            logging.basicConfig(stream=log, level=logging.INFO, format='[%(asctime)s] %(message)s')
            try:
                summarize_from_raw(run)
                (run/'execution_completion.marker').write_text('execution_status=success\n')
                logging.info('task completed / done; execution_status=success; please check separately layout and biological interpretation')
            except Exception:
                logging.exception('execution_status=failed')
                sys.exit(1)
