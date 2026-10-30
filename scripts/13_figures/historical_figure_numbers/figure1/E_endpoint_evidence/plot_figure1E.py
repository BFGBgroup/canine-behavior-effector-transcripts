#!/path/to/software/miniconda3/envs/dogs/bin/python

import os
for thread_key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                   'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[thread_key] = '1'
import logging
import contextlib
import sys
import argparse
import csv
import json
from datetime import datetime
from pathlib import Path

ROOT = Path('/path/to/repository/scripts/13_figures/historical_figure_numbers/figure1/E_endpoint_evidence')
BASE = Path('/path/to/project/analysis_20260709/5-dogs/9.new-gene')
TSS = BASE / '2.2_candidate_TSS_promoter_support/20260921_184512/all_candidate_transcript_TSS_promoter_support.tsv'
TES = BASE / '2.3_candidate_TES_polyA_support/20260921_204439/all_candidate_transcript_TES_polyA_support.tsv'
MAP = Path('/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.transcript_map.tsv')
METRICS = (
    ('5prime_cage_proximity', '5′', 'CAGE proximity (≤500 bp)'),
    ('5prime_context', '5′', 'Promoter context'),
    ('3prime_external', '3′', 'Poly(A)-site proximity (≤500 bp)'),
    ('3prime_motif', '3′', 'PAS motif (terminal 50 nt; selected11)'),
)


WIDTH_MM = 87.0
HEIGHT_MM = 60.0
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


def indexed(path, key):
    with path.open(newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    result = {}
    for row in rows:
        if not row[key] or row[key] in result:
            raise ValueError(f'Empty or duplicate {key} in {path}')
        result[row[key]] = row
    if not result:
        raise ValueError(f'Empty input: {path}')
    return result


def status(row, field, available):
    value = row[field]
    if value not in (('yes', 'no') if available else ('NA',)):
        raise ValueError(f'Invalid {field}={value} for available={available}')
    return value


def collect(mapping, tss, tes):
    records, seen = [], set()
    for tid, m in mapping.items():
        old = m['original_transcript_id']
        if old in seen or m['source_dataset'] != 'five-dogs':
            raise ValueError(f'Nonunique source or wrong dataset: {tid}')
        seen.add(old)
        a, b = tss[old], tes[old]
        for source in (a, b):
            for mf, sf in (('chr', 'transcript_chr'), ('strand', 'transcript_strand'),
                           ('transcript_start_0based', 'transcript_span_start0'),
                           ('transcript_end_0based', 'transcript_span_end0')):
                if m[mf] != source[sf]:
                    raise ValueError(f'Model mismatch: {tid} {mf}')
        if a['reference_assessable'] not in ('yes', 'no') or b['reference_assessable'] not in ('yes', 'no') or b['spliced_pas_scan_available'] not in ('0', '1'):
            raise ValueError(f'Invalid assessability: {tid}')
        av, bv, mv = a['reference_assessable'] == 'yes', b['reference_assessable'] == 'yes', b['spliced_pas_scan_available'] == '1'
        row = dict(final_transcript_id=tid, final_gene_id=m['final_gene_id'], original_transcript_id=old)
        for key, source, field, available in (
            ('5prime_cage_proximity', a, 'has_cage_proximity_support_500bp', av),
            ('5prime_context', a, 'has_promoter_context_support', av),
            ('3prime_external', b, 'external_within_500bp', bv),
            ('3prime_motif', b, 'spliced_pas_extended11_found', mv),
        ):
            row[key] = status(source, field, available)
        if av:
            for out, fields in (
                ('5prime_cage_proximity', ('doga_promoter_within_500bp',)),
                ('5prime_context', ('accessible_context_supported', 'promoter_like_context_supported')),
            ):

                components = [status(a, field, True) for field in fields]
                union = any(value == 'yes' for value in components)
                if (row[out] == 'yes') != union:
                    raise ValueError(f'Inconsistent support union: {tid} {out}')

        for field in ('external_within_50bp', 'external_within_100bp',
                      'polyaseq_within_500bp', 'animal_apadb_within_500bp'):
            row[field] = status(b, field, bv)
        row['spliced_pas_selected9_found'] = status(b, 'spliced_pas_selected9_found', mv)
        row['has_tss_position_support'] = status(a, 'has_tss_position_support', av)
        for field in ('source_support_class_500bp', 'spliced_pas_scan_failure_reason',
                      'spliced_pas_hit_list', 'downstream_20nt_A_percent'):
            row[field] = b[field]
        records.append(row)
    metrics = []
    for key, end, label in METRICS:
        n = sum(r[key] != 'NA' for r in records)
        k = sum(r[key] == 'yes' for r in records)
        if n == 0:
            raise ValueError(f'No assessable transcripts: {key}')
        metrics.append(dict(metric=key, endpoint=end, label=label, supported=k,
                            assessable=n, not_assessed=len(records)-n, percent=100*k/n))
    joint = {}
    for prefix, x, y in (('5prime', '5prime_cage_proximity', '5prime_context'), ('3prime', '3prime_external', '3prime_motif')):
        joint[prefix] = dict(supported=sum(r[x] == r[y] == 'yes' for r in records),
                             assessable=sum(r[x] != 'NA' and r[y] != 'NA' for r in records))
    return records, metrics, joint


def write_tsv(path, rows):
    with path.open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def draw(run, metrics):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    font_paths = configure_font()
    matplotlib.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': [FONT],
                                'font.size': 7, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
                                'axes.linewidth': 0.6, 'savefig.bbox': None,
                                'savefig.facecolor': 'white', 'savefig.transparent': False})
    ink, gray, point_color = '#102030', '#687583', '#3B73B9'
    fig, ax = plt.subplots(figsize=(WIDTH_MM/25.4, HEIGHT_MM/25.4))
    fig.subplots_adjust(left=0.47, right=0.97, bottom=0.20, top=0.90)
    fig.text(2/WIDTH_MM, 1-2/HEIGHT_MM, 'E', va='top', weight='bold', fontsize=10, color=ink)
    ys = [4.5, 3.5, 1.5, 0.5]
    names = ['CAGE proximity\n(≤500 bp)', 'Promoter context',
             'Poly(A)-site proximity\n(≤500 bp)', 'PAS motif\n(terminal 50 nt)']
    values = [m['percent'] for m in metrics]

    ax.scatter(values, ys, s=22, color=point_color, linewidths=0,
               zorder=3, clip_on=False)
    ax.set_xlim(0, max(80, max(values) + 20))
    ax.set_ylim(-0.3, 5.8)
    ax.set_xticks([0, 20, 40, 60])
    ax.set_yticks(ys, names)
    ax.tick_params(axis='x', length=2.5, width=0.6, pad=3, labelsize=7, colors=gray)
    ax.tick_params(axis='y', length=0, pad=7, labelsize=7, colors=ink)
    for spine in ('top', 'right', 'left'):
        ax.spines[spine].set_visible(False)
    ax.spines['bottom'].set_color(gray)
    for y, value in zip(ys, values):
        ax.annotate(f'{value:.2f}%', (value, y), xytext=(5, 0),
                    textcoords='offset points', ha='left', va='center', fontsize=7, color=ink)
    for y, label in ((5.3, '5′ evidence'), (2.3, '3′ evidence')):
        ax.text((2/WIDTH_MM - 0.47)/0.50, y, label, transform=ax.get_yaxis_transform(),
                ha='left', va='center', fontsize=7, weight='bold', color=ink)
    ax.set_xlabel('Supported transcripts (%)', fontsize=7.5, color=ink, labelpad=5)
    fig.savefig(run/'Figure1E_endpoint_evidence.pdf')
    fig.savefig(run/'Figure1E_endpoint_evidence.svg')
    fig.savefig(run/'Figure1E_endpoint_evidence.png', dpi=DPI)

    from PIL import Image
    with Image.open(run/'Figure1E_endpoint_evidence.png') as raster:
        rgb = raster.convert('RGB')
        rgb.save(run/'Figure1E_endpoint_evidence.tiff', compression='tiff_lzw',
                 dpi=(DPI, DPI))
        rgb.close()
    plt.close(fig)
    return font_paths


def main(run):
    parser = argparse.ArgumentParser(description='Run plot_figure1E.py; see the workflow documentation.')
    parser.add_argument('--map', type=Path, default=MAP)
    parser.add_argument('--tss', type=Path, default=TSS)
    parser.add_argument('--tes', type=Path, default=TES)
    args = parser.parse_args()
    logging.info('[stage 1/3] read identifier mapping and TSS/TES input; concurrent_tasks=1 threads_per_task=1')
    mapping = indexed(args.map, 'final_transcript_id')
    records, metrics, joint = collect(mapping, indexed(args.tss, 'query_transcript_id'), indexed(args.tes, 'query_transcript_id'))
    logging.info('[stage 2/3] summary%d records transcript of evidence and write source table: %s', len(records), run)
    run.mkdir()
    write_tsv(run/'Figure1E.transcripts.tsv', records)
    write_tsv(run/'Figure1E.support_counts.tsv', metrics)
    provenance = dict(inputs={k: {'path': str(p.resolve()), 'size_bytes': p.stat().st_size} for k, p in vars(args).items()},
                      total_transcripts=len(records), within_endpoint_intersections=joint,
                      unit='transcript', width_mm=WIDTH_MM, height_mm=HEIGHT_MM,
                      font_family=FONT, raster_dpi=DPI, tiff_mode='RGB',
                      tiff_compression='LZW',
                      font_sizes_pt={'panel': 10, 'group': 7, 'labels': 7,
                                     'percent': 7, 'ticks': 7, 'axis_title': 7.5},
                      panel_anchor_mm={'left': 2, 'top': 2},
                      quality_review={'pdf_font_embedding': 'NOT_CHECKED',
                                      'text_collisions': 'NOT_CHECKED',
                                      'assembled_figure1': 'NOT_CHECKED',
                                      'figure_and_caption_fit': 'NOT_CHECKED'},
                      tss_proximity_bp=500, tes_proximity_bp=500,
                      pas_definition='spliced terminal up-to-50 nt; selected11',
                      plot_type='grouped horizontal dot plot',
                      panel_alignment='Single quantitative panel; visual review required after plotting.')
    (run/'inputs_and_parameters.json').write_text(json.dumps(provenance, indent=2)+'\n')

    comparisons = []
    for key in ('has_tss_position_support', 'external_within_50bp', 'external_within_100bp',
                'polyaseq_within_500bp', 'animal_apadb_within_500bp', 'spliced_pas_selected9_found'):
        n = sum(r[key] != 'NA' for r in records)
        k = sum(r[key] == 'yes' for r in records)
        comparisons.append(dict(metric=key, supported=k, assessable=n,
                                not_assessed=len(records)-n, percent=100*k/n if n else None))
    write_tsv(run/'Figure1E.comparison_counts.tsv', comparisons)
    caption = (
        f'E. Endpoint evidence for {len(records):,} retained candidate transcript models. '
        'Points show the percentage supported among transcripts assessable for each metric; isoforms are counted separately. '
        'CAGE proximity: inferred TSS within 500 bp of a covered base of a same-strand DoGA robust promoter interval, using the all-promoter reference converted to CanFam3.1. '
        'Promoter context: ATAC or EpiC-Dog promoter-like states 1–4 overlapping the strand-oriented −1,000 to +100 bp window; interval strand is not required. '
        'Poly(A)-site proximity: a same-strand GSE30198 PolyA-seq or Animal-APAdb reference site within 500 bp of the inferred TES. Experimental and inferred reference sources are combined, not independent validations. '
        'PAS motif: a complete selected hexamer in the final up-to-50 nt of the spliced transcript, including its terminal base. '
        'The selected11 panel is AATAAA, ATTAAA, TATAAA, AGTAAA, AAGAAA, AATATA, AATACA, CATAAA, GATAAA, AATGAA and AATAGA; it is not an exhaustive or dog-validated panel. '
        'Non-primary contigs and unavailable sequence/model assessments remain unassessed, not negative. '
        + ' '.join(f"{m['label']}: {m['supported']:,}/{m['assessable']:,} ({m['percent']:.2f}%); {m['not_assessed']:,} not assessed." for m in metrics)
        + f" Within-5′ intersection: {joint['5prime']['supported']:,}/{joint['5prime']['assessable']:,}; within-3′ intersection: {joint['3prime']['supported']:,}/{joint['3prime']['assessable']:,}. These are not cross-endpoint intersections. "
        + 'Distance/source comparisons: '
        + '; '.join(f"{m['metric']}={m['supported']:,}/{m['assessable']:,}" for m in comparisons)
        + '. Categories overlap and must not be summed. Proximity and motif evidence do not validate precise cleavage sites or full-length structure; no internal-priming filter was applied. '
        'Percentages do not compare endpoint accuracy. No significance tests or sampling uncertainty estimates are shown.\n'
    )
    (run/'Figure1E.caption.txt').write_text(caption)
    logging.info('[stage 3/3] export 87 × 60 mm PDF/SVG, 600 dpi PNG and RGB LZW TIFF')
    provenance['resolved_font_files'] = draw(run, metrics)
    (run/'inputs_and_parameters.json').write_text(json.dumps(provenance, indent=2)+'\n')
    (run/'execution_completion.marker').write_text('execution_status=success\n')
    logging.info('task completed / done execution_status=success output=%s', run)
    logging.info('output completeness, statistical denominator and plot layout still requires separately check.')


if __name__ == '__main__':

    run_id = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    log_dir = ROOT/'log'
    log_dir.mkdir(exist_ok=True)
    with (log_dir/f'{run_id}.master.log').open('x') as log_handle:
        with contextlib.redirect_stdout(log_handle), contextlib.redirect_stderr(log_handle):
            logging.basicConfig(stream=log_handle, level=logging.INFO,
                                format='[%(asctime)s ] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
            try:
                main(ROOT/run_id)
            except BaseException:
                logging.exception('execution_status=failed; this run not completed')
                raise
