#!/path/to/software/miniconda3/envs/dogs/bin/python

import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'

import csv
import logging
import sys
from collections import Counter
from contextlib import redirect_stdout, redirect_stderr
from datetime import datetime
from pathlib import Path

ROOT = Path('/path/to/repository/scripts/13_figures/historical_figure_numbers/figure1/D_predicted_biotype_composition')
MAP = Path('/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.transcript_map.tsv')
CLASSES = ('lncRNA-only', 'mRNA-only', 'Mixed')
COLORS = ('#719FC1', '#DFA08F', '#BDBDBD')
INK, GRAY = '#102030', '#687583'
WIDTH_MM, HEIGHT_MM, DPI = 87.0, 60.0, 600
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


def read_records():
    records, groups, seen = [], {}, set()
    with MAP.open(newline='', encoding='utf-8') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            tid, gid, biotype = row['final_transcript_id'], row['final_gene_id'], row['FEELnc_biotype']
            if row['source_dataset'] != 'five-dogs':
                raise ValueError(f'non-five dogs input record: {tid}')
            if not tid or not gid or tid in seen:
                raise ValueError(f'blank identifier or duplicate transcript: {tid}')
            if biotype not in ('lncRNA', 'mRNA'):
                raise ValueError(f'unknown or missing FEELnc type, cannot automatically classify as Mixed: {tid} {biotype}')
            seen.add(tid)
            records.append(row)
            groups.setdefault(gid, []).append(row)
    loci = []
    for gid, members in groups.items():
        types = {r['FEELnc_biotype'] for r in members}
        category = 'lncRNA-only' if types == {'lncRNA'} else (
            'mRNA-only' if types == {'mRNA'} else 'Mixed')
        loci.append(dict(
            final_gene_id=gid, predicted_biotype_class=category,
            retained_transcript_count=len(members),
            lncRNA_transcript_count=sum(r['FEELnc_biotype'] == 'lncRNA' for r in members),
            mRNA_transcript_count=sum(r['FEELnc_biotype'] == 'mRNA' for r in members),
            noORF_transcript_count=sum(r['FEELnc_codpot_category'] == 'noORF' for r in members),
            final_transcript_id_list='|'.join(r['final_transcript_id'] for r in members)))
    return records, loci


def make_outputs(run, records, loci):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    if not loci:
        raise ValueError('candidate site set is empty, cannot compute proportion')
    totals = Counter(r['predicted_biotype_class'] for r in loci)
    n = len(loci)

    noorf_transcripts = sum(r['noORF_transcript_count'] for r in loci)
    noorf_loci = sum(r['noORF_transcript_count'] > 0 for r in loci)
    logging.info('transcripts=%d loci=%d classes=%s', len(records), n, dict(totals))
    logging.info('upstream_noORF_transcripts=%d affected_loci=%d; retained_under_lncRNA',
                 noorf_transcripts, noorf_loci)
    for name, rows in (('transcripts', records), ('loci', loci)):
        with (run/f'Figure1D.biotype_{name}.tsv').open('x', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
    with (run/'Figure1D.biotype_counts.tsv').open('x', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['predicted_biotype_class', 'candidate_loci', 'total_candidate_loci', 'percent'])
        for category in CLASSES:
            writer.writerow([category, totals[category], n, f'{100*totals[category]/n:.6f}'])

    fonts = configure_font()
    logging.info('font_family=%s normal=%s bold=%s', FONT, fonts['normal'], fonts['bold'])
    matplotlib.rcParams.update({'font.family': 'sans-serif',
                               'font.sans-serif': [FONT],
                               'font.size': 7, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                               'svg.fonttype': 'none', 'savefig.bbox': None})
    from math import cos, sin, radians
    fig = plt.figure(figsize=(WIDTH_MM/25.4, HEIGHT_MM/25.4), facecolor='white')


    ax = fig.add_axes([0.27, 0.10, 0.70, 0.78])
    wedges, _ = ax.pie([totals[category] for category in CLASSES], colors=COLORS,
                       radius=1, startangle=90, counterclock=True, normalize=True,
                       wedgeprops={'width': 0.43, 'edgecolor': 'white', 'linewidth': 0.6})
    for category, wedge in zip(CLASSES, wedges):
        theta = radians((wedge.theta1 + wedge.theta2) / 2)
        label = f'{100*totals[category]/n:.1f}%'
        if 100*totals[category]/n < 5:

            ax.annotate(label, xy=(0.59*cos(theta), 0.59*sin(theta)),
                        xytext=(0.30*cos(theta), 0.30*sin(theta)),
                        ha='center', va='center', fontsize=7, color=INK,
                        linespacing=1.3,
                        arrowprops={'arrowstyle': '-', 'color': GRAY,
                                    'linewidth': 0.6, 'shrinkA': 3, 'shrinkB': 0})
        else:
            ax.text(0.79*cos(theta), 0.79*sin(theta), label,
                    ha='center', va='center', fontsize=7, color=INK, linespacing=1.3)
    ax.set_xlim(-1.06, 1.06)
    ax.set_ylim(-1.06, 1.06)
    ax.set_aspect('equal')
    ax.axis('off')
    for category, color, y in zip(CLASSES, COLORS, (0.27, 0.18, 0.09)):
        fig.add_artist(Rectangle((0.04, y-0.014), 0.024, 0.028,
                                 transform=fig.transFigure, facecolor=color, edgecolor='none'))
        fig.text(0.085, y, category, fontsize=7, color=INK, ha='left', va='center')
    fig.text(2/WIDTH_MM, 1-2/HEIGHT_MM, 'D', fontsize=10, weight='bold', va='top', color=INK)
    fig.text(1-2/WIDTH_MM, 1-2/HEIGHT_MM, f'n = {n:,} loci',
             fontsize=7, ha='right', va='top', color=INK)
    for ext in ('pdf', 'svg', 'png'):
        fig.savefig(run/f'Figure1D_predicted_biotype_composition.{ext}', dpi=DPI,
                    facecolor='white', transparent=False, bbox_inches=None)

    import io
    from PIL import Image
    with io.BytesIO() as buffer:
        fig.savefig(buffer, format='png', dpi=DPI, facecolor='white',
                    transparent=False, bbox_inches=None)
        buffer.seek(0)
        with Image.open(buffer) as image:
            image.convert('RGB').save(run/'Figure1D_predicted_biotype_composition.tiff',
                                      format='TIFF', compression='tiff_lzw', dpi=(DPI, DPI))
    plt.close(fig)
    (run/'Figure1D.biotype_caption.txt').write_text(
        f'D. FEELnc-predicted biotype composition of {n:,} candidate loci represented by '
        f'{len(records):,} retained candidate transcript models from the five-dog dataset. '
        'Transcript annotations were taken from FEELnc_biotype in the final transcript mapping '
        'and grouped by final_gene_id. Each locus is counted once. '
        'The upstream annotation maps FEELnc_codpot_category lncRNA and noORF to '
        'FEELnc_biotype lncRNA, and mRNA to mRNA. This mapping is retained here without '
        'rerunning coding-potential prediction. '
        f'The lncRNA-labeled transcript set includes {noorf_transcripts:,} noORF transcript(s) '
        f'across {noorf_loci:,} locus/loci; these are not an additional category. '
        'lncRNA-only denotes loci for which all retained candidate transcripts are annotated '
        'as lncRNA by FEELnc; mRNA-only denotes loci for which all retained candidate transcripts '
        'are annotated as mRNA; Mixed denotes loci containing both predicted transcript types. '
        'These labels describe predictions for retained candidate models, not all possible isoforms '
        'at a genomic locus, experimentally verified coding status, or dual function of one transcript. '
        + ' '.join(f'{category}: {totals[category]:,}/{n:,} ({100*totals[category]/n:.1f}%).'
                   for category in CLASSES)
        + ' Donut sector angles show locus proportions. Labels give percentages of candidate loci. Small sectors (<5%) use an inward leader to a label in the donut hole; the total number of loci is displayed in the panel header. The complete final candidate set is used without '
        'filtering on 5′/3′ endpoint reference assessability. Biotype composition does not establish '
        'transcript completeness or endpoint accuracy. Counts are descriptive, without significance testing.\n',
        encoding='utf-8')
    (run/'inputs_and_parameters.txt').write_text(
        f'map={MAP}\noutput={run}\npython={sys.executable}\nmatplotlib={matplotlib.__version__}\n'
        f'total_transcripts={len(records)}\ntotal_loci={n}\n'
        f'noORF_transcripts={noorf_transcripts}\nnoORF_affected_loci={noorf_loci}\n'
        'noORF_policy=retain_upstream_lncRNA_biotype\n'
        'panel=predicted_biotype_composition\nunit=gene_locus\ngroup_key=final_gene_id\n'
        'biotype_field=FEELnc_biotype\nendpoint_assessability_filter=none\n'
        'chart=donut_with_percentages_bottom_left_legend\nwidth_cm=8.7\nheight_cm=6\npng_dpi=600\nconcurrent_tasks=1\nthreads_per_task=1\n',
        encoding='utf-8')

    with (run/'inputs_and_parameters.txt').open('a', encoding='utf-8') as handle:
        handle.write(f'font_family={FONT}\nfont_normal={fonts["normal"]}\nfont_bold={fonts["bold"]}\n')
        handle.write('panel_font_pt=10\nbody_font_pt=7\npanel_label_left_mm=2\npanel_label_top_mm=2\n')
        handle.write(f'tiff_dpi={DPI}\ntiff_color_mode=RGB\ntiff_compression=LZW\n')
        handle.write('intended_placement_width_mm=87\nrendered_collision_qc=NOT_CHECKED\n')
        handle.write('assembled_figure1_qc=NOT_CHECKED\ncaption_page_fit=NOT_CHECKED\n')
        handle.write(f'script={Path(__file__).resolve()}\n')
        for ext in ('pdf', 'svg', 'png', 'tiff'):
            handle.write(f'figure_{ext}={run / ("Figure1D_predicted_biotype_composition."+ext)}\n')


def main():
    run_id = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    run = ROOT/run_id
    run.mkdir(exist_ok=False)
    logdir = ROOT/'log'
    logdir.mkdir(exist_ok=True)

    with (logdir/f'{run_id}.master.log').open('x', buffering=1) as log, \
            (logdir/f'{run_id}.success.tsv').open('x') as success, \
            (logdir/f'{run_id}.failed.tsv').open('x') as failed, \
            redirect_stdout(log), redirect_stderr(log):
        logging.basicConfig(level=logging.INFO, stream=log, force=True,
                            format='[%(asctime)s ] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        success.write('task\tstatus\toutput_dir\n')
        failed.write('task\tstatus\treason\n')
        try:
            logging.info('[1/3] read final five dogs transcript mapping; default concurrent_tasks=1 threads_per_task=1')
            logging.info('map=%s output=%s', MAP, run)
            records, loci = read_records()
            os.environ['MPLCONFIGDIR'] = str(run/'.mplconfig')
            logging.info('[2/3] by site summary FEELnc type, output composition figure, source data and figure legend')
            make_outputs(run, records, loci)
            logging.info('[3/3] writing completed')
            (run/'execution_completion.marker').write_text('execution_status=success\n')
            success.write(f'Figure1D\tsuccess\t{run}\n')
            success.flush()
            logging.info('total=1 success=1 failed=0 skipped=0')
            logging.info('task completed / done execution_status=success')
            logging.info('please check separately figure spacing, output completeness, denominator and biological interpretation; successful execution does not establish output quality.')
            return 0
        except Exception as exc:
            logging.exception('execution_status=failed')
            reason = str(exc).replace('\t', ' ').replace('\n', ' ')
            failed.write(f'Figure1D\tfailed\t{reason}\n')
            failed.flush()
            logging.info('total=1 success=0 failed=1 skipped=0')
            return 1


if __name__ == '__main__':
    raise SystemExit(main())
