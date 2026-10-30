#!/path/to/software/miniconda3/envs/dogs/bin/python

import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'

import csv
import logging
import math
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

ROOT = Path('/path/to/repository/scripts/13_figures/historical_figure_numbers/figure1/C_transcript_architecture')
REFERENCE_ROOT = Path('/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722')
GTF = REFERENCE_ROOT / 'five_only_new_gene_loci.new_only.gtf'
MAP = REFERENCE_ROOT / 'five_only_new_gene_loci.transcript_map.tsv'
COLOR = '#3B73B9'
WIDTH_MM, HEIGHT_MM = 87.0, 72.0
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


def read_architecture():
    transcripts = {}
    with MAP.open(newline='', encoding='utf-8') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            tid = row['final_transcript_id']
            if row['source_dataset'] != 'five-dogs' or tid in transcripts:
                raise ValueError(f'transcript input duplicate or non-five-dogs: {tid}')
            transcripts[tid] = row['final_gene_id']
    exons = defaultdict(list)
    with GTF.open(encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9 or fields[2] != 'exon':
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', fields[8]))
            tid = attrs['transcript_id']
            if tid in transcripts:
                if attrs['gene_id'] != transcripts[tid]:
                    raise ValueError(f'Gene identifiers differ between the GTF and mapping table: {tid}')
                exons[tid].append((int(fields[3]), int(fields[4])))
    records = []
    for tid, gene in sorted(transcripts.items()):
        intervals = sorted(exons[tid])
        if not intervals:
            raise ValueError(f'transcript missing exon: {tid}')

        length = sum(end-start+1 for start, end in intervals)
        records.append((gene, tid, len(intervals), length))
    if not records:
        raise ValueError('input no candidate transcript')
    return records


def make_figure(run, records):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fonts = configure_font()
    logging.info('font_family=%s normal=%s bold=%s', FONT, fonts['normal'], fonts['bold'])
    matplotlib.rcParams.update({'font.family': 'sans-serif',
                                'font.sans-serif': [FONT], 'font.size': 7,
                                'pdf.fonttype': 42, 'ps.fonttype': 42,
                                'svg.fonttype': 'none', 'savefig.bbox': None})
    lengths = [r[3] for r in records]
    total = len(records)
    counts = Counter(min(r[2], 5) for r in records)
    median = float(np.median(lengths))

    lower, upper = min(200, min(lengths)), max(20000, max(lengths))
    edges = np.geomspace(lower, upper, 31)
    frequencies, edges = np.histogram(lengths, bins=edges)
    with (run / 'Figure1C.transcripts.tsv').open('x') as handle:
        handle.write('final_gene_id\tfinal_transcript_id\texon_count\tspliced_length_bp\n')
        for row in records:
            handle.write('\t'.join(map(str, row))+'\n')
    with (run / 'Figure1C.exon_counts.tsv').open('x') as handle:
        handle.write('exon_group\ttranscript_count\ttotal_transcripts\tpercent\n')
        for k in range(1, 6):
            label = str(k) if k < 5 else '>=5'
            handle.write(f'{label}\t{counts[k]}\t{total}\t{100*counts[k]/total:.6f}\n')
    with (run / 'Figure1C.length_bins.tsv').open('x') as handle:
        handle.write('left_bp\tright_bp\tright_inclusive\ttranscript_count\n')
        for i, n in enumerate(frequencies):
            handle.write(f'{edges[i]:.17g}\t{edges[i+1]:.17g}\t{int(i==29)}\t{n}\n')
    with (run / 'inputs_and_parameters.txt').open('x') as handle:
        handle.write(f'gtf={GTF}\nmap={MAP}\noutput={run}\npython={sys.executable}\npython_version={sys.version}\n')
        handle.write(f'numpy={np.__version__}\nmatplotlib={matplotlib.__version__}\n')
        handle.write(f'unit=transcript\nn={total}\nmedian_bp={median:g}\nlength_min_bp={min(lengths)}\nlength_max_bp={max(lengths)}\n')
        handle.write('width_cm=8.7\nheight_cm=7.2\nlength_bins=30\nbin_spacing=log10\nconcurrent_tasks=1\nthreads_per_task=1\n')
        handle.write(f'font_family={FONT}\nfont_normal={fonts["normal"]}\nfont_bold={fonts["bold"]}\n')
        handle.write(f'raster_dpi={DPI}\ntiff_color_mode=RGB\ntiff_compression=LZW\n')
        handle.write('panel_font_pt=10\nsubtitle_font_pt=8\naxis_title_pt=7.5\nbody_font_pt=7\n')
        handle.write('exon_title_pad_pt=10\nlength_title_pad_pt=4\nlength_axis_label=Length (bp; log10 scale)\n')
        handle.write('intended_placement_width_mm=87\nrendered_collision_qc=NOT_CHECKED\n')
        handle.write('assembled_figure1_qc=NOT_CHECKED\ncaption_page_fit=NOT_CHECKED\n')
        handle.write(f'script={Path(__file__).resolve()}\n')
        for ext in ('pdf', 'svg', 'png', 'tiff'):
            handle.write(f'figure_{ext}={run / ("Figure1C_transcript_architecture."+ext)}\n')

    fig = plt.figure(figsize=(WIDTH_MM/25.4, HEIGHT_MM/25.4))
    top = fig.add_axes([18/WIDTH_MM, 43/HEIGHT_MM, 65/WIDTH_MM, 18/HEIGHT_MM])
    bottom = fig.add_axes([18/WIDTH_MM, 10/HEIGHT_MM, 65/WIDTH_MM, 18/HEIGHT_MM])
    for ax in (top, bottom):
        ax.spines[['top', 'right']].set_visible(False)
        for side in ('left', 'bottom'):
            ax.spines[side].set_linewidth(0.6)
        ax.tick_params(labelsize=7, width=0.6, length=2)
        ax.set_ylabel('Number of transcripts', fontsize=7.5, labelpad=4)
    values = [counts[k] for k in range(1, 6)]
    top.bar(range(5), values, width=0.65, color=COLOR)
    top.set_xticks(range(5), ['1', '2', '3', '4', '≥5'])
    top.set_ylim(0, math.ceil(max(values)*1.10/500)*500)
    for i, n in enumerate(values):
        top.text(i, n+max(values)*0.06, f'{n:,}', ha='center', va='bottom', fontsize=7)
    top.set_xlabel('Number of exons', fontsize=7.5, labelpad=3)


    top.set_title('Exon number distribution', fontsize=8, loc='left', pad=10)
    bottom.bar(edges[:-1], frequencies, width=np.diff(edges), align='edge', color=COLOR,
               edgecolor='white', linewidth=0.3)
    bottom.set_xscale('log')

    from matplotlib.ticker import FuncFormatter, NullFormatter
    bottom.xaxis.set_major_formatter(FuncFormatter(lambda value, position: f'{value:,.0f}'))
    bottom.xaxis.set_minor_formatter(NullFormatter())
    bottom.set_xlim(lower, upper)
    bottom.set_ylim(0, max(frequencies)*1.25)
    bottom.axvline(median, color='#637080', lw=0.8, ls='--')
    bottom.text(0.98, 0.96, f'Median = {median:,.0f} bp', transform=bottom.transAxes,
                ha='right', va='top', fontsize=7)


    bottom.set_xlabel('Length (bp; log10 scale)', fontsize=7.5, labelpad=3)
    bottom.set_title('Spliced transcript length', fontsize=8, loc='left', pad=4)

    fig.text(2/WIDTH_MM, 1-2/HEIGHT_MM, 'C', fontsize=10, weight='bold', va='top')
    fig.text(1-2/WIDTH_MM, 1-2/HEIGHT_MM, f'n = {total:,} transcripts',
             ha='right', va='top', fontsize=7)

    sys.path.insert(0, '/path/to/repository/support/figure_alignment')
    from audit_panel_alignment import require_matplotlib_panel_alignment
    require_matplotlib_panel_alignment(
        fig, axes=[top, bottom], panel_ids=['exon_count', 'spliced_length'],
        column_groups=[['exon_count', 'spliced_length']],
        json_out=run / 'Figure1C.alignment.json', strict=True)
    for ext in ('pdf', 'svg', 'png'):
        fig.savefig(run / f'Figure1C_transcript_architecture.{ext}', dpi=DPI,
                    facecolor='white', transparent=False, bbox_inches=None)

    import io
    from PIL import Image
    with io.BytesIO() as buffer:
        fig.savefig(buffer, format='png', dpi=DPI, facecolor='white',
                    transparent=False, bbox_inches=None)
        buffer.seek(0)
        with Image.open(buffer) as image:
            image.convert('RGB').save(run / 'Figure1C_transcript_architecture.tiff',
                                      format='TIFF', compression='tiff_lzw', dpi=(DPI, DPI))
    plt.close(fig)
    with (run / 'Figure1C.caption.txt').open('x') as handle:
        handle.write(f'C. Architecture of {total:,} retained BarkBase candidate transcripts. Each final transcript ID is counted once; multiple isoforms from the same locus contribute separately.\n')
        handle.write('Top: exon counts per transcript, grouped as 1, 2, 3, 4, and at least 5 exons. Bar labels give transcript counts.\n')
        handle.write('The length x-axis uses a base-10 logarithmic scale; the y-axis shows untransformed transcript counts per bin. The median is reported in base pairs, not log-transformed units.\n')
        handle.write(f'Bottom: spliced lengths calculated as the sum of exon lengths from the final CanFam3.1 candidate GTF, excluding introns. Histogram counts use 30 equal-width log10 bins from {lower:g} to {upper:g} bp, not equal-width base-pair bins or density normalization. Bins are left-inclusive and right-exclusive except the final bin, which includes both edges.\n')
        handle.write(f'The dashed line indicates the median ({median:,.0f} bp). These distributions describe reconstructed models and do not establish transcript completeness or coding function.\n')


def main():
    run_id = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    run = ROOT / run_id
    run.mkdir(exist_ok=False)
    logdir = ROOT / 'log'
    logdir.mkdir(exist_ok=True)
    with (logdir / f'{run_id}.master.log').open('x', buffering=1) as log, \
            (logdir / f'{run_id}.success.tsv').open('x') as success, \
            (logdir / f'{run_id}.failed.tsv').open('x') as failed:
        sys.stdout = sys.stderr = log
        logging.basicConfig(level=logging.INFO, stream=log, force=True,
                            format='[%(asctime)s ] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        success.write('task\tstatus\toutput_dir\n')
        failed.write('task\tstatus\treason\n')
        try:
            logging.info('[1/3] read mapping table and GTF; concurrent_tasks=1 threads_per_task=1(default)')
            logging.info('gtf=%s map=%s output=%s', GTF, MAP, run)
            records = read_architecture()
            logging.info('[2/3] compute distribution and plot; transcripts=%d', len(records))
            os.environ['MPLCONFIGDIR'] = str(run / '.mplconfig')
            make_figure(run, records)
            logging.info('[3/3] figure, source data and figure legends have been written')
            success.write(f'Figure1C\tsuccess\t{run}\n')
            success.flush()
            (run / 'execution_completion.marker').write_text('execution_status=success\n')
            logging.info('total=1 success=1 failed=0 skipped=0')
            logging.info('task completed / done execution_status=success')
            logging.info('please separately manual check output completeness, text spacing, binning and biological interpretation; successful execution does not establish output quality.')
            return 0
        except Exception as exc:
            logging.exception('execution_status=failed')
            reason = str(exc).replace('\t', ' ').replace('\n', ' ')
            failed.write(f'Figure1C\tfailed\t{reason}\n')
            logging.info('total=1 success=0 failed=1 skipped=0')
            return 1


if __name__ == '__main__':
    raise SystemExit(main())
