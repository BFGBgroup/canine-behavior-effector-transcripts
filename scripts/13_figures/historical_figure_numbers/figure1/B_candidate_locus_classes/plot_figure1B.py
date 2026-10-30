#!/path/to/software/miniconda3/envs/dogs/bin/python

import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'

import csv
import logging
import math
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

ROOT = Path('/path/to/repository/scripts/13_figures/historical_figure_numbers/figure1/B_candidate_locus_classes')
INPUT = Path('/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.transcript_map.tsv')
CLASSES = [
    ('retained_no_known_gene_overlap', 'No reference-gene\noverlap'),
    ('retained_opposite_strand_complete_overlap', 'Opposite-strand\nnear-containment'),
    ('retained_opposite_strand_partial_overlap', 'Opposite-strand\npartial overlap'),
    ('retained_same_strand_partial_overlap', 'Same-strand\npartial overlap'),
]
GRAY, CANDIDATE = '#92928F', '#B96C5D'
BAR_COLOR = '#3B73B9'
FIGURE_CM = (8.7, 7.2)
TEXT_PT, PANEL_PT = 7, 10
AXIS_PT, DPI = 7.5, 600
LEGEND_CENTER_BOTTOM_MM = 5.0
LABEL_LINESPACING = 1.30
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
        raise ValueError('Arial regular and bold must resolve to different files')
    return resolved


def read_loci():

    loci = {}
    allowed = {key for key, _ in CLASSES}
    with INPUT.open(newline='', encoding='utf-8') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['source_dataset'] != 'five-dogs':
                raise ValueError('input contains non-five-dogs record')
            gene, reason = row['final_gene_id'], row['gene_locus_filter_reason']
            if not gene or reason not in allowed:
                raise ValueError(f'cannot classify of input: {gene} {reason}')
            if gene in loci and loci[gene] != reason:
                raise ValueError(f'the same locus classification conflict: {gene}')
            loci[gene] = reason
    if not loci:
        raise ValueError('input has no candidate loci')
    return loci


def draw(run, counts, total):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, FancyArrow
    fig = plt.figure(figsize=tuple(value/2.54 for value in FIGURE_CM))


    labels = fig.add_axes([2/87, 0.23, 28/87, 0.66])
    schematic = fig.add_axes([32/87, 0.23, 14/87, 0.66])
    bars = fig.add_axes([49/87, 0.23, 34/87, 0.66])
    for ax in (labels, schematic, bars):
        ax.set_ylim(3.5, -0.5)
        ax.set_yticks([])
    for ax in (labels, schematic):
        ax.set_xlim(0, 1)
        ax.set_axis_off()

    def span(a, b, y, color, direction):

        origin, length = (a, b-a) if direction == '>' else (b, a-b)
        schematic.add_patch(FancyArrow(
            origin, y, length, 0, width=0.065,
            head_width=0.12, head_length=0.065,
            length_includes_head=True, linewidth=0, facecolor=color))

    for i, (key, label) in enumerate(CLASSES):
        labels.text(0, i, label, fontsize=TEXT_PT, va='center', linespacing=LABEL_LINESPACING)

        if i == 0:
            span(0.04, 0.37, i, GRAY, '>')
            span(0.63, 0.96, i, CANDIDATE, '>')
            schematic.text(0.50, i, '//', fontsize=TEXT_PT, va='center', ha='center', color=GRAY)
        elif i == 1:
            span(0.08, 0.92, i-0.12, GRAY, '>')
            span(0.27, 0.73, i+0.12, CANDIDATE, '<')
        else:
            span(0.08, 0.64, i-0.12, GRAY, '>')
            span(0.39, 0.95, i+0.12, CANDIDATE, '<' if i == 2 else '>')
        n = counts[key]
        bars.barh(i, n, height=0.40, color=BAR_COLOR)

        bars.annotate(f'{n:,}', xy=(n, i), xytext=(2, 0),
                      textcoords='offset points', fontsize=TEXT_PT,
                      ha='left', va='center', annotation_clip=False)

    axis_max = math.ceil(max(counts.values())*1.12/1000)*1000
    bars.set_xlim(0, axis_max)
    ticks = list(range(0, axis_max+1, 1000))
    bars.set_xticks(ticks, [f'{value:,}' for value in ticks])
    bars.spines[['top', 'right']].set_visible(False)
    bars.spines['left'].set_linewidth(0.6)
    bars.spines['left'].set_color(GRAY)
    bars.spines['bottom'].set_linewidth(0.6)
    bars.tick_params(axis='x', labelsize=TEXT_PT, length=2, width=0.6)
    bars.set_xlabel('Number of candidate loci', fontsize=AXIS_PT, labelpad=4)
    fig.text(2/87, 1-2/72, 'B', fontsize=PANEL_PT, weight='bold', va='top')

    fig.text(85/87, 1-2/72, f'n = {total:,} loci', fontsize=TEXT_PT,
             ha='right', va='top')


    legend_y = LEGEND_CENTER_BOTTOM_MM / (FIGURE_CM[1] * 10)
    for x, color, label in ((0.25, GRAY, 'Reference'), (0.55, CANDIDATE, 'Candidate')):
        schematic.add_patch(Rectangle((x, legend_y-0.011), 0.024, 0.022, transform=fig.transFigure,
                                      clip_on=False, color=color))
        fig.text(x+0.035, legend_y, label, fontsize=TEXT_PT, va='center')

    for ext in ('pdf', 'svg', 'png'):
        fig.savefig(run / f'Figure1B_candidate_locus_classes.{ext}', dpi=DPI,
                    facecolor='white', transparent=False, bbox_inches=None)

    import io
    from PIL import Image
    with io.BytesIO() as buffer:
        fig.savefig(buffer, format='png', dpi=DPI, facecolor='white',
                    transparent=False, bbox_inches=None)
        buffer.seek(0)
        with Image.open(buffer) as image:
            image.convert('RGB').save(run / 'Figure1B_candidate_locus_classes.tiff',
                                      format='TIFF', compression='tiff_lzw', dpi=(DPI, DPI))
    plt.close(fig)


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
        logging.basicConfig(level=logging.INFO, stream=log,
                            format='[%(asctime)s ] %(message)s', datefmt='%Y-%m-%d %H:%M:%S', force=True)
        success.write('task\tstatus\toutput_dir\n')
        failed.write('task\tstatus\treason\n')
        try:
            logging.info('[1/3] read final five dogs mapping table; concurrent_tasks=1 threads_per_task=1(default)')
            logging.info('input=%s output=%s', INPUT, run)
            logging.info('log=%s', logdir / f'{run_id}.master.log')
            loci = read_loci()
            counts, total = Counter(loci.values()), len(loci)
            logging.info('[2/3] by gene deduplicate and output source data; loci=%d', total)
            with (run / 'Figure1B.loci.tsv').open('x') as handle:
                handle.write('final_gene_id\tgene_locus_filter_reason\n')
                for gene, reason in sorted(loci.items()):
                    handle.write(f'{gene}\t{reason}\n')
            with (run / 'Figure1B.counts.tsv').open('x') as handle:
                handle.write('gene_locus_filter_reason\tlabel\tlocus_count\ttotal_loci\tpercent\n')
                for reason, label in CLASSES:
                    label = label.replace('\n', ' ')
                    handle.write(f'{reason}\t{label}\t{counts[reason]}\t{total}\t{100*counts[reason]/total:.6f}\n')
            os.environ['MPLCONFIGDIR'] = str(run / '.mplconfig')
            import matplotlib
            matplotlib.use('Agg')
            fonts = configure_font()
            logging.info('font_family=%s normal=%s bold=%s', FONT, fonts['normal'], fonts['bold'])
            matplotlib.rcParams.update({'font.family': 'sans-serif',
                                        'font.sans-serif': [FONT],
                                        'pdf.fonttype': 42, 'ps.fonttype': 42,
                                        'font.size': TEXT_PT, 'text.usetex': False,
                                        'svg.fonttype': 'none', 'savefig.bbox': None})
            with (run / 'inputs_and_parameters.txt').open('x') as handle:
                handle.write(f'input={INPUT}\noutput={run}\npython={sys.executable}\n')
                handle.write(f'python_version={sys.version}\nmatplotlib={matplotlib.__version__}\n')
                handle.write(f'total_loci={total}\nunit=final_gene_id\ncategory=gene_locus_filter_reason\nwidth_cm={FIGURE_CM[0]}\nheight_cm={FIGURE_CM[1]}\nlayout=category_schematic_count_columns\n')
                handle.write('concurrent_tasks=1\nthreads_per_task=1\nNo additional filtering by chromosome, tissue, ASE, or endpoint support.\n')
            with (run / 'inputs_and_parameters.txt').open('a') as handle:
                handle.write(f'body_font_pt={TEXT_PT}\npanel_font_pt={PANEL_PT}\nfont_family={FONT}\ncolor_mode=RGB\nprimary_format=vector_PDF\npreview_png_dpi={DPI}\ntiff_dpi={DPI}\ntiff_compression=LZW\naxis_title_pt={AXIS_PT}\n')
            with (run / 'inputs_and_parameters.txt').open('a') as handle:
                handle.write(f"font_normal={fonts['normal']}\nfont_bold={fonts['bold']}\n")
                handle.write('intended_placement_width_mm=87\nrendered_layout_qc=NOT_CHECKED\n')
                handle.write('assembled_figure1_qc=NOT_CHECKED\ncaption_page_fit=NOT_CHECKED\n')
                handle.write(f"script={Path(__file__).resolve()}\nlog={logdir / (run_id+'.master.log')}\n")
                handle.write('panel_label_left_mm=2\npanel_label_top_mm=2\ncount_label_offset_pt=2\n')
                handle.write(f'legend_center_bottom_mm={LEGEND_CENTER_BOTTOM_MM}\nlabel_linespacing={LABEL_LINESPACING}\n')
                for extension in ('pdf', 'svg', 'png', 'tiff'):
                    handle.write(f"figure_{extension}={run / ('Figure1B_candidate_locus_classes.'+extension)}\n")
            logging.info('[3/3] plot B figure and figure legend')
            draw(run, counts, total)
            with (run / 'Figure1B.caption.txt').open('x') as handle:
                handle.write(f'B. Reference-relative classes of {total:,} retained BarkBase candidate loci. Each final gene ID is counted once, regardless of transcript number.\n')
                handle.write('Bars and their labels show locus counts; percentages of all retained loci are provided in the source-data table. Gray and brick-red arrow bars denote reference gene spans and candidate locus spans, respectively, and are not exon models or drawn to scale. Reference models are oriented left to right; arrowheads on both spans indicate relative strand. The no-overlap example illustrates one possible orientation; this category does not require matching strands. The double slash indicates schematic separation, not a quantitative distance. All count bars use the same blue; color does not encode confidence or a category ranking.\n')
                handle.write('Classes use the upstream gene-span assignments relative to CanFam3.1 annotation. Near-containment means overlap covers at least 95% of either span; either the reference or candidate may contain the other (one arrangement is illustrated). It does not imply matching exon structures.\n')
                handle.write('Assignments are mutually exclusive using upstream priority: opposite-strand near-containment, same-strand partial overlap, opposite-strand partial overlap, then no reference-gene overlap. Display order differs from assignment priority.\n')
                handle.write('No reference-gene overlap is a locus-span category, not a synonym for transcript-level GffCompare class u. Retained candidates are not all independently validated new genes.\n')
            (run / 'execution_completion.marker').write_text('execution_status=success\n')
            success.write(f'Figure1B\tsuccess\t{run}\n')
            logging.info('total=1 success=1 failed=0 skipped=0')
            logging.info('task completed / done execution_status=success')
            logging.info('please another outside manual check figure completeness, text spacing, counting scope with biological interpretation; successful execution does not establish output quality.')
            return 0
        except Exception as exc:
            logging.exception('execution_status=failed')
            reason = str(exc).replace('\t', ' ').replace('\n', ' ')
            failed.write(f'Figure1B\tfailed\t{reason}\n')
            logging.info('total=1 success=0 failed=1 skipped=0')
            return 1


if __name__ == '__main__':
    raise SystemExit(main())
