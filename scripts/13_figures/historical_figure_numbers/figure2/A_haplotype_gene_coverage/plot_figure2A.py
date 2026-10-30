#!/path/to/software/miniconda3/envs/dogs/bin/python

import os

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[name] = '1'

import csv
from io import BytesIO
import logging
import json
import math
import sys
from contextlib import redirect_stdout, redirect_stderr
from datetime import datetime
from pathlib import Path
from statistics import median

ROOT = Path('/path/to/repository/scripts/13_figures/historical_figure_numbers/figure2/A_haplotype_gene_coverage')
MANIFEST = Path('/path/to/repository/scripts/08_expression_ase_integration/preparation/1.5ase-cnv-clean-gene-id-haplotype-table/log/1.5.ase_cnv_clean_gene_id_haplotype_table.20260821_161805.success.tsv')
STEM = 'Figure2A_haplotype_gene_coverage'
DOGS = ['Ace_396226', 'CodRou_346627', 'Crak_358336', 'Dozer_386519', 'SanMcg_349891']
TISSUES = ['cerebellum', 'cortex_frontal', 'cortex_occipital', 'pituitary_gland']
LABELS = ['Cerebellum', 'Frontal\ncortex', 'Occipital\ncortex', 'Pituitary']
COLORS = ['#D55E00', '#0072B2', '#009E73', '#E69F00', '#CC79A7']
MARKERS = ['o', 's', '^', 'D', 'v']
COUNT_UNIT = 'ASEReadCounter_site_allele_support_sum'


WIDTH_MM = 87.0
HEIGHT_MM = 75.0
CONTENT_HEIGHT_MM = 65.0
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


def load_samples():

    with MANIFEST.open(encoding='utf-8', newline='') as handle:
        entries = list(csv.DictReader(handle, delimiter='\t'))
    expected = {(dog, tissue) for dog in DOGS for tissue in TISSUES}
    expected.remove(('Crak_358336', 'pituitary_gland'))
    keys = [(entry['animal_id'], entry['tissue']) for entry in entries]
    if len(keys) != len(set(keys)) or set(keys) != expected:
        raise ValueError('Success list does not match the current 19-sample design or contains duplicates; check the input.')
    return sorted(entries, key=lambda entry: (
        TISSUES.index(entry['tissue']), DOGS.index(entry['animal_id'])))


def count_genes(entry):

    path = Path(entry['output_tsv'])
    qc_path = Path(entry['qc_tsv'])
    if not path.is_absolute() or not qc_path.is_absolute():
        raise ValueError('input paths in the success list must be absolute.')
    with qc_path.open(encoding='utf-8', newline='') as handle:
        qc = {row['metric']: row['value'] for row in csv.DictReader(handle, delimiter='\t')}
    chromosomes = {}
    selected = {'known_gene': set(), 'manual_new_gene': set()}
    input_rows = {'known_gene': 0, 'manual_new_gene': 0}
    unannotated_rows = 0
    with path.open(encoding='utf-8', newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['sample'] != entry['animal_id'] or row['tissue'] != entry['tissue']:
                raise ValueError(f'Sample identity does not match the list: {path}')
            if row['record_type'] == 'unannotated_variant':
                unannotated_rows += 1
                continue
            kind = row['record_type']
            if kind not in selected:
                raise ValueError(f"unknown record_type: {row['record_type']}")
            if row['count_unit'] != COUNT_UNIT:
                raise ValueError(f'counting units do not match the current counting scope: {path}')
            gene = row['gene_id']
            chromosome = row['chromosome']
            if not gene or gene.startswith('(NA)'):
                raise ValueError(f'annotation gene record use invalid ID: {path}')
            if gene in chromosomes and chromosomes[gene] != chromosome:
                raise ValueError(f'the same gene_id across chromosome: {gene}, please first verify annotation.')
            chromosomes[gene] = chromosome
            input_rows[kind] += 1

            if int(row['variant_count']) >= 1 and int(row['total_count_sum']) > 0:
                selected[kind].add(gene)
    return {
        'sample': entry['animal_id'], 'tissue': entry['tissue'],

        'n_genes': len(selected['known_gene'] | selected['manual_new_gene']),
        'n_known_genes': len(selected['known_gene']),
        'n_manual_new_genes': len(selected['manual_new_gene']),
        'known_gene_input_rows': input_rows['known_gene'],
        'manual_new_gene_input_rows': input_rows['manual_new_gene'],
        'excluded_unannotated_variant_rows': unannotated_rows,
        'upstream_known_gene_groups': int(qc['known_gene_groups']),
        'upstream_manual_new_gene_groups': int(qc['manual_new_gene_groups']),
        'cnv_input_policy': qc['cnv_input_policy'],
        'count_rule': 'unique_gene_id;known_gene_OR_manual_new_gene;variant_count>=1;total_count_sum>0',
        'count_unit': COUNT_UNIT, 'input_tsv': str(path), 'input_qc_tsv': str(qc_path),
    }


def draw(rows, run_dir):

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import StrMethodFormatter
    from PIL import Image
    sys.path.insert(0, '/path/to/repository/support/figure_alignment')
    from audit_panel_alignment import require_matplotlib_panel_alignment

    font_paths = configure_font()
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': [FONT],
        'font.size': 7, 'axes.labelsize': 7.5, 'xtick.labelsize': 7,
        'ytick.labelsize': 7, 'axes.linewidth': 0.7,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
        'figure.facecolor': 'white', 'savefig.facecolor': 'white',
        'savefig.bbox': None, 'savefig.transparent': False,
    })
    fig, ax = plt.subplots(figsize=(WIDTH_MM/25.4, HEIGHT_MM/25.4))

    def content_y(fraction):
        return 1 - CONTENT_HEIGHT_MM * (1 - fraction) / HEIGHT_MM

    fig.subplots_adjust(left=0.19, right=0.98,
                        bottom=content_y(0.25), top=content_y(0.78))
    fig.text(2/WIDTH_MM, 1-2/HEIGHT_MM, 'A', fontsize=10, fontweight='bold', va='top')
    for dog_index, dog in enumerate(DOGS):
        values = [row for row in rows if row['sample'] == dog]
        ax.scatter([TISSUES.index(row['tissue']) + (dog_index - 2) * 0.065
                    for row in values], [row['n_genes'] for row in values],
                   s=22, color=COLORS[dog_index], marker=MARKERS[dog_index],
                   linewidths=0.4, edgecolors='white', zorder=3,
                   label=dog.split('_')[0])
    for x, tissue in enumerate(TISSUES):
        counts = [row['n_genes'] for row in rows if row['tissue'] == tissue]

        ax.hlines(median(counts), x - 0.22, x + 0.22,
                  colors='#202020', linewidth=1.0, zorder=2)


        label_x = 0.19 + (0.98 - 0.19) * (x + 0.5) / len(TISSUES)
        fig.text(label_x, content_y(0.065), f'n = {len(counts)}',
                 ha='center', va='bottom', fontsize=7)
    ax.set_xticks(range(len(TISSUES)), LABELS)
    ax.set_xlim(-0.5, 3.5)

    upper = max(2000, math.ceil(max(row['n_genes'] for row in rows) * 1.04 / 2000) * 2000)
    ax.set_ylim(0, upper)
    ax.set_yticks(range(0, upper + 1, 4000 if upper > 12000 else 2000))
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    ax.set_ylabel('Number of genes')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', width=0.7, length=3)
    fig.legend(*ax.get_legend_handles_labels(), loc='upper center',
               bbox_to_anchor=(0.53, content_y(0.90)), ncol=5, frameon=False,
               fontsize=7, handletextpad=0.15, columnspacing=0.5)
    require_matplotlib_panel_alignment(
        fig,
        json_out=run_dir / f'{STEM}.alignment.json',
        tolerance_pt=1.5,
        strict=True,
    )
    fig.savefig(run_dir / f'{STEM}.pdf')
    fig.savefig(run_dir / f'{STEM}.svg')
    fig.savefig(run_dir / f'{STEM}.png', dpi=DPI)

    with BytesIO() as buffer:
        fig.savefig(buffer, format='png', dpi=DPI, facecolor='white')
        buffer.seek(0)
        with Image.open(buffer) as raster:
            with raster.convert('RGB') as rgb:
                rgb.save(run_dir / f'{STEM}.tiff', compression='tiff_lzw', dpi=(DPI, DPI))
    plt.close(fig)
    return font_paths


def main(run_dir):
    logging.info('start; concurrent_tasks=1; threads_per_task=1(new script default)')
    logging.info('input list: %s; output directory: %s', MANIFEST, run_dir)
    entries = load_samples()
    rows = []
    for index, entry in enumerate(entries, start=1):
        result = count_genes(entry)
        rows.append(result)
        logging.info('read %s/%s: %s / %s; n_genes=%s', index, len(entries),
                     result['sample'], result['tissue'], result['n_genes'])
    logging.info('save statistics source table with figure legend')
    with (run_dir / f'{STEM}.source_data.tsv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    caption = (
        'A. Genes with haplotype-resolved counts. Each point represents one dog–tissue '
        'sample (19 samples from five dogs; Crak has no pituitary sample). Colors and '
        'shapes identify dogs; black bars indicate tissue medians. Within each sample, '
        'unique gene IDs from reference genes (known_gene) and candidate gene-loci '
        '(manual_new_gene), with at least one variant and positive total '
        'allele-support counts were counted after the upstream phasing and specified '
        'CNV filtering procedure (known_cnv_removed: removal of single-state consensus '
        'CNV at 7 of 9 scales; flagged records may remain). Inputs are the 1.5 '
        'gene-level tables upstream of Allelome, without its subsequent bias selection. '
        'Unannotated variant records were excluded. Either '
        'haplotype may have zero counts. No additional coverage threshold or '
        'significance test was applied. Counts summarize SNP-level allele support, '
        'not independent RNA molecules. This panel describes data availability, '
        'not the number of significantly imbalanced genes.\n\n'
        f'Input manifest: {MANIFEST}\n\n'
        'Counts are deduplicated by gene_id and combine reference genes and candidate gene-loci; the source table separately retains the eligible counts for both categories, '
        'along with the upstream QC count of gene+chromosome groups and the CNV filtering policy.'
        'manual_new_gene is a candidate annotation category, not a validated new gene. '
        'Unannotated sites are not treated as candidate new genes, and passing the current CNV policy is not interpreted as proof of CNV absence.\n'
    )
    (run_dir / f'{STEM}.caption.md').write_text(caption, encoding='utf-8')
    logging.info('plot scatter points, medians and save PDF / SVG / PNG / RGB TIFF(600 dpi, LZW)')
    font_paths = draw(rows, run_dir)

    layout = dict(width_mm=WIDTH_MM, height_mm=HEIGHT_MM, raster_dpi=DPI,
                  content_height_mm=CONTENT_HEIGHT_MM,
                  bottom_padding_added_mm=HEIGHT_MM-CONTENT_HEIGHT_MM,
                  layout_role="first_row_75mm_trial_top_aligned",
                  plot_area_mm=dict(left=WIDTH_MM*0.19, top=CONTENT_HEIGHT_MM*0.22,
                                    width=WIDTH_MM*0.79, height=CONTENT_HEIGHT_MM*0.53),
                  font_family=FONT, resolved_font_files=font_paths,
                  font_sizes_pt=dict(panel=10, axis_title=7.5, tick=7, legend=7),
                  inputs=dict(manifest=str(MANIFEST),
                              gene_tables=[e['output_tsv'] for e in entries],
                              qc_tables=[e['qc_tsv'] for e in entries]),
                  panel_anchor_mm=dict(left=2, top=2),
                  review=dict(pdf_font_embedding='NOT_CHECKED', collisions='NOT_CHECKED',
                              figure2_assembly='NOT_CHECKED', figure_and_caption='NOT_CHECKED'))
    (run_dir / 'inputs_and_layout.json').write_text(json.dumps(layout, indent=2)+'\n', encoding='utf-8')
    (run_dir / 'execution_completion.marker').write_text('execution_status=success\n', encoding='utf-8')
    logging.info('task completed / done; execution_status=success')
    logging.info('please check separately source table, plot layout, font and biological interpretation; execution completed does not imply result quality passed.')


if __name__ == '__main__':

    run = ROOT / datetime.now().strftime('%Y%m%d_%H%M%S')
    run.mkdir(parents=True, exist_ok=False)
    os.environ['MPLCONFIGDIR'] = str(run / '.matplotlib')
    with (run / 'plot_figure2A.log').open('x', encoding='utf-8', buffering=1) as log_handle:
        with redirect_stdout(log_handle), redirect_stderr(log_handle):
            logging.basicConfig(stream=log_handle, level=logging.INFO,
                                format='[%(asctime)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
            try:
                main(run)
            except Exception:
                logging.exception('execution_status=failed; not completed, please check this log.')
                sys.exit(1)
