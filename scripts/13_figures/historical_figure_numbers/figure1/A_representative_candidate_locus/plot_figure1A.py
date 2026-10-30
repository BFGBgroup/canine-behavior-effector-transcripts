#!/path/to/software/venvs/figure1a-pygenometracks-3.9/bin/python

import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import argparse
import csv
import re
from collections import Counter
import json
import logging
from contextlib import redirect_stdout, redirect_stderr
from datetime import datetime
from importlib.metadata import version
from pathlib import Path

ROOT = Path('/path/to/repository/scripts/13_figures/historical_figure_numbers/figure1/A_representative_candidate_locus')
MANIFEST = Path('/path/to/project/legacy/five-dogs/9.new-gene/1.5_stringtie3_bulk_transcript_discovery-v6/stringtie3_input_manifest.tsv')
GTF = Path('/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.new_only.gtf')
REFERENCE = Path('/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.gtf')
DOGS = ['Ace_396226', 'CodRou_346627', 'Crak_358336', 'Dozer_386519', 'SanMcg_349891']
TISSUES = {'cortex_occipital': 'Occipital cortex', 'cortex_frontal': 'Frontal cortex',
           'cerebellum': 'Cerebellum', 'pituitary_gland': 'Pituitary'}
WIDTH_MM = 178.0
HEIGHT_MM = 55.0
FONT = 'Arial'

FONT_FILES = [Path('/path/to/fonts/Arial/Arial.TTF'),
              Path('/path/to/fonts/Arial/Arialbd.TTF')]
DPI = 600


def configure_font(font_files):

    from matplotlib import font_manager
    for path in font_files:
        if not path.is_absolute() or not path.is_file():
            raise ValueError('--font-file must point to existing font file of absolute path')
        font_manager.fontManager.addfont(str(path))
    resolved = {}
    try:
        for weight in ('normal', 'bold'):
            props = font_manager.FontProperties(family=FONT, weight=weight)
            path = font_manager.findfont(props, fallback_to_default=False)
            if font_manager.FontProperties(fname=path).get_name() != FONT:
                raise ValueError('matching to of font is not Arial')
            resolved[weight] = str(Path(path).resolve())
        if resolved['normal'] == resolved['bold']:
            raise ValueError('missing independent of Arial bold font')
    except ValueError as exc:
        raise RuntimeError(
            'Arial regular/bold unavailable.please installation obtained through permitted access of Arial, or duplicate use '
            '--font-file provide regular and bold file; do not automatically substitute another font.') from exc
    return resolved


def models(path, gene=None, region=None):

    result = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            f = line.rstrip().split('\t')
            if len(f) != 9 or f[2] not in ('transcript', 'exon'):
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', f[8]))
            if gene is not None and attrs.get('gene_id') != gene:
                continue
            if region and f[0] != region[0]:
                continue
            lo, hi = int(f[3]) - 1, int(f[4])
            tid = attrs['transcript_id']
            m = result.setdefault(tid, dict(gene=None, chrom=f[0], strand=f[6], exons=[]))
            if attrs.get('gene_id'):
                m['gene'] = attrs['gene_id']
            if f[2] == 'transcript':
                m['span'] = (lo, hi)
            else:
                m['exons'].append((lo, hi))
    for tid, m in result.items():
        if m['gene'] is None:
            raise ValueError(f'{path}: {tid} missing gene_id')
        m['exons'].sort()
        if 'span' not in m:
            m['span'] = (m['exons'][0][0], m['exons'][-1][1])
    if region:
        result = {tid: m for tid, m in result.items()
                  if m['span'][0] < region[2] and m['span'][1] > region[1]}
    return result


def extract(sample, chrom, start, end, strand):

    import numpy as np
    import pysam
    delta = np.zeros(end - start + 1, dtype=np.int64)
    junctions = Counter()
    with pysam.AlignmentFile(sample['bam_path'], 'rb', threads=1) as bam:
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
                continue
            if not read.is_paired or not (read.is_read1 or read.is_read2):
                raise ValueError('rf library encountered cannot determine orientation of read')
            if read.get_tag('NH') != 1:
                continue
            plus = read.is_reverse if read.is_read1 else not read.is_reverse
            if plus != (strand == '+'):
                continue
            pos = read.reference_start
            for op, length in read.cigartuples:
                if op in (0, 7, 8):
                    lo, hi = max(start, pos), min(end, pos + length)
                    if lo < hi:
                        delta[lo - start] += 1
                        delta[hi - start] -= 1
                    pos += length
                elif op == 3:
                    if start <= pos < pos + length <= end:
                        junctions[(pos, pos + length)] += 1
                    pos += length
                elif op == 2:
                    pos += length
    return np.cumsum(delta)[:-1], junctions


def write_bed12(path, group, gene_prefix=None):

    with path.open('w') as out:
        for tid, m in sorted(group.items(), key=lambda item: (item[1]['chrom'], item[1]['span'], item[0])):
            lo, hi = m['span']
            label = tid.removeprefix(gene_prefix + '.') if gene_prefix else tid
            sizes = ','.join(str(b-a) for a, b in m['exons'])
            starts = ','.join(str(a-lo) for a, _ in m['exons'])
            fields = [m['chrom'], lo, hi, label, 0, m['strand'], lo, lo, 0,
                      len(m['exons']), sizes, starts]
            out.write('\t'.join(map(str, fields)) + '\n')


def write_tracks_config(path, samples, run, gene, ceiling, chosen, refs):

    sections = ['# Generated from BAM and GTF; 0-based half-open inputs.',
                '[spacer]', 'height = 0.40', '']
    for sample in samples:
        sections.extend([f"[coverage_{sample['sample']}]",
                         f"file = {run / (gene+'.'+sample['sample_tissue']+'.bedgraph')}",
                         'file_type = bedgraph', 'height = 0.65',
                         f"title = {sample['sample'].split('_')[0]}",
                         'color = #3B73B9', 'min_value = 0', f'max_value = {ceiling}',
                         'type = fill', 'show_data_range = true',
                         'use_middle = false', 'rasterize = false',
                         '', '[spacer]', 'height = 0.08', ''])
    for name, group, color in [('reference', refs, '#90979F'), ('candidate', chosen, '#4D4D4D')]:
        if not group:
            continue
        sections.extend([f'[{name}]', f'file = {run / (gene+"."+name+".bed")}',
                         'file_type = bed', f'title = {name.capitalize()}',
                         f'height = {0.20*len(group):.2f}', 'style = flybase',
                         f'labels = {"false" if name == "candidate" else "true"}',
                         'fontsize = 7', 'arrowhead_included = true',
                         f'color = {color}', f'color_utr = {color}',
                         f'color_backbone = {color}', f'color_arrow = {color}',
                         'border_color = none', 'height_utr = 1', 'line_width = 0.5',
                         'global_max_row = true', '', '[spacer]', 'height = 0.08', ''])
    sections.extend(['[x-axis]', 'fontsize = 7', 'where = bottom', ''])
    path.write_text('\n'.join(sections), encoding='utf-8')


def run_plot(args, run):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from pygenometracks import tracksClass
    from PIL import Image
    import io

    import numpy as np
    logging.info('[1/4] check font with original input; concurrent_tasks=1 threads_per_task=1')
    logging.info('output_dir=%s', run)
    logging.info('manifest=%s；candidate_gtf=%s；reference_gtf=%s', MANIFEST, GTF, REFERENCE)
    fonts = configure_font(FONT_FILES if args.font_file is None else args.font_file)
    logging.info('font_family=%s；normal=%s；bold=%s', FONT, fonts['normal'], fonts['bold'])
    with MANIFEST.open() as handle:
        samples = [r for r in csv.DictReader(handle, delimiter='\t') if r['tissue'] == args.tissue]
    expected = set(DOGS) - ({'Crak_358336'} if args.tissue == 'pituitary_gland' else set())
    if len(samples) != len(expected) or {r['sample'] for r in samples} != expected:
        raise ValueError('sample list does not satisfy five dogs tissue design')
    if any(r['strandness'] != 'rf' for r in samples):
        raise ValueError('list is not rf library')
    samples.sort(key=lambda r: DOGS.index(r['sample']))
    gene = args.gene
    chosen = models(GTF, gene=gene)
    if args.transcript:
        chosen = {tid: m for tid, m in chosen.items() if tid == args.transcript}
    if not chosen or any(not m['exons'] for m in chosen.values()):
        raise ValueError('candidate does not exist or no exon')
    identities = {(m['chrom'], m['strand']) for m in chosen.values()}
    if len(identities) != 1 or next(iter(identities))[1] not in ('+', '-'):
        raise ValueError('candidate locus of chromosome/strand mismatch')
    chrom, strand = next(iter(identities))
    start = max(0, min(m['exons'][0][0] for m in chosen.values()) - args.flank)
    end = max(m['exons'][-1][1] for m in chosen.values()) + args.flank
    refs = models(REFERENCE, region=(chrom, start, end))
    tracks = []
    track_inputs = []
    for sample in samples:
        logging.info('[2/4] extract %s %s:%s-%s; bam=%s',
                     sample['sample_tissue'], chrom, start+1, end, sample['bam_path'])
        depth, junctions = extract(sample, chrom, start, end, strand)
        tracks.append((depth, junctions))
        bounds = np.r_[0, np.flatnonzero(depth[1:] != depth[:-1])+1, len(depth)]
        bg_path = run / f"{gene}.{sample['sample_tissue']}.bedgraph"
        track_inputs.append(str(bg_path))
        with bg_path.open('w') as bg:
            for a, b in zip(bounds[:-1], bounds[1:]):
                bg.write(f'{chrom}\t{start+a}\t{start+b}\t{depth[a]}\n')
    with (run / f'{gene}.junctions.tsv').open('w') as out:
        out.write('sample_tissue\tchrom\tintron_start0\tintron_end0\tstrand\tunique_alignment_count\tmeets_support_threshold\n')
        for sample, (_, junctions) in zip(samples, tracks):
            for (a, b), count in sorted(junctions.items()):
                out.write(f"{sample['sample_tissue']}\t{chrom}\t{a}\t{b}\t{strand}\t{count}\t{int(count >= args.min_junction_reads)}\n")
    with (run / f'{gene}.models.tsv').open('w') as out:
        out.write('track\tgene_id\ttranscript_id\tchrom\tstrand\tfeature\tstart0\tend0\n')
        for name, group in [('reference', refs), ('candidate', chosen)]:
            for tid, m in sorted(group.items()):
                for feature, a, b in [('transcript', *m['span'])]+[('exon', a, b) for a, b in m['exons']]:
                    out.write(f"{name}\t{m['gene']}\t{tid}\t{m['chrom']}\t{m['strand']}\t{feature}\t{a}\t{b}\n")
            if group:
                bed_path = run / f'{gene}.{name}.bed'
                write_bed12(bed_path, group, gene_prefix=gene if name == 'candidate' else None)
                track_inputs.append(str(bed_path))
    ceiling = max(1, max(int(d.max()) for d, _ in tracks))
    logging.info('shared_depth_range=0–%s；transcripts=%s；reference_transcripts=%s',
                 ceiling, ','.join(chosen), len(refs))
    provenance = dict(gene=gene, transcripts=list(chosen), tissue=args.tissue,
                      flank=args.flank, min_junction_reads=args.min_junction_reads,
                      manifest=str(MANIFEST), candidate_gtf=str(GTF), reference_gtf=str(REFERENCE),
                      samples=samples, coordinate_system='0-based half-open', region=[chrom, start, end],
                      reference_transcript_count=len(refs), shared_depth_ceiling=ceiling,
                      renderer_source=str(Path(tracksClass.__file__).resolve()))

    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': [FONT],
                         'font.size': 7, 'axes.labelsize': 7.5,
                         'xtick.labelsize': 7, 'ytick.labelsize': 7,
                         'axes.linewidth': 0.6, 'lines.linewidth': 0.6,
                         'pdf.fonttype': 42, 'ps.fonttype': 42,
                         'svg.fonttype': 'none', 'savefig.bbox': None})
    config = run / 'tracks.ini'
    write_tracks_config(config, samples, run, gene, ceiling, chosen, refs)
    logging.info('[3/4] pyGenomeTracks plot; region=%s:%s-%s; size=%s×%s mm',
                 chrom, start + 1, end, WIDTH_MM, HEIGHT_MM)
    renderer = tracksClass.PlotTracks(
        str(config), fig_width=WIDTH_MM / 10, fig_height=HEIGHT_MM / 10,
        fontsize=7, dpi=DPI, track_label_width=0.12,
        plot_regions=[(chrom, start, end)])
    fig = None
    try:

        with io.BytesIO() as draft:
            fig = renderer.plot(draft, chrom, start, end, h_align_titles='left')
        fig.text(2 / WIDTH_MM, 0.985, 'A', ha='left', va='top',
                 fontsize=10, fontweight='bold')
        fig.text(0.052, 0.985, 'Read depth', ha='left', va='top', fontsize=7.5)

        fig.text(0.50, 0.985, gene, ha='center', va='top', fontsize=8)

        for ax in fig.axes:
            for label in ax.texts:
                if label.get_text() == chrom:
                    label.set_fontsize(7.5)
        logging.info('[4/4] export PDF/SVG/PNG and RGB, LZW compress TIFF')
        for ext in ('pdf', 'svg', 'png'):
            fig.savefig(run / f'{gene}.{ext}', dpi=DPI, facecolor='white',
                        transparent=False, bbox_inches=None)

        buffer = io.BytesIO()
        fig.savefig(buffer, format='png', dpi=DPI, facecolor='white',
                    transparent=False, bbox_inches=None)
        buffer.seek(0)
        with Image.open(buffer) as image:
            image.convert('RGB').save(run / f'{gene}.tiff', format='TIFF',
                                      compression='tiff_lzw', dpi=(DPI, DPI))
    finally:
        if fig is not None:
            plt.close(fig)
        renderer.close_files()

    reference_note = ('No reference transcripts or exons overlap this window in the specified reference GTF. '
                      if not refs else 'Grey models show reference transcripts overlapping the window. ')
    caption = (f'A. RNA-seq coverage and candidate transcript models at {gene} '
               f'in {TISSUES[args.tissue].lower()} samples from {len(samples)} BarkBase dogs. '
               f'Models: {", ".join(chosen)}. '
               f'The displayed interval is {chrom}:{start+1:,}–{end:,} '
               '(CanFam3.1; 1-based inclusive coordinates). '
               f'Blue tracks show candidate-strand raw read depth on a shared linear scale of 0–{ceiling}. '
               'Dark grey boxes and connecting lines indicate candidate exons and introns; '
               'arrows indicate transcriptional direction. Coverage is not normalized expression. '
               + reference_note + 'Coverage alone does not validate full-length transcript structure.\n')
    (run / f'{gene}.caption.txt').write_text(caption, encoding='utf-8')
    methods = (
        'Read-counting and plotting details\n'
        'Coverage and junctions use rf paired-end primary NH=1 alignments; unmapped, secondary, '
        'supplementary and QC-fail alignments are excluded. Candidate-strand reads are inferred '
        'as reverse read 1 or forward read 2 for a positive-strand model, and the opposite '
        'orientations for a negative-strand model. Duplicate flags are retained and overlapping '
        'mates are counted separately. No mapping-quality threshold, base-quality threshold '
        'or library-size normalization is applied. CIGAR M, = and X contribute to depth; N '
        'defines splice junctions; deletions advance the reference coordinate without adding coverage. '
        'Only junctions fully contained in the displayed interval are tabulated. '
        f'The meets_support_threshold column marks alignment counts >= {args.min_junction_reads} '
        'per sample without pooling dogs or removing lower counts. Junction arcs are not plotted.\n'
        'GTF 1-based inclusive coordinates are converted to 0-based half-open BED/bedGraph intervals. '
        'BED12 thickStart equals thickEnd: CDS coordinates are unspecified; exon blocks do not '
        'establish CDS/UTR boundaries or protein-coding status. Arrowheads end within terminal '
        'exon boundaries. Full identifiers are retained in models.tsv and the caption; the full '
        'candidate gene ID is shown in the header. Reference absence is limited to the supplied '
        'annotation and displayed interval.\n'
        'BedGraphTrack plots the original interval values without bin aggregation, interpolation, '
        'smoothing or rasterization. Coverage uses BarkBase blue (#3B73B9); candidate structures '
        'use dark grey (#4D4D4D). Font: Arial. '
        f'Rendered by pyGenomeTracks {version("pyGenomeTracks")} using tracks.ini.\n')
    (run / f'{gene}.methods.txt').write_text(methods, encoding='utf-8')
    provenance['layout_run'] = {
        'track_inputs': track_inputs,
        'source_models': str(run / f'{gene}.models.tsv'),
        'source_junctions': str(run / f'{gene}.junctions.tsv'),
        'output_dir': str(run), 'font_family': FONT, 'font_files': fonts,
        'width_mm': WIDTH_MM, 'height_mm': HEIGHT_MM, 'raster_dpi': DPI,
        'text_sizes_pt': {'panel_label': 10, 'candidate_label': 8,
                          'axis_title': 7.5, 'track_labels_and_ticks': 7},
        'intended_placement_width_mm': WIDTH_MM,
        'tiff_color_mode': 'RGB', 'tiff_compression': 'LZW',
        'coverage_source': 'computed from manifest BAMs in this run',
        'software_versions': {name: version(name) for name in
                              ('pyGenomeTracks', 'matplotlib', 'Pillow', 'numpy', 'pysam')},
        'qc_status': {'rendered_font_size_and_collision': 'NOT_CHECKED',
                      'assembled_figure1': 'NOT_CHECKED',
                      'manuscript_placement_and_caption_fit': 'NOT_CHECKED'},
    }
    (run / 'inputs_and_parameters.json').write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False), encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description='Run plot_figure1A.py; see the workflow documentation.')
    parser.add_argument('--gene', default='DOG5_NEWGENE_01279')
    parser.add_argument('--transcript', default=None,
                        help='Plot one transcript of the specified gene; omit to plot all retained models for that gene.')
    parser.add_argument('--tissue', choices=TISSUES, default='cortex_occipital')
    parser.add_argument('--flank', type=int, default=1000)
    parser.add_argument('--min-junction-reads', type=int, default=3)
    parser.add_argument('--font-file', type=Path, action='append', default=None,
                        help='Replace the default Arial file list; specify this option repeatedly with the absolute paths to regular and bold font files.')
    args = parser.parse_args()
    if args.flank < 0 or args.min_junction_reads < 1:
        parser.error('flank must be nonnegative; min-junction-reads must be positive number')
    run = ROOT / datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    run.mkdir(exist_ok=False)
    os.environ['MPLCONFIGDIR'] = str(run / '.mplconfig')
    with (run / 'figure1A.log').open('x', buffering=1) as log:
        with redirect_stdout(log), redirect_stderr(log):
            logging.basicConfig(stream=log, level=logging.INFO,
                                format='[%(asctime)s] %(message)s', force=True)
            try:
                run_plot(args, run)
                (run / 'execution_completion.marker').write_text('execution_status=success\n')
                logging.info('task completed / done; execution_status=success; please check separately font, '
                             'collision, track align, vector embedding and page space occupied by the complete figure and legend; not yet completed publication validation')
                return 0
            except Exception:
                logging.exception('execution_status=failed')
                return 1


if __name__ == '__main__':
    raise SystemExit(main())
