#!/path/to/software/miniconda3/envs/dogs/bin/python

import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'

import csv
import logging
import json
import subprocess
import sys
from collections import Counter
from contextlib import redirect_stdout, redirect_stderr
from datetime import datetime
from pathlib import Path

ROOT = Path('/path/to/repository/scripts/13_figures/historical_figure_numbers/figure2/B_CNV_filtered_sites')
MANIFEST = Path('/path/to/repository/scripts/08_expression_ase_integration/preparation/1.4ase-cnv-qc-filter/log/1.4.ase_cnv_qc_filter.20260821_154914.success.tsv')
FAI = Path('/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna.fai')
STEM = 'Figure2B_CNV_filtered_sites'
DOGS = ['Ace_396226', 'CodRou_346627', 'Crak_358336', 'Dozer_386519', 'SanMcg_349891']
TISSUES = ['cerebellum', 'cortex_frontal', 'cortex_occipital', 'pituitary_gland']
SHORT = ['CB', 'FC', 'OC', 'Pit']
STATES = ['amp', 'del', 'cnloh']


COLORS = ['#D6604D', '#0072B2', '#009E73']
BIN_BP = 1_000_000
TISSUE_BACKGROUNDS = ['#F2F6F9', '#F7F2F8', '#F3F7F0', '#FAF5EC']


def read_inputs():

    with FAI.open() as handle:
        lengths_all = {parts[0]: int(parts[1]) for line in handle
                       if (parts := line.rstrip().split('\t'))}
    lengths = {f'chr{i}': lengths_all[f'chr{i}'] for i in range(1, 39)}
    with MANIFEST.open(encoding='utf-8', newline='') as handle:
        entries = list(csv.DictReader(handle, delimiter='\t'))
    expected = {(dog, tissue) for dog in DOGS for tissue in TISSUES}
    expected.remove(('Crak_358336', 'pituitary_gland'))

    keys = [(row['animal_id'], row['term']) for row in entries]
    if len(keys) != len(set(keys)) or set(keys) != expected:
        raise ValueError('Step 1.4 success list does not match the current 19-sample design.')
    entries.sort(key=lambda row: (TISSUES.index(row['term']), DOGS.index(row['animal_id'])))
    return lengths, entries


def aggregate(lengths, entries):

    rows, summaries = [], []
    for entry in entries:
        path = Path(entry['full_tsv'])
        if not path.is_absolute():
            raise ValueError(f'input must be an absolute path: {path}')
        seen, counts = {}, Counter()
        with path.open(encoding='utf-8', newline='') as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                if row['SAMPLE'] != entry['animal_id'] or row['TISSUE'] != entry['term']:
                    raise ValueError(f'Sample identity does not match the list: {path}')
                chrom = row['CHROM']
                if chrom not in lengths:
                    continue
                if row['PHASE_STATUS'] != 'PHASED_HET_EXACT_MATCH' or row['VCF_GT'] not in {'0|1', '1|0'}:
                    continue
                if int(row['TOTAL_COUNT']) <= 0:
                    continue
                if row['cnv_result_status'] != 'OK' or row['cnv_evaluable'] != 'TRUE':
                    continue
                pos = int(row['POS'])
                if not 1 <= pos <= lengths[chrom]:
                    raise ValueError(f'site exceeds FAI: {chrom}:{pos}')
                action, state = row['cnv_filter_action'], row['cnv_step7_consensus_state']
                key = (chrom, pos, row['REF'], row['ALT'])
                signature = (action, state, row['TOTAL_COUNT'], row['VCF_GT'],
                             row['cnv_step7_consensus_support'], row['cnv_step7_consensus_required_support'])
                if key in seen:
                    if seen[key] != signature:
                        raise ValueError(f'sample within duplicate site information conflict: {path} {key}')
                    continue
                seen[key] = signature
                if action != 'exclude_cnv_overlap':
                    continue
                support = int(row['cnv_step7_consensus_support'])
                if state not in STATES or int(row['cnv_step7_consensus_required_support']) != 7 or support < 7:
                    raise ValueError(f'exclude site does not satisfy current three type and 7/9 rule: {path} {key}')
                counts[state] += 1
                rows.append(dict(sample=entry['animal_id'], tissue=entry['term'],
                                 chromosome=chrom, pos1=pos, ref=row['REF'], alt=row['ALT'],
                                 state=state, consensus_support=support))
        summaries.append(dict(sample=entry['animal_id'], tissue=entry['term'],
                              eligible_unique_sites=len(seen), input_tsv=str(path),
                              **{state: counts[state] for state in STATES}))
        logging.info('%s/%s: passes requirements unique site=%s; exclude amp=%s del=%s cnloh=%s',
                     entry['animal_id'], entry['term'], len(seen), counts['amp'], counts['del'], counts['cnloh'])
    return rows, summaries


def bin_positions(rows, lengths):

    positions = {(r['tissue'], r['state'], r['chromosome'], r['pos1']) for r in rows}
    counts = Counter((tissue, state, chrom, (pos - 1) // BIN_BP)
                     for tissue, state, chrom, pos in positions)
    return [dict(tissue=tissue, state=state, chromosome=chrom,
                 start0=start, end0=min(start + BIN_BP, length),
                 width_bp=min(BIN_BP, length - start),
                 unique_position_count=counts[tissue, state, chrom, start // BIN_BP])
            for tissue in TISSUES for state in STATES
            for chrom, length in lengths.items() for start in range(0, length, BIN_BP)]


R_RENDER = r"""
suppressPackageStartupMessages(library(circlize))
args <- commandArgs(trailingOnly=TRUE)
bins <- read.delim(args[1], check.names=FALSE)
chroms <- read.delim(args[2], check.names=FALSE)
prefix <- args[3]
width_mm <- 87
height_mm <- 75

legend_height_mm <- 12
circle_box_mm <- height_mm - legend_height_mm - 2
circle_left_mm <- (width_mm - circle_box_mm) / 2

font_paths <- systemfonts::match_fonts(c("Arial", "Arial"), weight=c("normal", "bold"))$path
expected_fonts <- c("/path/to/fonts/Arial/Arial.TTF",
                    "/path/to/fonts/Arial/Arialbd.TTF")
if (!all(file.exists(expected_fonts)) || !identical(normalizePath(font_paths), normalizePath(expected_fonts)))
    stop("Arial regular/bold font resolution mismatch; check the font installation.")
tissues <- c("cerebellum", "cortex_frontal", "cortex_occipital", "pituitary_gland")
short <- c("CB", "FC", "OC", "Pit")
states <- c("amp", "del", "cnloh")
cols <- c("#D6604D", "#0072B2", "#009E73")
tissue_cols <- c("#BCD2E2", "#D9C5DE", "#CCDBBF", "#E7D3AD")
root_max <- max(1, ceiling(sqrt(max(bins$unique_position_count))))

render <- function() {

    par(fig=c(circle_left_mm/width_mm, (circle_left_mm+circle_box_mm)/width_mm,
              legend_height_mm/height_mm, (height_mm-2)/height_mm),
        mar=c(0,0,0,0), family="Arial",
        ps=10, cex=1, xpd=NA, lend="butt")
    circos.clear()
    circos.par(start.degree=90, gap.degree=0, cell.padding=c(0,0,0,0),
               track.margin=c(0.001,0.001), canvas.xlim=c(-1.18,1.18),
               canvas.ylim=c(-1.18,1.18))
    circos.genomicInitialize(chroms, plotType=NULL, sector.width=rep(1,nrow(chroms)))
    circos.genomicTrack(chroms, ylim=c(0,1), track.height=0.025,
        bg.border="#9AA1A7", bg.lwd=0.35,
        panel.fun=function(region,value,...) {
            circos.genomicRect(region, value, ybottom=0, ytop=1,
                               col="#E5E8EA", border=NA)
        })
    for (ti in seq_along(tissues)) {
        circos.genomicTrack(chroms, ylim=c(0,1), track.height=0.008,
            bg.border=NA, bg.col=tissue_cols[ti], panel.fun=function(region,value,...) {})
        for (si in seq_along(states)) {
            d <- bins[bins$tissue==tissues[ti] & bins$state==states[si],
                      c("chromosome","start0","end0","unique_position_count")]
            d[[4]] <- sqrt(d[[4]])
            circos.genomicTrack(d, ylim=c(0,root_max), numeric.column=4,
                track.height=0.058, bg.border=NA, bg.col="white",
                panel.fun=function(region,value,...) {
                    circos.lines(CELL_META$xlim,c(0,0),col="#DFE2E4",lwd=0.4)
                    keep <- value[[1]] > 0
                    if (any(keep)) circos.genomicLines(region[keep,,drop=FALSE],
                        value[keep,,drop=FALSE], numeric.column=1, type="h",
                        baseline=0, col=cols[si], lwd=0.5)
                })
        }
    }

    for (i in seq_len(nrow(chroms))) {
        degrees <- 90 - (i - 0.5) * 360 / nrow(chroms)
        angle <- degrees * pi / 180
        rotation <- (degrees - 90 + 180) %% 360 - 180
        if (rotation < -90) rotation <- rotation + 180
        if (rotation > 90) rotation <- rotation - 180
        text(1.08*cos(angle),1.08*sin(angle),sub("chr","",chroms[[1]][i]),
             srt=rotation,adj=c(0.5,0.5),cex=0.65)
    }

    label <- function(x, y, text, size=7, face="plain")
        grid::grid.text(text, x=grid::unit(x,"mm"), y=grid::unit(y,"mm"),
                        just="left", gp=grid::gpar(fontfamily="Arial", fontsize=size, fontface=face))
    grid::grid.text("B", x=grid::unit(2,"mm"),
                    y=grid::unit(height_mm-2,"mm"), just=c("left","top"),
                    gp=grid::gpar(fontfamily="Arial",fontsize=10,fontface="bold"))
    label(3, 9, "Tissue", face="bold")
    for (i in seq_along(short)) {
        x <- c(23,39,55,71)[i]
        grid::grid.rect(x=grid::unit(x,"mm"), y=grid::unit(9,"mm"),
                        width=grid::unit(2,"mm"), height=grid::unit(2,"mm"),
                        gp=grid::gpar(fill=tissue_cols[i],col=NA))
        label(x+2, 9, short[i])
    }
    label(3, 3.5, "CNV", face="bold")
    for (i in seq_along(states)) {
        x <- c(15,46,69)[i]
        grid::grid.lines(x=grid::unit(c(x,x+3),"mm"), y=grid::unit(c(3.5,3.5),"mm"),
                         gp=grid::gpar(col=cols[i],lwd=0.8))
        label(x+4, 3.5, c("Amplification","Deletion","cnLOH")[i])
    }
    circos.clear()
}

cairo_pdf(paste0(prefix,".pdf"),width=width_mm/25.4,height=height_mm/25.4,
          family="Arial",pointsize=10,bg="white")
render(); dev.off()
svglite::svglite(paste0(prefix,".svg"),width=width_mm/25.4,height=height_mm/25.4,pointsize=10,bg="white")
render(); dev.off()
png(paste0(prefix,".png"),width=width_mm,height=height_mm,units="mm",res=600,pointsize=10,bg="white",type="cairo")
render(); dev.off()
tiff(paste0(prefix,".tiff"),width=width_mm,height=height_mm,units="mm",res=600,pointsize=10,bg="white",compression="lzw",type="cairo")
render(); dev.off()
writeLines(capture.output(sessionInfo()),paste0(prefix,".R_session.txt"))
"""


def draw(bins, lengths, run):

    rscript = Path('/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/Rscript')
    renderer = run / f'{STEM}.render.R'
    renderer.write_text(R_RENDER, encoding='utf-8')
    chrom_path = run / f'{STEM}.chromosomes.tsv'
    write_table(chrom_path, [dict(chromosome=c, start0=0, end0=n) for c,n in lengths.items()],
                ['chromosome', 'start0', 'end0'])
    logging.info('invoke circlize: twelve records categories with white backgrounds track; four records narrow tissue color band; 1 Mb count')
    subprocess.run([str(rscript), '--vanilla', str(renderer),
                    str(run / f'{STEM}.bins.tsv'), str(chrom_path), str(run / STEM)],
                   check=True, stdout=sys.stdout, stderr=sys.stderr)


def write_table(path, rows, fields):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main(run):
    logging.info('start; concurrent_tasks=1; threads_per_task=1; list=%s; FAI=%s', MANIFEST, FAI)
    lengths, entries = read_inputs()
    rows, summaries = aggregate(lengths, entries)
    write_table(run / f'{STEM}.source_data.tsv', rows,
                ['sample', 'tissue', 'chromosome', 'pos1', 'ref', 'alt', 'state', 'consensus_support'])
    write_table(run / f'{STEM}.sample_summary.tsv', summaries,
                ['sample', 'tissue', 'eligible_unique_sites', 'amp', 'del', 'cnloh', 'input_tsv'])
    bins = bin_positions(rows, lengths)
    write_table(run / f'{STEM}.bins.tsv', bins,
                ['tissue', 'state', 'chromosome', 'start0', 'end0', 'width_bp', 'unique_position_count'])

    maximum = max((row['unique_position_count'] for row in bins), default=0)
    root_limit = max(1, int(maximum ** 0.5))
    if root_limit ** 2 < maximum:
        root_limit += 1
    caption = (
        'B. Binned counts of CNV-filtered eligible ASE positions. '
        'Inputs are the complete 1.4 CNV-annotated variant tables from run '
        '20260821_154914, upstream of the 1.5 gene-level tables used in panel A. '
        'No selection by gene annotation or candidate-gene membership is applied. '
        'CanFam3.1 autosomes 1–38 occupy equal angular sectors in a closed circle. '
        'Within chromosomes, positions use the actual FAI length. Therefore Mb per '
        'display angle differs between chromosomes; arc widths and apparent areas '
        'must not be interpreted as directly comparable genomic spans or densities. '
        'Tissues outside inward: CB, FC, OC, Pit, each identified by a narrow coloured band. '
        'Within each group, white-background tracks outside inward represent amplification '
        '(red), deletion (blue), and cnLOH (green). '
        'Peaks count distinct CHROM/POS within each tissue/class pooled across dogs '
        'and alleles in nonoverlapping 1 Mb windows. Final chromosome windows may be shorter '
        'and are not length-normalized. All twelve tracks share one square-root height scale. '
        'Independent radial lines are drawn at actual bin midpoints using circlize '
        'circos.genomicLines(type="h", baseline=0); width is 0.5 R line-width units '
        'and does not encode genomic span. Midpoints are display positions, not necessarily '
        'observed variant coordinates. Zero-count bins have no line. Lines are not connected '
        'or smoothed. Zero is at each track inner edge; counts increase outward. '
        f'The common radial scale is 0 to {root_limit ** 2:,} positions per bin '
        '(square-root transformed); this scale is reported here rather than inside the panel. '
        'Peaks encode position counts, not copy number, CNV segments, recurrence '
        'across dogs, or the fraction of sites removed. A coordinate may appear in '
        'multiple tissue/class groups; summing groups is not a genome-wide unique count. '
        'Eligibility: PHASED_HET_EXACT_MATCH, VCF_GT 0|1 or 1|0, TOTAL_COUNT > 0, '
        'cnv_result_status=OK and cnv_evaluable=TRUE. Plotted records have '
        'cnv_filter_action=exclude_cnv_overlap and an upstream amp/del/cnloh '
        'consensus supported by at least seven of nine scales; CNV is not recalculated. '
        'Source and sample summary tables deduplicate CHROM/POS/REF/ALT within each '
        'sample, retaining sample and allele identity, including intergenic records. '
        'The binned drawing additionally deduplicates CHROM/POS within tissue/class. '
        'Source POS is 1-based; bin coordinates are 0-based half-open. '
        'CB/FC/OC pool five dogs each; Pit pools four (Crak absent). Raw counts '
        'are not normalized for sample number, coverage or the number of eligible sites. '
        'Empty background means no plotted excluded eligible positions, not proven '
        'normal copy number. Labels are upstream RNA/BAF-inferred classes, not '
        'independently validated absolute DNA copy numbers.\n\n'
        f'Input manifest (full_tsv column): {MANIFEST}\nFAI: {FAI}\n'
        'Colour reference (red amplification / blue deletion): '
        'https://datacommons.cancer.gov/media/253 ; green cnLOH follows the original figure in directory 1.\n'
    )
    (run / f'{STEM}.caption.md').write_text(caption, encoding='utf-8')
    logging.info('plot the complete across-dog summary; circlize four tissue each three category; export PDF/SVG/PNG/TIFF')
    draw(bins, lengths, run)

    layout = dict(width_mm=87, height_mm=75, raster_dpi=600, font_family='Arial',
                  layout_role="first_row_75mm_trial",
                  circle_box_mm=dict(left=13, top=2, width=61, height=61),
                  legend_region_height_mm=12,
                  font_sizes_pt=dict(panel=10, chromosome=6.5, legend=7),
                  inputs=dict(manifest=str(MANIFEST), fai=str(FAI),
                              variant_tables=[entry['full_tsv'] for entry in entries]),
                  statistics='unchanged; 1 Mb bins; shared square-root radial scale',
                  review=dict(font_embedding='NOT_CHECKED', collisions='NOT_CHECKED',
                              figure2_assembly='NOT_CHECKED', figure_and_caption='NOT_CHECKED'))
    (run / 'inputs_and_layout.json').write_text(json.dumps(layout, indent=2)+'\n', encoding='utf-8')
    (run / 'execution_completion.marker').write_text('execution_status=success\n', encoding='utf-8')
    logging.info('task completed / done; execution_status=success')
    logging.info('please separately review exclude site, plot layout and biological interpretation; execution completed does not imply quality passed.')


if __name__ == '__main__':
    run_dir = ROOT / datetime.now().strftime('%Y%m%d_%H%M%S')
    run_dir.mkdir(parents=True, exist_ok=False)
    with (run_dir / 'plot_figure2B.log').open('x', encoding='utf-8', buffering=1) as handle:
        with redirect_stdout(handle), redirect_stderr(handle):
            logging.basicConfig(stream=handle, level=logging.INFO,
                                format='[%(asctime)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
            try:
                main(run_dir)
            except Exception:
                logging.exception('execution_status=failed; successful completion marker was not written.')
                sys.exit(1)
