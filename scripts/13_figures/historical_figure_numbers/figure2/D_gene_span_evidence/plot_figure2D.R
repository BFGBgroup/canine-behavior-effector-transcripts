#!/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/Rscript


Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           NUMEXPR_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
ROOT <- "/path/to/repository/scripts/13_figures/historical_figure_numbers/figure2/D_gene_span_evidence"
MANIFEST <- "/path/to/repository/scripts/08_expression_ase_integration/preparation/1.5ase-cnv-clean-gene-id-haplotype-table/log/1.5.ase_cnv_clean_gene_id_haplotype_table.20260821_161805.success.tsv"
REFERENCE <- "https://github.com/charitylaw/Intron-reads/blob/master/Analyses/analysis_of_intron_exploration.R"
KNOWN_ANNOTATION <- "/path/to/project/analysis_20260709/5-dogs/3.ASE-process/3.7_build_canfam31_gene_span_reference/CanFam3.1.autosomes.gene_span.annotation.tsv"
NEW_GTF <- "/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.new_only.gtf"
DOGS <- c("Ace_396226", "CodRou_346627", "Crak_358336", "Dozer_386519", "SanMcg_349891")
TISSUES <- c("cerebellum", "cortex_frontal", "cortex_occipital", "pituitary_gland")
CLASSES <- c("exon_only", "exon_and_intron", "intron_only")

COLORS <- c("#0072B2", "#009E73", "#D55E00")
RECORD_TYPES <- c("known_gene", "manual_new_gene")
WIDTH_MM <- 87
HEIGHT_MM <- 75
DPI <- 600
FONT <- "Arial"

configure_font <- function() {
    paths <- systemfonts::match_fonts(c(FONT,FONT), weight=c("normal","bold"))$path
    expected <- c("/path/to/fonts/Arial/Arial.TTF",
                  "/path/to/fonts/Arial/Arialbd.TTF")
    if (!all(file.exists(expected)) || !identical(normalizePath(paths),normalizePath(expected)))
        stop("Arial regular or bold font was not matched correctly, do not silently substitute font")
    paths
}

read_tsv <- function(path) {
    read.delim(path, colClasses="character", na.strings=character(), quote="",
               comment.char="", check.names=FALSE, stringsAsFactors=FALSE)
}
write_tsv <- function(x, path) {
    write.table(x, path, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
}
log_stage <- function(x) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), x, "\n")
    flush.console()
}

load_structure <- function() {
    k <- read_tsv(KNOWN_ANNOTATION)
    known <- data.frame(gene_id=k$GeneID, chromosome=k$CHROM,
        exon_union_blocks=as.numeric(k$ExonUnionBlockN),
        exon_union_bp=as.numeric(k$ExonUnionBp), span_bp=as.numeric(k$SpanBp))

    gtf <- read.delim(NEW_GTF, header=FALSE, sep="\t", quote="", comment.char="#",
                      colClasses="character", na.strings=character(), check.names=FALSE)
    gtf <- gtf[gtf$V3 %in% c("transcript", "exon"), ]

    gtf$gene_id <- vapply(strsplit(gtf$V9, ";", fixed=TRUE), function(attributes) {
        fields <- trimws(attributes)
        field <- fields[grepl('^gene_id[[:space:]]+"[^"]+"$', fields)]
        if (length(field) != 1L) stop("candidate GTF missing unique gene_id attribute")
        sub('^gene_id[[:space:]]+"([^"]+)"$', '\\1', field)
    }, character(1))
    manual <- lapply(split(gtf, gtf$gene_id), function(x) {
        if (length(unique(x$V1)) != 1L)
            stop(paste("candidate gene across chromosome:", x$gene_id[1], paste(unique(x$V1), collapse=",")))
        ex <- x[x$V3 == "exon", ]
        if (!nrow(ex)) stop("candidate gene missing exon")
        starts <- as.numeric(ex$V4); ends <- as.numeric(ex$V5)
        ord <- order(starts, ends); starts <- starts[ord]; ends <- ends[ord]
        left <- starts[1]; right <- ends[1]; blocks <- 1L; bp <- 0
        if (length(starts) > 1L) for (j in 2:length(starts)) {
            if (starts[j] <= right + 1) right <- max(right, ends[j]) else {
                bp <- bp + right - left + 1
                left <- starts[j]; right <- ends[j]; blocks <- blocks + 1L
            }
        }
        data.frame(gene_id=x$gene_id[1], chromosome=x$V1[1],
                   exon_union_blocks=blocks, exon_union_bp=bp + right - left + 1,
                   span_bp=max(as.numeric(x$V5)) - min(as.numeric(x$V4)) + 1)
    })
    structure <- rbind(known, do.call(rbind, manual))
    structure$intronic_region_bp <- structure$span_bp - structure$exon_union_bp
    if (anyDuplicated(structure$gene_id) || any(structure$intronic_region_bp < 0))
        stop("structure annotation contains duplicates ID or invalid exon union length")
    structure
}

sample_evidence <- function(entry, structure) {
    g <- read_tsv(entry$output_tsv)
    qc <- read_tsv(entry$qc_tsv)
    clean_path <- qc$value[qc$metric == "input_tsv"]
    if (length(clean_path) != 1L || !startsWith(clean_path, "/"))
        stop("1.5 QC not provided unique absolute input_tsv path")
    g <- g[g$record_type %in% c("known_gene", "manual_new_gene") &
           as.numeric(g$variant_count) >= 1 & as.numeric(g$total_count_sum) > 0, ]
    if (any(g$sample != entry$animal_id | g$tissue != entry$tissue) ||
        anyDuplicated(g$gene_id)) stop("1.5 sample identity or unique gene_id does not satisfy A of count rule")
    v <- read_tsv(clean_path)
    if (any(v$SAMPLE != entry$animal_id | v$TISSUE != entry$tissue))
        stop("1.4 site table sample identity mismatch")

    v <- v[v$GeneID %in% g$gene_id & v$cnv_clean_pass == "TRUE" &
           v$PHASE_STATUS == "PHASED_HET_EXACT_MATCH" &
           v$VCF_GT %in% c("0|1", "1|0") & as.numeric(v$TOTAL_COUNT) > 0, ]
    if (any(v$PHASE_SET != v$CHROM) ||
        any(!v$GeneFeature %in% c("exon_union", "intron")))
        stop("target gene recorded phase or GeneFeature cannot interpretation")

    v <- v[!duplicated(v[c("GeneID", "CHROM", "POS", "REF", "ALT")]), ]
    evidence <- cbind(n_exon_sites=as.integer(v$GeneFeature == "exon_union"),
                      n_intron_sites=as.integer(v$GeneFeature == "intron"),
                      n_multi_gene_sites=as.integer(as.numeric(v$GeneHitN) > 1),
                      n_flagged_sites=as.integer(v$cnv_filter_action != "keep"))
    if (!nrow(v)) stop("A of gene set no can classification of positive count site")
    counts <- rowsum(evidence, group=v$GeneID, reorder=FALSE)
    ix <- match(g$gene_id, rownames(counts))
    if (anyNA(ix)) stop("subset A gene cannot in same source clean in the table find support, cannot silently reduce denominator")
    cts <- counts[ix, , drop=FALSE]
    category <- ifelse(cts[, 1] > 0 & cts[, 2] > 0, "exon_and_intron",
                       ifelse(cts[, 1] > 0, "exon_only", "intron_only"))
    si <- match(g$gene_id, structure$gene_id)
    if (anyNA(si)) stop("A gene set missing structure annotation")
    annotation <- structure[si, ]
    if (any(annotation$chromosome != g$chromosome)) stop("structure annotation chromosome mismatch")
    no_intron <- annotation$intronic_region_bp == 0
    if (any(no_intron & cts[, 2] > 0)) stop("without intronic regions gene encountered intron support")


    data.frame(sample=entry$animal_id, tissue=entry$tissue, gene_id=g$gene_id,
               record_type=g$record_type, chromosome=g$chromosome,
               category=category,
               exon_union_blocks=annotation$exon_union_blocks,
               exon_union_bp=annotation$exon_union_bp, span_bp=annotation$span_bp,
               intronic_region_bp=annotation$intronic_region_bp,
               cts, upstream_variant_count=g$variant_count,
               upstream_total_support=g$total_count_sum,
               input_gene_tsv=entry$output_tsv, input_clean_tsv=clean_path,
               row.names=NULL, check.names=FALSE)
}

draw_panel <- function(summary) {

    old <- par(no.readonly=TRUE)
    on.exit(par(old))
    spacing <- rep(0.22,19)
    spacing[c(1,6,11,16)] <- 1.1
    groups <- list(1:5,6:10,11:15,16:19)
    text_mm <- function(label,x,y,size=7,face="plain",rot=0,just="left") {
        grid::grid.text(label,x=grid::unit(x,"mm"),y=grid::unit(y,"mm"),
            just=just,rot=rot,gp=grid::gpar(fontfamily=FONT,fontsize=size,fontface=face))
    }
    for (i in seq_along(RECORD_TYPES)) {
        d <- summary[summary$record_type == RECORD_TYPES[i], ]
        if (nrow(d) != 19*length(CLASSES) ||
            !identical(as.character(d$category),rep(CLASSES,19)))
            stop("plot summary table missing sample or category order does not match")
        mat <- matrix(d$percentage,nrow=length(CLASSES))
        bottom <- if(i==1) 51 else 27
        top <- bottom+18
        par(fig=c(15/WIDTH_MM,85/WIDTH_MM,bottom/HEIGHT_MM,top/HEIGHT_MM),
            new=i>1,mar=c(0,0,0,0),mgp=c(1.6,0.35,0),
            tcl=-0.15,ps=10,cex=1,cex.axis=0.7,family=FONT,las=1,lwd=0.6,
            yaxs="i",xpd=NA)
        centers <- as.vector(barplot(mat,beside=FALSE,space=spacing,col=COLORS,
            border=NA,ylim=c(0,100),axes=FALSE,axisnames=FALSE))
        axis(2,at=seq(0,100,25),cex.axis=0.7); box(bty="l")

        xs <- grconvertX(centers,from="user",to="ndc")*WIDTH_MM
        text_mm(c("Reference genes","Candidate added loci")[i],15,top+2.5,8,"bold")
        if(i==2) {
            labels <- sub("_.*","",d$sample[seq(1,nrow(d),length(CLASSES))])
            for(j in seq_along(labels)) text_mm(labels[j],xs[j],25.5,6.5,rot=90,just="right")
            for(j in seq_along(groups)) text_mm(c("CB","FC","OC","Pit")[j],
                mean(xs[groups[[j]]]),13,7,"bold",just="centre")
        }
    }

    text_mm("Within-group genes (%)",3,(27+69)/2,7.5,rot=90,just="centre")

    grid::grid.text("D",x=grid::unit(2,"mm"),y=grid::unit(HEIGHT_MM-2,"mm"),
        just=c("left","top"),gp=grid::gpar(fontfamily=FONT,fontsize=10,fontface="bold"))
    legends <- c("Exon only","Exon + intron","Intron only")

    legend_x <- c(16.75,35.10,57.20)
    for(j in seq_along(legends)) {
        x <- legend_x[j]
        y <- 5
        grid::grid.rect(x=grid::unit(x,"mm"),y=grid::unit(y,"mm"),
            width=grid::unit(2,"mm"),height=grid::unit(2,"mm"),
            gp=grid::gpar(fill=COLORS[j],col=NA))
        text_mm(legends[j],x+2.5,y,7)
    }
}

main <- function(run) {
    log_stage("start; concurrent_tasks=1; threads_per_task=1")
    font_paths <- configure_font()
    entries <- read_tsv(MANIFEST)
    expected <- expand.grid(animal_id=DOGS, tissue=TISSUES, stringsAsFactors=FALSE)
    expected <- expected[!(expected$animal_id == DOGS[3] & expected$tissue == TISSUES[4]), ]
    if (nrow(entries) != 19 || anyDuplicated(entries[c("animal_id", "tissue")]) ||
        !setequal(paste(entries$animal_id, entries$tissue), paste(expected$animal_id, expected$tissue)))
        stop("manifest and 19 items actual sample combination does not match")
    entries <- entries[order(match(entries$tissue, TISSUES), match(entries$animal_id, DOGS)), ]
    log_stage("read reference gene structure and candidate GTF; compute exon union with intron region")
    structure <- load_structure()
    details <- summaries <- vector("list", nrow(entries))
    for (i in seq_len(nrow(entries))) {
        log_stage(paste("read", i, "/19", entries$animal_id[i], entries$tissue[i]))
        d <- sample_evidence(entries[i, ], structure)
        details[[i]] <- d
        summaries[[i]] <- do.call(rbind, lapply(RECORD_TYPES, function(type) {
            group <- d[d$record_type == type, ]
            if (!nrow(group)) stop(paste("sample missing evaluable genome, cannot plot zero proportion:",
                                        entries$animal_id[i], entries$tissue[i], type))
            n <- as.integer(table(factor(group$category, levels=CLASSES)))
            data.frame(sample=entries$animal_id[i], tissue=entries$tissue[i],
                record_type=type, category=CLASSES, n_genes=n,
                denominator=nrow(group), percentage=100*n/nrow(group))
        }))
    }
    detail <- do.call(rbind, details)
    summary <- do.call(rbind, summaries)
    write_tsv(detail, file.path(run, "Figure2D.gene_evidence.tsv"))
    write_tsv(summary, file.path(run, "Figure2D.sample_composition.tsv"))
    write_tsv(entries, file.path(run, "Figure2D.input_manifest.tsv"))
    writeLines(c("D. Location-based evidence composition of gene-span haplotype counts.",
        "The upper facet shows reference genes and the lower facet candidate added loci. Each facet contains one three-category stacked bar per dog-tissue sample. CB/FC/OC have five dogs and Pit four (Crak absent).",
        "The source gene set matches panel A: unique known_gene and manual_new_gene IDs with variant_count >=1 and positive total support in the 1.5 tables. Each bar is normalized separately within its sample and record_type; the two denominators sum to the panel A sample count. Denominators are provided in the source table.",
        "Reference/candidate identity is inherited from record_type in the 1.5 annotation contract, not inferred from evidence category or coding potential. Candidate added loci may include reference-related novel isoforms and are not all independent novel genes.",
        "Exon only, exon + intron, and intron only refer to the presence of retained positive-count phased heterozygous variant sites in gene-specific exon unions and the remaining gene span.",
        "Sites are deduplicated within each gene by CHROM/POS/REF/ALT. Multi-gene assignments and flags retained by the upstream CNV policy are retained here; they are recorded in the gene source table.",
        "These are proportions of genes by support location, not proportions of RNA reads, intron retention, or significant ASE. No additional depth or biallelic-SNP threshold is imposed.",
        "Exon-only genes are not subdivided by intronic-region length. Gene-span and exon-union structural measurements remain in the source table. Candidate annotations are not independently validated genes.",
        "Adapted classification concept from Lee et al., Covering all your bases (2020), GENE SIGNAL code. The original uses RNA counts >=3, flags genes without introns, and plots category points with all input genes as denominator. This independent implementation uses positive phased variant support, the panel A gene set, separate reference/candidate denominators and three stacked categories; it is not a reproduction of the original figure.",
        REFERENCE, paste("Manifest:", MANIFEST), paste("Known annotation:", KNOWN_ANNOTATION),
        paste("Candidate GTF:", NEW_GTF)), file.path(run, "Figure2D.caption.md"))
    log_stage("export 87 × 75 mm PDF/SVG and 600 dpi PNG/RGB TIFF; statistical denominator remains unchanged")
    cairo_pdf(file.path(run,"Figure2D_gene_span_evidence.pdf"),width=WIDTH_MM/25.4,
              height=HEIGHT_MM/25.4,pointsize=10,family=FONT,bg="white")
    draw_panel(summary); dev.off()
    svglite::svglite(file.path(run,"Figure2D_gene_span_evidence.svg"),width=WIDTH_MM/25.4,
                    height=HEIGHT_MM/25.4,pointsize=10,bg="white")
    draw_panel(summary); dev.off()
    png(file.path(run,"Figure2D_gene_span_evidence.png"),width=WIDTH_MM,height=HEIGHT_MM,
        units="mm",res=DPI,pointsize=10,bg="white",type="cairo")
    draw_panel(summary); dev.off()
    tiff(file.path(run,"Figure2D_gene_span_evidence.tiff"),width=WIDTH_MM,height=HEIGHT_MM,
         units="mm",res=DPI,pointsize=10,bg="white",type="cairo",compression="lzw")
    draw_panel(summary); dev.off()
    writeLines(c(paste("width_mm",WIDTH_MM),paste("height_mm",HEIGHT_MM),
        paste("raster_dpi",DPI),paste("font_regular",font_paths[1]),paste("font_bold",font_paths[2]),
        "font_pt: panel=10; titles=8; axis_title=7.5; ticks/legend/tissue=7; dog=6.5",
        "plot_area_mm: left=15; right=85; upper=51..69; lower=27..45",
        "shared_y_title_mm: x=3; y=48; font_pt=7.5",
        "legend: rows=1; columns=3; y_mm=5; font_pt=7; labels=Exon only|Exon + intron|Intron only",
        "font_embedding=NOT_CHECKED; collisions=NOT_CHECKED; assembled_figure2=NOT_CHECKED",
        "figure_and_caption_fit=NOT_CHECKED"),file.path(run,"Figure2D.layout_parameters.txt"))
    writeLines(capture.output(sessionInfo()), file.path(run, "R_session.txt"))
    writeLines("execution_status=success", file.path(run, "execution_completion.marker"))
    log_stage("task completed / done; execution_status=success")
    log_stage("please check separately denominator, classification, low-support and multi-gene-assignment records and the final layout; execution completed does not imply quality passed.")
}

run <- file.path(ROOT, format(Sys.time(), "%Y%m%d_%H%M%S"))
if (dir.exists(run) || !dir.create(run, recursive=TRUE)) stop("cannot create independent run directory")
log_con <- file(file.path(run, "plot_figure2D.log"), open="wt")
sink(log_con); sink(log_con, type="message")
status <- tryCatch({main(run); 0L}, error=function(e) {
    message("execution_status=failed: ", conditionMessage(e)); 1L
})
sink(type="message"); sink(); close(log_con)
quit(save="no", status=status)
