#!/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/Rscript


Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           NUMEXPR_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")

ROOT <- "/path/to/repository/scripts/13_figures/historical_figure_numbers/figure3/D_regional_finemapping"
FINEMAP <- "/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.10_run_dog_behavior_susie_rss_finemapping/20260802_205352/main/primary/score_energy/DOG_GWAS_PRIMARY__score_energy__chr15_001"
REGIONS <- "/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.7_define_dog_behavior_ld_loci/20260802_171634/all_locus_regions.tsv"
CANDIDATES <- "/path/to/project/analysis_20260817/results/5-dogs/8.new-gene-allelome-eqtl-gwas-finemap-integration/summary/results/20260910_201222/dog5_candidate_records.tsv"
GTF <- "/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/GCF_000002285.3_CanFam3.1_genomic_chr.plus_five_only_new_gene.gtf"
LOCUS <- "DOG_GWAS_PRIMARY__score_energy__chr15_001"
WINDOW <- c(40850000, 41360000)
GENE_IDS <- c("gene-LOC111090065", "DOG5_NEWGENE_01873", "gene-IGF1")
GENE_NAMES <- c("LOC111090065", "DOG5_NEWGENE_01873", "IGF1")
FONT <- "Arial"
WIDTH_MM <- 178
HEIGHT_MM <- 84
BLUE <- "#21618C"
ORANGE <- "#D55E00"

read_tsv <- function(path) {
    utils::read.delim(path, quote="", comment.char="", check.names=FALSE,
                      stringsAsFactors=FALSE)
}
stage <- function(message) {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), message, "\n")
    flush(log_connection)
}


read_selected_exons <- function() {
    connection <- file(GTF, "r")
    on.exit(close(connection))
    parts <- list()
    repeat {
        lines <- readLines(connection, n=50000L, warn=FALSE)
        if (!length(lines)) break
        lines <- lines[startsWith(lines, "chr15\t") & grepl("\texon\t", lines, fixed=TRUE)]
        if (!length(lines)) next
        ids <- rep("", length(lines))

        for (id in GENE_IDS)
            ids[grepl(paste0('(?:^|[;\t])[[:space:]]*gene_id "', id, '"'),
                      lines, perl=TRUE)] <- id
        keep <- ids != ""
        if (!any(keep)) next
        fields <- strsplit(lines[keep], "\t", fixed=TRUE)
        parts[[length(parts)+1L]] <- data.frame(
            gene_id=ids[keep],
            start=as.integer(vapply(fields, `[`, character(1), 4L)),
            end=as.integer(vapply(fields, `[`, character(1), 5L)),
            strand=vapply(fields, `[`, character(1), 7L),
            stringsAsFactors=FALSE)
    }
    if (!length(parts)) stop("GTF in not found three items specified gene of chr15 exon")
    exons <- do.call(rbind, parts)
    if (!all(GENE_IDS %in% exons$gene_id)) stop("GTF missing specified gene of exon")
    exons
}


union_exons <- function(exons) {
    rows <- lapply(seq_along(GENE_IDS), function(i) {
        gene <- exons[exons$gene_id == GENE_IDS[i], ]
        gene <- gene[order(gene$start, gene$end), ]
        intervals <- list()
        left <- gene$start[1]; right <- gene$end[1]
        if (nrow(gene) > 1L) for (j in 2:nrow(gene)) {
            if (gene$start[j] <= right + 1L) right <- max(right, gene$end[j]) else {
                intervals[[length(intervals)+1L]] <- c(left, right)
                left <- gene$start[j]; right <- gene$end[j]
            }
        }
        intervals[[length(intervals)+1L]] <- c(left, right)
        spans <- do.call(rbind, intervals)
        strands <- unique(gene$strand)
        data.frame(gene_id=GENE_IDS[i], gene_name=GENE_NAMES[i],
                   strand=if (length(strands)==1L) strands else "mixed",
                   start=spans[,1], end=spans[,2], model="gene_exon_union")
    })
    do.call(rbind, rows)
}

draw_figure <- function(variants, genes, lead_id, focus_id, cs_size) {
    xlim <- WINDOW / 1e6
    lead <- variants[variants$ld_variant_id == lead_id, ]
    focus <- variants[variants$ld_variant_id == focus_id, ]
    is_cs <- variants$in_credible_set == "yes"

    graphics::layout(matrix(1:4, ncol=1), heights=c(1.30, 1.35, 1.15, 0.32))


    graphics::par(family=FONT, cex=1, cex.axis=7/8, cex.lab=7.5/8,
                  mar=c(0.4, 4.0, 2.0, 1.0), mgp=c(2.5, 0.5, 0),
                  las=1, xaxs="i", yaxs="i", tcl=-0.2)
    y_assoc <- -log10(variants$p_lrt)
    graphics::plot(variants$pos/1e6, y_assoc, type="n", xlim=xlim,
                   ylim=c(0, max(y_assoc)*1.17), xaxt="n", bty="l", xlab="",
                   ylab="-log10(P_LRT)")
    graphics::abline(h=-log10(c(2.97e-6, 1e-5)), col="#B9B9B9",
                     lty=c(2, 3), lwd=0.6)
    graphics::points(variants$pos/1e6, y_assoc, pch=16, cex=0.53, col="#B2B2B2")
    graphics::points(variants$pos[is_cs]/1e6, y_assoc[is_cs],
                     pch=16, cex=0.66, col=BLUE)
    graphics::points(lead$pos/1e6, -log10(lead$p_lrt), pch=18, cex=1.03, col="#202020")
    graphics::points(focus$pos/1e6, -log10(focus$p_lrt),
                     pch=1, cex=1.17, lwd=1, col=ORANGE)

    tag_x <- graphics::grconvertX(5/72, from="inches", to="user")
    graphics::mtext("D", side=3, at=tag_x, adj=0, line=0.35,
                    font=2, cex=10/8, xpd=NA)
    graphics::mtext("Energy | chr15", side=3, adj=1, line=0.35, cex=1)

    graphics::par(mar=c(0.4, 4.0, 0.2, 1.0))
    graphics::plot(variants$pos/1e6, variants$pip, type="n", xlim=xlim,
                   ylim=c(0, 0.20), xaxt="n", yaxt="n", bty="l", xlab="",
                   ylab="PIP")
    graphics::axis(2, at=seq(0, 0.20, 0.05))
    graphics::points(variants$pos/1e6, variants$pip, pch=16, cex=0.53, col="#B2B2B2")
    graphics::points(variants$pos[is_cs]/1e6, variants$pip[is_cs],
                     pch=16, cex=0.66, col=BLUE)
    graphics::points(lead$pos/1e6, lead$pip, pch=18, cex=1.03, col="#202020")
    graphics::points(focus$pos/1e6, focus$pip, pch=1, cex=1.17, lwd=1, col=ORANGE)
    graphics::text(lead$pos/1e6, lead$pip+0.014,
                   labels=sprintf("%s  PIP %.3f", lead_id, lead$pip),
                   cex=7/8, adj=c(0.5, 0), col="#202020")

    label_x <- focus$pos/1e6 + 0.018
    label_y <- 0.096
    graphics::segments(focus$pos/1e6, focus$pip+0.006,
                       label_x-0.003, label_y-0.006, col=ORANGE, lwd=0.65)
    graphics::text(label_x, label_y,
                   labels=sprintf("%s  PIP %.3f", focus_id, focus$pip),
                   cex=7/8, adj=c(0, 0), col=ORANGE)


    graphics::par(mar=c(2.7, 4.0, 0.1, 1.0))
    graphics::plot.new()
    graphics::plot.window(xlim=xlim, ylim=c(0.35, 3.65), xaxs="i", yaxs="i")
    for (i in seq_along(GENE_IDS)) {
        gene <- genes[genes$gene_id == GENE_IDS[i], ]
        y <- 4-i
        color <- if (GENE_IDS[i] == "DOG5_NEWGENE_01873") ORANGE else "#626262"
        start <- max(min(gene$start), WINDOW[1])/1e6
        end <- min(max(gene$end), WINDOW[2])/1e6
        graphics::segments(start, y, end, y, col=color, lwd=0.7)
        visible <- gene[gene$start <= WINDOW[2] & gene$end >= WINDOW[1], ]
        if (nrow(visible)) graphics::rect(
            pmax(visible$start, WINDOW[1])/1e6, y-0.12,
            pmin(visible$end, WINDOW[2])/1e6, y+0.12,
            col=color, border=NA)
        graphics::text((start+end)/2, y+0.36, labels=GENE_NAMES[i],
                       cex=7/8, col=color)
    }
    ticks <- seq(40.9, 41.3, by=0.1)
    graphics::axis(1, at=ticks, labels=sprintf("%.1f", ticks), cex.axis=7/8)
    graphics::mtext("Chr15 position (Mb), CanFam3.1", side=1, line=1.8, cex=7.5/8)

    graphics::par(mar=c(0, 4.0, 0, 1.0))
    graphics::plot.new()
    graphics::legend("center",
        legend=c(sprintf("95%% CS (n=%d)", cs_size), "Other modeled",
                 "Lead / top PIP", "Matched SNP"),
        pch=c(16, 16, 18, 1), col=c(BLUE, "#B2B2B2", "#202020", ORANGE),
        cex=7/8, bty="n", horiz=TRUE, x.intersp=0.7)
}

main <- function(out) {

    font_paths <- systemfonts::match_fonts(c(FONT,FONT),weight=c("normal","bold"))$path
    expected_fonts <- c("/path/to/fonts/Arial/Arial.TTF",
                        "/path/to/fonts/Arial/Arialbd.TTF")
    if (!all(file.exists(expected_fonts)) ||
        !identical(normalizePath(font_paths),normalizePath(expected_fonts)))
        stop("Arial regular or bold font was not resolved correctly, please check font installation")
    stage("stage 1/3 read main branch fine-mapping with candidate source; concurrent_tasks=1 threads_per_task=1")
    variants <- read_tsv(file.path(FINEMAP, "variant_pip.tsv.gz"))
    cs <- read_tsv(file.path(FINEMAP, "credible_sets.tsv"))
    summary <- read_tsv(file.path(FINEMAP, "locus_summary.tsv"))
    regions <- read_tsv(REGIONS)
    candidates <- read_tsv(CANDIDATES)
    cs <- cs[cs$locus_id == LOCUS & cs$branch == "main" & cs$cs_id == "L1", ]
    summary <- summary[summary$locus_id == LOCUS & summary$branch == "main", ]
    region <- regions[regions$locus_id == LOCUS, ]
    candidate <- candidates[candidates$new_gene_id == "DOG5_NEWGENE_01873" &
                            candidates$gwas_locus_id == LOCUS &
                            candidates$main_result == "yes", ]
    if (nrow(cs)!=1L || nrow(summary)!=1L || nrow(region)!=1L || nrow(candidate)!=1L)
        stop("Energy region, main branch L1 or 01873 candidate does not have a unique match")
    focus_id <- strsplit(candidate$gwas_overlap_snp[1], " | ", fixed=TRUE)[[1]][1]
    lead_id <- summary$top_variant_id[1]
    if (any(!is.finite(variants$p_lrt) | variants$p_lrt<=0) ||
        any(!is.finite(variants$pip) | variants$pip<0 | variants$pip>0.20))
        stop("input P value or PIP not suitable for the current axes")
    shown <- variants[variants$pos >= WINDOW[1] & variants$pos <= WINDOW[2], ]
    if (nrow(shown[shown$ld_variant_id==lead_id, ])!=1L ||
        nrow(shown[shown$ld_variant_id==focus_id, ])!=1L ||
        sum(shown$in_credible_set=="yes") != cs$cs_size[1] ||
        WINDOW[1] < region$analysis_start[1] || WINDOW[2] > region$analysis_end[1])
        stop("the focused window does not fully cover the credible set, specified variant, or exceeds original analysis region")

    stage("stage 2/3 read CanFam3.1 merge GTF and compute three items gene of exon union")
    genes <- union_exons(read_selected_exons())
    utils::write.table(shown, file.path(out, "Figure3D.plotted_variants.tsv"),
                       sep="\t", quote=FALSE, row.names=FALSE, na="NA")
    utils::write.table(genes, file.path(out, "Figure3D.gene_exon_union.tsv"),
                       sep="\t", quote=FALSE, row.names=FALSE, na="NA")

    stage("stage 3/3 plot 178×84 mm focus region figure with figure legend; no model refitting")
    pdf <- file.path(out, "Figure3D_energy_regional_finemapping.pdf")
    png <- file.path(out, "Figure3D_energy_regional_finemapping.png")
    grDevices::cairo_pdf(pdf, width=WIDTH_MM/25.4, height=HEIGHT_MM/25.4,
                         pointsize=8, family=FONT)
    draw_figure(shown, genes, lead_id, focus_id, cs$cs_size[1])
    grDevices::dev.off()
    grDevices::png(png, width=WIDTH_MM/25.4, height=HEIGHT_MM/25.4, units="in",
                   res=600, pointsize=8, type="cairo")
    draw_figure(shown, genes, lead_id, focus_id, cs$cs_size[1])
    grDevices::dev.off()


    svglite::svglite(file.path(out,"Figure3D_energy_regional_finemapping.svg"),
                    width=WIDTH_MM/25.4,height=HEIGHT_MM/25.4,pointsize=8,bg="white")
    draw_figure(shown, genes, lead_id, focus_id, cs$cs_size[1])
    grDevices::dev.off()
    ragg::agg_tiff(file.path(out,"Figure3D_energy_regional_finemapping.tiff"),
                   width=WIDTH_MM,height=HEIGHT_MM,units="mm",res=600,
                   pointsize=8,background="white",compression="lzw")
    draw_figure(shown, genes, lead_id, focus_id, cs$cs_size[1])
    grDevices::dev.off()
    utils::write.table(data.frame(
        role=c("variant_pip","credible_sets","locus_summary","regions","candidate_records","merged_gtf"),
        path=c(file.path(FINEMAP,c("variant_pip.tsv.gz","credible_sets.tsv","locus_summary.tsv")),
               REGIONS,CANDIDATES,GTF)),file.path(out,"Figure3D.input_manifest.tsv"),
        sep="\t",quote=FALSE,row.names=FALSE)
    writeLines(c("width_mm=178; height_mm=84; font=Arial",
        "font_pt: panel=10; title=8; axis_title=7.5; ticks/variant/gene/legend=7",
        "association_axis_label=plain_text; no_small_math_subscripts",
        "panel_tag_left_margin=5 pt (matches A canvas margin); tag_font=10 pt bold",
        "track_height_ratios=1.30,1.35,1.15,0.32; genomic_window=40850000-41360000",
        "PDF/SVG vector; PNG/TIFF 600 dpi; TIFF RGB LZW",
        "final_size_collision_and_assembly_QA=NOT_CHECKED",
        paste("font_paths",paste(font_paths,collapse=";"))),
        file.path(out,"Figure3D.layout_parameters.txt"))
    writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))

    candidate_gene <- genes[genes$gene_id=="DOG5_NEWGENE_01873", ]
    distance_bp <- shown$pos[shown$ld_variant_id==focus_id] - max(candidate_gene$end)
    gene_strands <- vapply(GENE_IDS, function(id) unique(genes$strand[genes$gene_id==id]),
                           character(1))
    caption <- c(
        "Figure 3D. Energy chr15 main-branch GWAS regional fine-mapping (CanFam3.1).",
        sprintf("This expands the Energy L1 credible set marked in panel C. Focused window: chr15:%d-%d; %d of %d modeled variants and all %d L1 95%% credible-set members are shown.",
                WINDOW[1], WINDOW[2], nrow(shown), nrow(variants), cs$cs_size[1]),
        "Top: -log10(P_LRT) of modeled variants; dashed and dotted lines denote P_LRT=2.97e-6 and 1e-5, respectively. Middle: GWAS PIP from the same main-branch model. Blue points are exported credible-set members; gray points are other modeled variants. The black diamond marks the GWAS lead and highest-PIP variant.",
        sprintf("The orange ring marks %s, selected from downstream integration, with GWAS P_LRT=%.6g and PIP=%.3f; it is not the highest-PIP variant.",
                focus_id, shown$p_lrt[shown$ld_variant_id==focus_id], shown$pip[shown$ld_variant_id==focus_id]),
        sprintf("The matched variant is %d bp beyond the DOG5_NEWGENE_01873 interval. Orange denotes a candidate locus and downstream matched SNP, not demonstrated regulation or causality.", distance_bp),
        "Bottom: selected LOC111090065, DOG5_NEWGENE_01873 and IGF1 exon-union tracks; other genes are omitted from this annotation track, not absent from the region. Horizontal positions and exon widths remain to scale.",
        sprintf("GTF strand fields: LOC111090065 %s, DOG5_NEWGENE_01873 %s, IGF1 %s. Strand arrows are not drawn.",
                gene_strands[1], gene_strands[2], gene_strands[3]),
        "Layout: 178 x 84 mm; Arial; vector PDF/SVG and 600 dpi PNG/LZW TIFF. The top axis uses plain text -log10(P_LRT) to retain readable lettering.",
        "The 95% figure refers to target credible-set coverage, not each variant's PIP. This focused view omits modeled variants outside its window; it does not represent all variants tested in the original GWAS."
    )
    writeLines(caption, file.path(out, "Figure3D.caption.txt"))
}


run_id <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
log_dir <- file.path(ROOT, "log")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive=TRUE)
log_connection <- file(file.path(log_dir, paste0("plot_figure3D.", run_id, ".log")), "wt")
sink(log_connection)
sink(log_connection, type="message")
status <- tryCatch({
    out <- file.path(ROOT, run_id)
    if (!dir.create(out, recursive=FALSE)) stop("cannot create independent result directory")
    stage(paste("start; output directory", out))
    main(out)
    writeLines("execution_status=success", file.path(out, "execution_completion.marker"))
    stage("task completed / done; execution_status=success")
    stage("please check separately output completeness, font layout, coordinates and biological interpretation; execution success does not imply QC passed.")
    0L
}, error=function(e) {
    stage(paste("execution_status=failed", conditionMessage(e)))
    1L
})
sink(type="message")
sink()
close(log_connection)
quit(save="no", status=status)
