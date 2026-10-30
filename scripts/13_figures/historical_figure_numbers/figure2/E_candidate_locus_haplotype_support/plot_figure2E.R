#!/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/Rscript


Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           NUMEXPR_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
ROOT <- "/path/to/repository/scripts/13_figures/historical_figure_numbers/figure2/E_candidate_locus_haplotype_support"
PREP <- "/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation"
SITE_ROOT <- paste0(PREP,"/1.4ase-cnv-qc-filter/20260821_154914")
GENE_ROOT <- paste0(PREP,"/1.5ase-cnv-clean-gene-id-haplotype-table/20260821_161805")
GTF <- "/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.new_only.gtf"
PHASE_QC <- "/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.10_reconcile_dog11_whatsHap_shapeit5/block_qc/chr38.whatsHap_shapeit5.block_qc.tsv.gz"
GENE <- "DOG5_NEWGENE_03504"
DOG <- "Dozer_386519"
TISSUES <- c("cerebellum","cortex_frontal","cortex_occipital","pituitary_gland")
LABELS <- c("CB","FC","OC","Pit")
WIDTH_MM <- 87

HEIGHT_MM <- 60
FONT <- "Arial"
DPI <- 600

POINT_MIN_MM <- 0.6
POINT_MAX_MM <- 2.0
support_size_mm <- function(support, maximum)
    pmax(POINT_MIN_MM, POINT_MAX_MM * sqrt(support / maximum))
read_tsv <- function(p) read.delim(p, colClasses="character", quote="", comment.char="",
                                  na.strings=character(), check.names=FALSE)
write_tsv <- function(x,p) write.table(x,p,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
log_stage <- function(s) {cat(format(Sys.time(),"[%Y-%m-%d %H:%M:%S ]"),s,"\n"); flush.console()}


read_gene <- function(path, field) {
    con <- file(path,"rt"); on.exit(close(con))
    header <- readLines(con,n=1L)
    selected <- character()
    repeat {
        z <- readLines(con,n=50000L)
        if (!length(z)) break
        selected <- c(selected,z[grepl(GENE,z,fixed=TRUE)])
    }
    d <- read.delim(text=paste(c(header,selected),collapse="\n"),colClasses="character",
                    quote="",comment.char="",na.strings=character(),check.names=FALSE)
    d[d[[field]]==GENE,,drop=FALSE]
}
summary_counts <- function(d) data.frame(n_sites=nrow(d), hap1=sum(d$HAP1_COUNT),
    hap2=sum(d$HAP2_COUNT), total=sum(d$support),
    hap1_fraction=if(nrow(d)) sum(d$HAP1_COUNT)/sum(d$support) else NA_real_)

configure_font <- function() {
    paths <- systemfonts::match_fonts(c(FONT,FONT),weight=c("normal","bold"))$path
    expected <- c("/path/to/fonts/Arial/Arial.TTF",
                  "/path/to/fonts/Arial/Arialbd.TTF")
    if(!all(file.exists(expected)) || !identical(normalizePath(paths),normalizePath(expected)))
        stop("Arial regular/bold parse does not match, please check font installation")
    paths
}

draw_figure <- function(sites, exons, summary, out) {
    grid::grid.newpage()

    txt <- function(s,x,y,size=7,just="centre",face="plain",rot=0)
        grid::grid.text(s,grid::unit(x,"mm"),grid::unit(y,"mm"),just=just,rot=rot,
            gp=grid::gpar(fontfamily=FONT,fontsize=size,fontface=face,col="#222222"))
    line <- function(x0,y0,x1,y1,col="#B8BDC3",lty=1)
        grid::grid.segments(grid::unit(x0,"mm"),grid::unit(y0,"mm"),
            grid::unit(x1,"mm"),grid::unit(y1,"mm"),gp=grid::gpar(col=col,lwd=0.55,lty=lty))
    grid::grid.rect(gp=grid::gpar(fill="white",col=NA))
    txt("E",2,HEIGHT_MM-2,10,c("left","top"),"bold")
    txt(paste(GENE,"- Dozer"),12,57,7.5,"left")
    lo <- min(exons$start)-120; hi <- max(exons$end)+120
    left <- 17; right <- 59
    xpos <- function(pos) left+(pos-lo)/(hi-lo)*(right-left)
    txt("CanFam3.1, chr38 (Mb)",38,52.5,7)
    ticks <- pretty(c(lo,hi),n=4); ticks <- ticks[ticks>=lo & ticks<=hi]
    line(left,47.5,right,47.5)
    for(p in ticks) {
        line(xpos(p),47.5,xpos(p),48.2)
        txt(sprintf("%.3f",p/1e6),xpos(p),49.5,6.5)
    }
    txt("Candidate",15,44.5,6.5,"right")
    line(xpos(min(exons$start)),44.5,xpos(max(exons$end)),44.5,"#4879A6")
    for(i in seq_len(nrow(exons))) {
        grid::grid.rect(x=grid::unit(xpos((exons$start[i]+exons$end[i])/2),"mm"),
            y=grid::unit(44.5,"mm"),
            width=grid::unit((exons$end[i]-exons$start[i]+1)/(hi-lo)*(right-left),"mm"),
            height=grid::unit(0.8,"mm"),gp=grid::gpar(fill="#4879A6",col=NA))
    }

    ax <- xpos((exons$end[1]+exons$start[2])/2)
    line(ax-0.7,44.1,ax,44.5,"#4879A6")
    line(ax-0.7,44.9,ax,44.5,"#4879A6")

    txt("HAP1 (%)",74.25,43,7)
    txt("All",68.5,39.5,6.5); txt("Common",80,39.5,6.5)
    bottoms <- c(31,23,15,7); h <- 6
    txt("HAP1 fraction",2.5,22,7.5,rot=90)
    geometry <- list()
    max_support <- max(sites$support)
    for(j in seq_along(TISSUES)) {
        d <- sites[sites$TISSUE==TISSUES[j],]; b <- bottoms[j]

        txt(LABELS[j],8,b+h/2,7,face="bold")
        for(v in c(0,0.5,1)) {
            line(left,b+v*h,right,b+v*h,if(v==0.5) "#AEB5BC" else "#E4E7EB",if(v==0.5) 2 else 1)
            txt(format(v),15,b+v*h,6.5,"right")
        }

        d <- d[order(d$support,decreasing=TRUE),]
        diam <- support_size_mm(d$support,max_support)
        for(k in seq_len(nrow(d))) grid::grid.points(x=grid::unit(xpos(d$POS[k]),"mm"),
            y=grid::unit(b+d$hap1_fraction[k]*h,"mm"),
            pch=if(d$GeneFeature[k]=="exon_union") 21 else 24,
            size=grid::unit(diam[k],"mm"),
            gp=grid::gpar(col="#28679A",fill="#6699BE",alpha=0.8,lwd=0.35))
        a <- summary[summary$tissue==TISSUES[j] & summary$set=="all",]
        c <- summary[summary$tissue==TISSUES[j] & summary$set=="common",]
        txt(sprintf("%.2f",100*a$hap1_fraction),68.5,b+h/2,7)
        txt(sprintf("%.2f",100*c$hap1_fraction),80,b+h/2,7)
        geometry[[j]] <- data.frame(track=LABELS[j],
            left_mm=grid::convertX(grid::unit(left,"mm"),"mm",TRUE),
            right_mm=grid::convertX(grid::unit(right,"mm"),"mm",TRUE),
            bottom_mm=grid::convertY(grid::unit(b,"mm"),"mm",TRUE),
            top_mm=grid::convertY(grid::unit(b+h,"mm"),"mm",TRUE))
    }

    grid::grid.points(x=grid::unit(c(19,37),"mm"),y=grid::unit(c(3,3),"mm"),
        pch=c(21,24),size=grid::unit(1.6,"mm"),gp=grid::gpar(col="#28679A",fill="#6699BE"))
    txt("Exonic",21,3,7,"left"); txt("Intronic",39,3,7,"left")
    txt("Area: support",57,3,6.5,"left")
    write_tsv(do.call(rbind,geometry),file.path(out,paste0("Figure2E.track_geometry.",names(grDevices::dev.cur()),".tsv")))
}

main <- function(out) {
    font_paths <- configure_font()
    log_stage("stage 1/4 read target site with 1.5 summary; concurrent_tasks=1 threads_per_task=1")
    all <- list(); refs <- list(); manifest <- list(); excluded <- list()
    for (t in TISSUES) {
        prefix <- paste(DOG,t,sep="_")
        p <- file.path(SITE_ROOT,DOG,t,paste0(prefix,".ase_haplotype.cnv_clean.tsv"))
        g <- file.path(GENE_ROOT,DOG,t,paste0(prefix,".cnv_clean.gene_id_haplotype.tsv"))
        log_stage(p)
        d <- read_gene(p,"GeneID"); refs[[t]] <- read_gene(g,"gene_id")
        d$included <- d$cnv_clean_pass=="TRUE" & d$PHASE_STATUS=="PHASED_HET_EXACT_MATCH" &
            d$VCF_GT %in% c("0|1","1|0")
        excluded[[t]] <- d[!d$included,,drop=FALSE]
        d <- d[d$included,,drop=FALSE]
        if(!nrow(d) || any(d$PHASE_SET!=d$CHROM)) stop("phase fields or target-site inputs do not satisfy the required contract: ",t)
        for(n in c("POS","HAP1_COUNT","HAP2_COUNT","TOTAL_COUNT")) d[[n]] <- as.numeric(d[[n]])
        d$support <- d$HAP1_COUNT+d$HAP2_COUNT
        if(any(!is.finite(d$support) | d$support<=0)) stop("Cannot calculate site proportions: ",t)
        d$key <- paste(d$CHROM,d$POS,d$REF,d$ALT,sep=":")

        d <- unique(d)
        if(anyDuplicated(d$key)) stop("target gene in the same variant exists conflict record: ",t)
        d$hap1_fraction <- d$HAP1_COUNT/d$support
        all[[t]] <- d
        manifest[[t]] <- data.frame(tissue=t,site_input=p,gene_input=g)
    }
    sites <- do.call(rbind,all)
    common <- Reduce(intersect,lapply(all,function(d) d$key[d$GeneFeature=="exon_union"]))
    if(!length(common)) stop("four tissue no shared exon site")
    sites$common_exonic <- sites$key %in% common
    log_stage("stage 2/4 append original PS harmonization basis; do not change existing HAP orientation")
    q <- read_tsv(PHASE_QC)
    q <- q[q$sample_id==DOG & q$chromosome=="chr38" & q$ps %in% sites$PS,,drop=FALSE]
    if(anyDuplicated(q$ps) || any(!sites$PS %in% q$ps) || any(!q$orientation %in% c("keep","flip")))
        stop("current PS missing interpretable of orientation harmonize record")
    sites$block_orientation <- q$orientation[match(sites$PS,q$ps)]
    sites$phase_note <- "Counts already oriented upstream; do not flip again"
    rows <- list(); compare <- list()
    for(t in TISSUES) {
        d <- sites[sites$TISSUE==t,]; c <- d[d$common_exonic,]
        a <- summary_counts(d)
        rows[[paste0(t,"_all")]] <- cbind(tissue=t,set="all",a)
        rows[[paste0(t,"_common")]] <- cbind(tissue=t,set="common",summary_counts(c))
        r <- refs[[t]]
        if(nrow(r)!=1L) stop("1.5 target summary is not unique row: ",t)
        compare[[t]] <- data.frame(tissue=t,observed_sites=a$n_sites,upstream_sites=as.numeric(r$variant_count),
            observed_h1=a$hap1,upstream_h1=as.numeric(r$hap1_count_sum),
            observed_h2=a$hap2,upstream_h2=as.numeric(r$hap2_count_sum),
            matches_1_5=a$n_sites==as.numeric(r$variant_count) & a$hap1==as.numeric(r$hap1_count_sum) &
                a$hap2==as.numeric(r$hap2_count_sum))
    }
    summary <- do.call(rbind,rows)
    common_sites <- sites[sites$common_exonic,]
    pit <- common_sites[common_sites$TISSUE=="pituitary_gland",]
    delta <- common_sites[common_sites$TISSUE!="pituitary_gland",]
    delta$pit_hap1_fraction <- pit$hap1_fraction[match(delta$key,pit$key)]
    delta$delta_vs_pit <- delta$hap1_fraction-delta$pit_hap1_fraction
    log_stage("stage 3/4 save all site, shared site, phasing basis, count comparison and per site difference")
    write_tsv(sites,file.path(out,"Figure2E.all_sites.tsv"))
    write_tsv(common_sites,file.path(out,"Figure2E.common_sites.tsv"))
    write_tsv(do.call(rbind,excluded),file.path(out,"Figure2E.excluded_phase_sites.tsv"))
    write_tsv(summary,file.path(out,"Figure2E.summary.tsv"))
    write_tsv(q,file.path(out,"Figure2E.phase_evidence.tsv"))
    write_tsv(do.call(rbind,compare),file.path(out,"Figure2E.upstream_count_comparison.tsv"))
    write_tsv(delta,file.path(out,"Figure2E.common_site_differences.tsv"))
    write_tsv(do.call(rbind,manifest),file.path(out,"Figure2E.input_manifest.tsv"))
    gtf <- readLines(GTF)
    gtf <- gtf[grepl(paste0('gene_id "',GENE,'";'),gtf,fixed=TRUE)]
    models <- read.delim(text=paste(gtf,collapse="\n"),header=FALSE,sep="\t",quote="",comment.char="")
    exons <- data.frame(start=as.numeric(models[models$V3=="exon","V4"]),
                        end=as.numeric(models[models$V3=="exon","V5"]))
    exons <- exons[order(exons$start),]
    if(nrow(exons)!=2 || any(models$V1!="chr38") || any(models$V7!="+")) stop("candidate model does not satisfy the contract for this example")
    write_tsv(exons,file.path(out,"Figure2E.exons.tsv"))

    all_summary <- summary[summary$set=="all",]
    all_summary <- all_summary[match(TISSUES,all_summary$tissue),]
    common_ps <- unique(common_sites$PS)
    common_phase_caption <- paste("Original PS value(s) for the common exonic sites:",
        paste(common_ps,collapse=", "),"(before upstream orientation).")
    common_allele_caption <- if(all(common_sites$HAP1_ALLELE==common_sites$ALT &
                                   common_sites$HAP2_ALLELE==common_sites$REF))
        "At all common exonic sites, HAP1 corresponds to ALT and HAP2 to REF." else
        "REF/ALT correspondence varies or is not uniform across common sites; see the site-level allele columns."
    writeLines(c("Figure E. Variant-level haplotype support at DOG5_NEWGENE_03504 in Dozer.",
        "Positions use CanFam3.1 genomic coordinates; exon lengths and the intron are not compressed.",
        "All retained phased sites are shown; circles: exon union; triangles: intron. No site-depth display cutoff.",
        "Point size shows site support = HAP1_COUNT + HAP2_COUNT, not independent molecules or numbers of sites.",
        "No numeric size legend is displayed; point sizes remain continuous rather than binned. Exact support counts are provided in the site tables.",
        sprintf("For each symbol shape, area scales continuously with support using a common scale across all four tissues: symbol size = max(%.1f, %.1f * sqrt(support / %g)) mm. The minimum visible size departs from strict area proportionality; no sites are removed by this display rule.",
            POINT_MIN_MM,POINT_MAX_MM,max(sites$support)),
        "All: all retained gene-span sites. Common: exonic sites observed in all four tissues.",
        paste("All retained site counts:",paste(paste(LABELS,all_summary$n_sites,sep="="),collapse=", "),"(site counts, not biological replicates)."),
        paste("Common set:",length(common),"exonic variants present in all four tissues, matched by CHROM/POS/REF/ALT."),
        "HAP1 fractions are ratios of summed support, not averages of site fractions. No significance tests are performed.",
        "Dashed lines mark a HAP1 fraction of 0.5; all tissues use the upstream-oriented HAP1 definition.",
        "The common-site comparison fixes the variant set, but site-support weights can differ between tissues.",
        "Original PS blocks were oriented against the SHAPEIT5 scaffold upstream; the counts are not flipped here.",
        common_phase_caption,
        common_allele_caption,
        "The common exonic variants lie in exon 2 and do not validate exon 1 or its splice junction.",
        "Symbols at nearby genomic positions may overlap; positions are not jittered.",
        "One dog is shown. Sites are not biological replicates, and overlapping reads can contribute to multiple sites.",
        "This panel shows site-level detail for a locus included in panel F; it is not an independent validation.",
        "Retained by current CNV rules does not establish normal diploidy. Mapping bias and WGS genotype quality need separate review.",
        "Layout reference only: https://github.com/weishwu/ASEplot/blob/main/R/snp_gene_ase_scatter.R",
        "Original R grid implementation; no exon-coordinate compression or site-median summary is reused."),
        file.path(out,"Figure2E.caption.txt"))
    writeLines(c(paste("GTF",GTF),paste("PHASE_QC",PHASE_QC),paste("R",R.version.string),
        sprintf("Dimensions: %d x %d mm; text: 6.5-10 pt; font: Arial; PNG/TIFF: 600 dpi; PDF: cairo.",WIDTH_MM,HEIGHT_MM),
        paste("Arial regular:",font_paths[1]),paste("Arial bold:",font_paths[2]),
        "Layout: genomic tracks x=17..59 mm; height=6 mm; bottom=31,23,15,7 mm; All/Common x=68.5/80 mm.",
        "Labels: tissue x=8 mm; y-axis ticks right=15 mm; summary title y=43 mm; column titles y=39.5 mm.",
        "Rendered font, collision and alignment QA remains separate; execution completion is not publication approval."),
        file.path(out,"Figure2E.parameters.txt"))
    log_stage("stage 4/4 export 87 × 60 mm PDF, SVG, PNG and RGB TIFF; render after still requires manual check")
    for(ext in c("pdf","svg","png","tiff")) {
        dest <- file.path(out,paste0("Figure2E_candidate_locus_haplotype_support.",ext))
        if(ext=="pdf") grDevices::cairo_pdf(dest,width=WIDTH_MM/25.4,height=HEIGHT_MM/25.4,family=FONT,pointsize=10,bg="white")
        else if(ext=="svg") svglite::svglite(dest,width=WIDTH_MM/25.4,height=HEIGHT_MM/25.4,pointsize=10,bg="white")
        else if(ext=="png") grDevices::png(dest,width=WIDTH_MM,height=HEIGHT_MM,units="mm",res=DPI,type="cairo",pointsize=10,bg="white")
        else grDevices::tiff(dest,width=WIDTH_MM,height=HEIGHT_MM,units="mm",res=DPI,type="cairo",pointsize=10,bg="white",compression="lzw")
        tryCatch(draw_figure(sites,exons,summary,out),finally=grDevices::dev.off())
    }
}

run <- function() {
    run_id <- paste0(format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid())
    out <- file.path(ROOT,run_id)
    if(!dir.create(out,recursive=FALSE)) stop("cannot create new of run directory: ",out)
    con <- file(file.path(out,"plot_figure2E.log"),"wt")
    sink(con,type="output"); sink(con,type="message")
    on.exit({sink(type="message"); sink(type="output"); close(con)})
    tryCatch({
        log_stage(paste("start",out))
        main(out)
        writeLines("execution_status=success",file.path(out,"execution_completion.marker"))
        log_stage("task completed / done execution_status=success")
        log_stage("please separately check count comparison, phase evidence, font, collision, completeness and biological interpretation; completion marker only indicates execution finished.")
    },error=function(e) {
        log_stage(paste("execution_status=failed",conditionMessage(e)))
        stop(e)
    })
}
run()
