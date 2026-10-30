#!/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/Rscript


Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           NUMEXPR_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
ROOT <- "/path/to/repository/scripts/13_figures/historical_figure_numbers/figure2/F_candidate_haplotype_imbalance"
MANIFEST <- "/path/to/repository/scripts/08_expression_ase_integration/preparation/1.5ase-cnv-clean-gene-id-haplotype-table/log/1.5.ase_cnv_clean_gene_id_haplotype_table.20260821_161805.success.tsv"
MAP <- "/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.transcript_map.tsv"
SOURCE_URL <- "https://github.com/JianpingQuan/ASEanalysis/blob/fc5c1229d16d024c2bcdfeae9e44ceca37459750/RNA-Seq/R/Figure%203e.R"
DOGS <- c("Ace_396226", "CodRou_346627", "Crak_358336", "Dozer_386519", "SanMcg_349891")
TISSUES <- c("cerebellum", "cortex_frontal", "cortex_occipital", "pituitary_gland")
MIN_SUPPORT <- 20
TOP_N <- 12L
WIDTH_MM <- 87
FONT <- "Arial"
DPI <- 600

args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 2L) stop("usage: plot_figure2F.R [min_support] [top_n]")
if (length(args) >= 1L) MIN_SUPPORT <- as.numeric(args[1])
if (length(args) >= 2L) TOP_N <- as.numeric(args[2])
if (!is.finite(MIN_SUPPORT) || MIN_SUPPORT < 1 || MIN_SUPPORT != floor(MIN_SUPPORT) ||
    !is.finite(TOP_N) || TOP_N < 1 || TOP_N > 30 || TOP_N != floor(TOP_N))
    stop("min_support must be positive integer, top_n must be 1 to 30 integer")
read_tsv <- function(p) read.delim(p, colClasses="character", quote="",
    comment.char="", na.strings=character(), check.names=FALSE)
write_tsv <- function(x,p) write.table(x,p,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
log_stage <- function(s) {cat(format(Sys.time(),"[%Y-%m-%d %H:%M:%S]"),s,"\n"); flush.console()}
collapse_values <- function(x) paste(sort(unique(x[!x %in% c("", "NA", "-")])),collapse=";")
safe_median <- function(x) if (length(x)) median(x) else NA_real_

configure_font <- function() {
    paths <- systemfonts::match_fonts(c(FONT,FONT),weight=c("normal","bold"))$path
    expected <- c("/path/to/fonts/Arial/Arial.TTF",
                  "/path/to/fonts/Arial/Arialbd.TTF")
    if(!all(file.exists(expected)) || !identical(normalizePath(paths),normalizePath(expected)))
        stop("Arial regular/bold parse does not match, please check font installation")
    paths
}

draw_heatmap <- function(m, status, selected, entries, height_mm) {

    grid::grid.newpage()
    gp <- function(size=7,face="plain",col="black")
        grid::gpar(fontfamily=FONT,fontsize=size,fontface=face,col=col)
    txt <- function(s,x,y,size=7,just="centre",face="plain",rot=0)
        grid::grid.text(s,grid::unit(x,"mm"),grid::unit(y,"mm"),just=just,rot=rot,gp=gp(size,face))
    grid::grid.rect(gp=grid::gpar(fill="white",col=NA))
    txt("F",2,height_mm-2,10,c("left","top"),"bold")
    n <- nrow(m); top <- height_mm-16; bottom <- 14; h <- (top-bottom)/n
    yy <- top-(seq_len(n)-0.5)*h
    dog_ix <- match(entries$animal_id,DOGS)

    xx <- 12+(seq_len(ncol(m))-0.5)*3.65+1.0*(dog_ix-1)
    pal <- grDevices::colorRampPalette(c("#FFF5F0","#FCAE91","#CB181D","#67000D"))(100)
    for(j in seq_len(ncol(m))) {
        v <- m[,j]; fills <- rep("#D9D9D9",n)
        ok <- status[,j]=="shown"
        fills[ok] <- pal[pmax(1,pmin(100,1+floor((v[ok]-0.5)*198)))]
        fills[status[,j]=="low_support"] <- "#FAFAFA"
        grid::grid.rect(x=grid::unit(xx[j],"mm"),y=grid::unit(yy,"mm"),
            width=grid::unit(3.5,"mm"),height=grid::unit(h*0.97,"mm"),
            gp=grid::gpar(fill=fills,col="white",lwd=0.35))
        low <- status[,j]=="low_support"
        if(any(low)) grid::grid.points(x=grid::unit(rep(xx[j],sum(low)),"mm"),
            y=grid::unit(yy[low],"mm"),pch=4,size=grid::unit(1.2,"mm"),gp=gp(6.5,col="#777777"))

        tissue_index <- match(entries$tissue[j],TISSUES)
        txt(c("CB","FC","OC","Pit")[tissue_index],xx[j],top+2.6,6.5)
    }
    for(d in seq_along(DOGS)) txt(sub("_.*","",DOGS[d]),mean(xx[dog_ix==d]),top+6.4,7,face="bold")

    txt("DOG5_NEWGENE_",10,height_mm-3.5,7,"left")
    txt(sub("^DOG5_NEWGENE_","",selected$gene_id),10.5,yy,6.5,"right")
    legend_center <- (min(xx)+max(xx))/2
    grid::grid.rect(x=grid::unit(seq(legend_center-20,legend_center+20,length.out=100),"mm"),
        y=grid::unit(9,"mm"),width=grid::unit(0.41,"mm"),height=grid::unit(1.2,"mm"),
        gp=grid::gpar(fill=pal,col=NA))
    txt(c("0.5","0.75","1.0"),legend_center+c(-20,0,20),6.7,6.5)
    txt("Major haplotype support fraction",legend_center,3.5,7)

    if(any(status=="low_support"))
        txt(paste0("x  Support < ",MIN_SUPPORT),22,11.8,6.5,"left")
    if(any(status=="not_evaluable")) {
        grid::grid.rect(x=grid::unit(59,"mm"),y=grid::unit(11.8,"mm"),
            width=grid::unit(1.5,"mm"),height=grid::unit(1.5,"mm"),gp=grid::gpar(fill="#D9D9D9",col=NA))
        txt("Not evaluable",61,11.8,6.5,"left")
    }
}

main <- function(run) {
    font_paths <- configure_font()
    log_stage("start; concurrent_tasks=1; threads_per_task=1; main count retain all site")
    if (!requireNamespace("svglite",quietly=TRUE)) stop("requires svglite package")
    entries <- read_tsv(MANIFEST)
    if (any(!entries$animal_id %in% DOGS) || any(!entries$tissue %in% TISSUES) ||
        nrow(entries)!=19 || anyDuplicated(entries[c("animal_id","tissue")]))
        stop("sample list does not satisfy 19 dogs tissue combination")
    entries <- entries[order(match(entries$animal_id,DOGS),match(entries$tissue,TISSUES)),]
    sample_keys <- paste(entries$animal_id,entries$tissue,sep="|")
    if (any(entries$animal_id==DOGS[3] & entries$tissue==TISSUES[4])) stop("Crak pituitary does not belong to this run sample")
    mp <- read_tsv(MAP)

    meta <- do.call(rbind,lapply(split(mp,mp$final_gene_id),function(z) {
        if (length(unique(z$chr))!=1 || length(unique(z$gene_start_0based))!=1)
            stop("candidate mapping exists mismatch of locus coordinates")
        data.frame(gene_id=z$final_gene_id[1],chromosome=z$chr[1],
            start0=as.numeric(z$gene_start_0based[1]),
            candidate_class=collapse_values(z$novelty_category),
            reference_gene_ids=collapse_values(z$reference_gene_id),
            predicted_biotypes=collapse_values(z$FEELnc_biotype))
    }))
    pieces <- vector("list",nrow(entries))
    numeric_fields <- c("variant_count","hap1_count_sum","hap2_count_sum","total_count_sum",
        "unambiguous_variant_count","multi_gene_variant_count","unambiguous_hap1_count_sum",
        "unambiguous_hap2_count_sum","multi_gene_hap1_count_sum","multi_gene_hap2_count_sum")
    for (i in seq_len(nrow(entries))) {
        log_stage(paste("read",i,"/19",sample_keys[i]))
        z <- read_tsv(entries$output_tsv[i]); z <- z[z$record_type=="manual_new_gene",]
        for (nm in numeric_fields) z[[nm]] <- as.numeric(z[[nm]])
        if (anyDuplicated(z$gene_id) || any(!z$gene_id %in% meta$gene_id) ||
            any(z$sample!=entries$animal_id[i] | z$tissue!=entries$tissue[i]))
            stop("candidate ID or sample identity does not satisfy input contract")
        if (any(!is.finite(z$total_count_sum)) || any(z$total_count_sum<0) ||
            any(z$hap1_count_sum+z$hap2_count_sum!=z$total_count_sum))
            stop("haplotype count cannot form a proportion")
        z$sample_key <- sample_keys[i]
        z$observed <- z$variant_count>=1 & z$total_count_sum>0
        z$hap1_fraction <- z$major_fraction <- NA_real_
        z$hap1_fraction[z$observed] <- z$hap1_count_sum[z$observed]/z$total_count_sum[z$observed]
        z$major_fraction[z$observed] <- pmax(z$hap1_fraction[z$observed],1-z$hap1_fraction[z$observed])
        z$display_state <- ifelse(!z$observed,"not_evaluable",
                                 ifelse(z$total_count_sum>=MIN_SUPPORT,"shown","low_support"))
        z$input_gene_tsv <- entries$output_tsv[i]
        pieces[[i]] <- z
    }
    all <- do.call(rbind,pieces)

    m <- matrix(NA_real_,nrow(meta),19,dimnames=list(meta$gene_id,sample_keys))
    state <- matrix("not_evaluable",nrow(meta),19,dimnames=dimnames(m))
    ix <- cbind(match(all$gene_id,meta$gene_id),match(all$sample_key,sample_keys))
    m[ix] <- all$major_fraction; state[ix] <- all$display_state
    stats <- do.call(rbind,lapply(meta$gene_id,function(id) {
        z <- all[all$gene_id==id,]; obs <- z[z$observed,]; shown <- z[z$display_state=="shown",]
        data.frame(gene_id=id,n_dogs_observed=length(unique(obs$sample)),
            n_dogs_meeting_display_criteria=length(unique(shown$sample)),
            n_dogs_with_multiple_evaluable_tissues=sum(table(shown$sample)>=2),
            n_evaluable_combinations=nrow(shown),
            min_observed_support=if(nrow(obs)) min(obs$total_count_sum) else NA_real_,
            median_total_support=safe_median(shown$total_count_sum),
            median_variant_count=safe_median(shown$variant_count))
    }))
    stats <- cbind(stats,meta[match(stats$gene_id,meta$gene_id),setdiff(names(meta),"gene_id")])
    chr_order <- match(stats$chromosome,paste0("chr",1:38))
    rank_order <- order(-stats$n_dogs_meeting_display_criteria,
        -stats$n_dogs_with_multiple_evaluable_tissues,-stats$n_evaluable_combinations,
        -stats$median_total_support,-stats$median_variant_count,chr_order,stats$chromosome,
        stats$start0,stats$gene_id,na.last=TRUE)
    stats <- stats[rank_order,]; stats$selection_rank <- seq_len(nrow(stats))
    stats$selected_for_F <- stats$n_evaluable_combinations>0 & stats$selection_rank<=TOP_N
    stats$selection_reason <- ifelse(stats$selected_for_F,
        "support-ranked subset; no imbalance-based selection",
        ifelse(stats$n_evaluable_combinations==0,"no cells meeting display threshold","outside top_n"))
    selected <- stats[stats$selected_for_F,]
    if (!nrow(selected)) stop("current display threshold under no candidate")

    selected <- selected[order(match(selected$chromosome,paste0("chr",1:38)),
                               selected$chromosome,selected$start0,selected$gene_id),]
    log_stage(paste("complete candidate",nrow(meta),"; display",nrow(selected),"；min_support",MIN_SUPPORT))
    write_tsv(all,file.path(run,"Figure2F.candidate_counts.tsv"))
    write_tsv(stats,file.path(run,"Figure2F.selection.tsv"))
    write_tsv(entries,file.path(run,"Figure2F.input_manifest.tsv"))
    write_tsv(data.frame(gene_id=rownames(m),m,check.names=FALSE),file.path(run,"Figure2F.major_fraction.all.tsv"))
    write_tsv(data.frame(gene_id=rownames(state),state,check.names=FALSE),file.path(run,"Figure2F.display_state.all.tsv"))
    write_tsv(selected,file.path(run,"Figure2F.selected_rows.tsv"))
    mm <- m[selected$gene_id,,drop=FALSE]; ss <- state[selected$gene_id,,drop=FALSE]

    state_caption <- if (all(ss=="shown"))
        paste("All",nrow(selected),"displayed candidates meet the total-support threshold of",
              MIN_SUPPORT,"in all",ncol(ss),"sampled dog-tissue combinations.") else
        c(if (any(ss=="low_support"))
              "Crosses indicate positive counts below the display threshold; they do not indicate absent expression.",
          if (any(ss=="not_evaluable"))
              "Gray cells indicate no evaluable counts, not necessarily absent expression.")
    writeLines(c("F. Gene-span haplotype imbalance at a support-selected subset of candidate added loci.",
        "Colors show max(HAP1,HAP2)/(HAP1+HAP2), fixed at 0.5-1; no row scaling or clustering. Haplotype direction and parental origin are not encoded.",
        "The shared prefix DOG5_NEWGENE_ is shown once on a single line at the top of the panel; for example, 03504 denotes DOG5_NEWGENE_03504. Full candidate IDs and the dog-count, median-variant and median-support summaries are retained in Figure2F.selected_rows.tsv; these summaries are not displayed in the panel.",
        paste("Exploratory display threshold: total support >=",MIN_SUPPORT,"; top_n =",TOP_N),
        state_caption,
        "Crak pituitary was not sampled and has no column.",
        "Selection prioritizes dogs meeting the display threshold, dogs with at least two evaluable tissues, evaluable combinations, median total support, and median variant count, then genomic position and ID. No imbalance statistic is used for selection. In the source table, medians use above-threshold cells only, and dog counts include dogs with at least one above-threshold tissue.",
        "All retained gene-span variant support is used, including intronic and multi-gene assignments. Summed site support is not independent RNA molecule counts, TPM, transcript-specific expression, or a significant ASE call. Candidate loci can include reference-related models.",
        "The complete annotation-universe matrices include unobserved loci; candidate_counts retains all available candidate records, including zero-support records. Missing cells have no imputed counts.",
        "Layout references Quan Fig.3e (fixed matrix/sample order and missingness), and pheatmap draw_matrix/find_coordinates. This is an independent R grid implementation, not a reproduction of parental-origin analysis.",SOURCE_URL),file.path(run,"Figure2F.caption.md"))

    height_mm <- max(60,30+2.5*nrow(selected))
    writeLines(c(paste("manifest",MANIFEST),paste("candidate_map",MAP),
        paste("min_support",MIN_SUPPORT),paste("top_n",TOP_N),
        paste("width_mm",WIDTH_MM),paste("height_mm",height_mm),
        paste("Arial regular",font_paths[1]),paste("Arial bold",font_paths[2]),
        "font_pt: panel=10; dog/prefix/colorbar_title=7; candidate/tissue/ticks=6.5",
        "layout_mm: matrix left=12; right=85.35; bottom=14; top=height-16; column_pitch=3.65; cell_width=3.5; dog_gap=1.0; candidate_right=10.5; tissue_rotation=0; tissue_y=top+2.6; dog_y=top+6.4; prefix_x=10; prefix_y=height-3.5; prefix_single_line=true",
        "raster_dpi 600; TIFF RGB LZW",
        "concurrent_tasks 1", "threads_per_task 1",
        "Rendered PDF text, collision and final-size visual QA remain separate from execution success."),
        file.path(run,"Figure2F.parameters.txt"))
    w <- WIDTH_MM/25.4; h <- height_mm/25.4
    log_stage(paste("export",WIDTH_MM,"×",height_mm,"mm PDF、SVG、600 dpi PNG/RGB TIFF"))
    grDevices::cairo_pdf(file.path(run,"Figure2F_candidate_haplotype_imbalance.pdf"),
        width=w,height=h,family=FONT,pointsize=10,bg="white")
    draw_heatmap(mm,ss,selected,entries,height_mm); grDevices::dev.off()
    svglite::svglite(file.path(run,"Figure2F_candidate_haplotype_imbalance.svg"),width=w,height=h,pointsize=10,bg="white")
    draw_heatmap(mm,ss,selected,entries,height_mm); grDevices::dev.off()
    grDevices::png(file.path(run,"Figure2F_candidate_haplotype_imbalance.png"),
        width=w,height=h,units="in",res=DPI,pointsize=10,bg="white",type="cairo")
    draw_heatmap(mm,ss,selected,entries,height_mm); grDevices::dev.off()
    grDevices::tiff(file.path(run,"Figure2F_candidate_haplotype_imbalance.tiff"),
        width=WIDTH_MM,height=height_mm,units="mm",res=DPI,pointsize=10,bg="white",type="cairo",compression="lzw")
    draw_heatmap(mm,ss,selected,entries,height_mm); grDevices::dev.off()
    writeLines(capture.output(sessionInfo()),file.path(run,"R_session.txt"))
    writeLines("execution_status=success",file.path(run,"execution_completion.marker"))
    log_stage("task completed / done; execution_status=success; please check separately support threshold, select record and final layout")
}
run <- file.path(ROOT,format(Sys.time(),"%Y%m%d_%H%M%S"))
if (dir.exists(run) || !dir.create(run,recursive=TRUE)) stop("cannot create independent output directory")
con <- file(file.path(run,"plot_figure2F.log"),"wt")
sink(con); sink(con,type="message")
exit_status <- tryCatch({main(run);0L},error=function(e) {
    message("execution_status=failed: ",conditionMessage(e));1L
})
sink(type="message"); sink(); close(con)
quit(save="no",status=exit_status)
