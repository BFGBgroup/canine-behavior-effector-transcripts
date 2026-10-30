#!/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/Rscript


Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
  NUMEXPR_NUM_THREADS="1",VECLIB_MAXIMUM_THREADS="1")
ROOT <- "/path/to/project/analysis_20260817"
OUT_ROOT <- file.path(ROOT,"article/figure/figure5/D_exploratory_positional_context")
GWAS_DIR <- "/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.10_run_dog_behavior_susie_rss_finemapping/20260802_205352/main/exploratory/score_separation_related_problems/DOG_GWAS_EXPLORATORY__score_separation_related_problems__chr4_001"
LINK_FILE <- file.path(ROOT,"results/5-dogs/5.allelome-compare-gwas-finemap/summary/20260910_182124/5dogs_gwas_finemap_allelome_main_focus.tsv")
GTF <- file.path(ROOT,"reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/GCF_000002285.3_CanFam3.1_genomic_chr.plus_five_only_new_gene.gtf")
CANDIDATE <- "DOG5_NEWGENE_00538"
PARTNER <- "gene-LOC611908"
FOCUS <- "Affx-206206334"
CHROM <- "chr4"
FONT <- "Arial"
WIDTH_MM <- 87
HEIGHT_MM <- 100
BLUE <- "#0072B2"
GREEN <- "#009E73"

LINK_COLOURS <- c(enhancing="#478577",repressing="#BC220D")
read_tsv <- function(p) read.delim(p,quote="\"",comment.char="",stringsAsFactors=FALSE,
  check.names=FALSE,encoding="UTF-8")
write_tsv <- function(d,p) write.table(d,p,sep="\t",quote=TRUE,row.names=FALSE,
  na="NA",fileEncoding="UTF-8",qmethod="double")
log_stage <- function(s) {cat(format(Sys.time(),"[%Y-%m-%d %H:%M:%S ]"),s,"\n");flush.console()}
attr_value <- function(x,key) {
  hits <- regmatches(x,regexec(paste0('(?:^|;)[[:space:]]*',key,' "([^"]+)"'),x,perl=TRUE))
  vapply(hits,function(a) if(length(a)>1) a[2] else NA_character_,character(1))
}
read_exons <- function() {
  con <- file(GTF,"r"); on.exit(close(con)); chunks <- list()
  repeat {
    x <- readLines(con,n=50000,warn=FALSE)
    if(!length(x)) break
    x <- x[grepl(paste0("^",CHROM,"\t[^\t]+\texon\t"),x)]
    if(!length(x)) next
    d <- read.delim(text=paste(x,collapse="\n"),header=FALSE,quote="",comment.char="")
    chunks[[length(chunks)+1]] <- data.frame(gene_id=attr_value(d$V9,"gene_id"),
      gene_name=attr_value(d$V9,"gene_name"),transcript_id=attr_value(d$V9,"transcript_id"),
      start=d$V4,end=d$V5,strand=d$V7)
  }
  unique(do.call(rbind,chunks))
}
spans <- function(e,key) do.call(rbind,lapply(split(e,e[[key]]),function(z)
  data.frame(id=z[[key]][1],gene_id=z$gene_id[1],gene_name=z$gene_name[1],
    start=min(z$start),end=max(z$end),strand=z$strand[1])))

exon_union <- function(e) {
  do.call(rbind,lapply(split(e,e$gene_id),function(z) {
    z <- z[order(z$start,z$end),]
    out <- list(); a <- z$start[1]; b <- z$end[1]
    emit <- function(a,b) data.frame(gene_id=z$gene_id[1],gene_name=z$gene_name[1],
      start=a,end=b,strand=z$strand[1])
    if(nrow(z)>1) for(i in 2:nrow(z)) {
      if(z$start[i]<=b+1) b <- max(b,z$end[i]) else {
        out[[length(out)+1]] <- emit(a,b); a <- z$start[i]; b <- z$end[i]
      }
    }
    out[[length(out)+1]] <- emit(a,b)
    do.call(rbind,out)
  }))
}
canvas <- function(w,ylim) plot(NA,xlim=w/1e6,ylim=ylim,xaxs="i",yaxs="i",
  axes=FALSE,xlab="",ylab="",bty="n")
coord_axis <- function(w) {
  ticks <- pretty(w/1e6,n=3); ticks <- ticks[ticks>=w[1]/1e6 & ticks<=w[2]/1e6]
  segments(w[1]/1e6,par("usr")[3],w[2]/1e6,par("usr")[3],lwd=0.6,col="#444444")
  axis(1,at=ticks,labels=format(ticks,trim=TRUE),lwd=0,lwd.ticks=0.5,tcl=-0.2)
  mtext("Chr4 position (Mb)",1,line=1.6,cex=7.5/7)
}
structure_row <- function(z,e,key,w,y,label) {
  col <- if(z$gene_id==CANDIDATE) GREEN else "#555555"
  a <- max(z$start,w[1]); b <- min(z$end,w[2])
  segments(a/1e6,y,b/1e6,y,col=col,lwd=0.7)
  ee <- e[e[[key]]==z$id & e$end>=w[1] & e$start<=w[2],,drop=FALSE]
  rect(pmax(ee$start,w[1])/1e6,y-0.09,pmin(ee$end,w[2])/1e6,y+0.09,col=col,border=NA)
  if(b-a>diff(w)*0.01) {
    mid <- mean(c(a,b)); dx <- min((b-a)/4,diff(w)*0.012)
    if(z$strand=="-") dx <- -dx
    arrows((mid-dx)/1e6,y,(mid+dx)/1e6,y,length=0.04,lwd=0.7,
      col=if(z$gene_id==CANDIDATE) "white" else col)
  }
  if(z$start<w[1]) text(w[1]/1e6,y,"/",cex=6.5/7,xpd=NA)
  if(z$end>w[2]) text(w[2]/1e6,y,"/",cex=6.5/7,xpd=NA)

  label_cex <- 6.5/7
  half <- strwidth(label,cex=label_cex)/2; pad <- diff(w)/1e6*0.008
  lx <- max(w[1]/1e6+half+pad,min(w[2]/1e6-half-pad,mean(c(a,b))/1e6))

  label_height <- strheight("Mg",cex=label_cex)
  label_y <- y+0.09+label_height*0.95
  rect(lx-half-pad/3,label_y-label_height*0.65,
    lx+half+pad/3,label_y+label_height*0.65,col="white",border=NA)
  text(lx,label_y,label,cex=label_cex,col=col)
}
quant_track <- function(gw,w,value,label,ylab,ymax,focus,lead,subtitle=NULL) {
  canvas(w,c(0,ymax)); abline(h=0,col="#DDDDDD")
  if(value=="neg_log10_p") {
    abline(h=5,col="#AAAAAA",lty=2,lwd=0.6)

    text((w[2]-diff(w)*0.012)/1e6,5+ymax*0.14,"P = 1e-5",
      adj=c(1,0.5),cex=6.5/7,col="#666666")
  }
  abline(v=focus/1e6,col="#777777",lty=3,lwd=0.6)
  cols <- ifelse(gw$is_cs,BLUE,"#B5B5B5")

  points(gw$pos/1e6,gw[[value]],pch=16,cex=0.65,col=cols)
  ii <- which(gw$is_focus); points(gw$pos[ii]/1e6,gw[[value]][ii],pch=1,cex=1.15,lwd=1)
  jj <- which(gw$ld_variant_id==lead); points(gw$pos[jj]/1e6,gw[[value]][jj],pch=18,cex=0.95)
  axis(2,at=if(value=="pip") seq(0,0.30,0.10) else c(0,2,4,6),las=1,lwd=0,lwd.ticks=0.5,tcl=-0.2)
  segments(w[1]/1e6,0,w[1]/1e6,ymax,col="#444444",lwd=0.6)
  mtext(ylab,2,line=2.6,cex=7.5/7)
  mtext(label,3,line=0.4,adj=0,font=1,cex=8/7)
  if(!is.null(subtitle)) mtext(subtitle,3,line=0.4,adj=1,cex=6.5/7,col="#555555")
}
draw_figure <- function(gw,info,e_union,overview_genes,local_genes,overview,zoom,link_annotation) {


  layout(matrix(1:4,ncol=1),heights=c(20,20,23,37))
  par(family=FONT,ps=7,cex=1,cex.axis=1,omi=c(0,0,4,0)/25.4,mgp=c(1.7,0.45,0))
  focus <- gw$pos[gw$is_focus]
  par(mai=c(1.5,13,5.5,2)/25.4)
  quant_track(gw,zoom,"neg_log10_p","Separation-related problems",
    "-log10(P)",6.5,focus,info$gwas_lead_variant_id,subtitle="GWAS / exploratory")
  quant_track(gw,zoom,"pip","GWAS fine-mapping",
    "PIP",0.30,focus,info$gwas_lead_variant_id,subtitle=paste0(sum(gw$is_cs)," CS variants"))

  f <- gw[gw$is_focus,]
  focus_label <- sprintf("%s\nPIP = %.3f",f$ld_variant_id,f$pip)
  label_width <- max(strwidth(strsplit(focus_label,"\n")[[1]],cex=6.5/7))
  label_x <- max((focus-diff(zoom)*0.015)/1e6,zoom[1]/1e6+label_width)
  text(label_x,0.275,focus_label,adj=c(1,1),cex=6.5/7,col="#444444")
  par(mai=c(9,13,4,2)/25.4); canvas(zoom,c(0,nrow(local_genes)+0.5))
  abline(v=focus/1e6,col="#777777",lty=3,lwd=0.6)
  for(i in seq_len(nrow(local_genes))) {

    structure_row(local_genes[i,],e_union,"gene_id",zoom,nrow(local_genes)-i+0.35,
      sub("^DOG5_NEWGENE_","",local_genes$gene_name[i]))
  }
  coord_axis(zoom); mtext("Gene exon unions",3,line=0.3,adj=0,cex=1)


  par(mai=c(9,13,4,2)/25.4); canvas(overview,c(0,6.2))
  rect(zoom[1]/1e6,0.1,zoom[2]/1e6,3.5,col="#EAF2F8",border="#9DB8CB",lwd=0.6)
  for(i in seq_len(nrow(overview_genes))) structure_row(overview_genes[i,],e_union,"gene_id",overview,
    3.0-(i-1)*0.95,sub("^DOG5_NEWGENE_","",overview_genes$gene_name[i]))
  centers <- vapply(c(CANDIDATE,PARTNER),function(id) {
    z<-overview_genes[overview_genes$id==id,];mean(c(z$start,z$end))/1e6
  },numeric(1))
  lines(c(centers[1],centers[1],centers[2],centers[2]),c(4.5,4.75,4.75,4.5),col=link_annotation$colour,lwd=0.7)
  text(mean(centers),5.55,link_annotation$label,cex=6.5/7,col=link_annotation$colour)
  coord_axis(overview); mtext("Gene context (shaded: local view)",3,line=0.3,adj=0,cex=1)

  grid::grid.text("D",x=grid::unit(2,"mm"),y=grid::unit(1,"npc")-grid::unit(1.5,"mm"),
    just=c("left","top"),gp=grid::gpar(fontfamily=FONT,fontsize=10,fontface="bold"))

}
main <- function() {
  out <- file.path(OUT_ROOT,format(Sys.time(),"%Y%m%d_%H%M%S"))
  if(!dir.create(out,showWarnings=FALSE)) stop("cannot create output directory: ",out)
  con <- file(file.path(out,"plot_figure5D.log"),"wt"); sink(con); sink(con,type="message")
  on.exit({sink(type="message");sink();close(con)},add=TRUE)
  tryCatch({

    script_arg <- grep("^--file=",commandArgs(trailingOnly=FALSE),value=TRUE)
    script_path <- normalizePath(sub("^--file=","",script_arg),mustWork=TRUE)
    snapshot <- file.path(out,"plot_figure5D.executed.R")
    if(!file.copy(script_path,snapshot,overwrite=FALSE)) stop("cannot save this run script copy")
    write_tsv(data.frame(source_path=script_path,snapshot_path=snapshot,
      md5=unname(tools::md5sum(snapshot))),file.path(out,"script_checksum.tsv"))
    log_stage("[1/4] read GWAS and specified pair; concurrent_tasks=1 threads_per_task=1 (default)")
    font_paths <- systemfonts::match_fonts(c(FONT,FONT),weight=c("normal","bold"))$path
    expected_fonts <- c("/path/to/fonts/Arial/Arial.TTF",
                        "/path/to/fonts/Arial/Arialbd.TTF")
    if(!all(file.exists(expected_fonts)) ||
       !identical(normalizePath(font_paths),normalizePath(expected_fonts)))
      stop("Arial regular or bold font is unavailable; do not silently substitute fonts.")
    log_stage(paste("Arial font",paste(font_paths,collapse=";")))
    gw <- read_tsv(file.path(GWAS_DIR,"variant_pip.tsv.gz"))
    info <- read_tsv(file.path(GWAS_DIR,"locus_summary.tsv"))
    cs <- read_tsv(file.path(GWAS_DIR,"credible_sets.tsv"))
    links <- read_tsv(LINK_FILE)
    links <- links[links$manual_new_gene==CANDIDATE & links$protein_coding_gene_id==PARTNER &
      links$gwas_locus_id==info$locus_id & startsWith(links$gwas_variant_display,paste0(FOCUS," |")),,drop=FALSE]
    gw$is_cs <- gw$in_credible_set=="yes"; gw$is_focus <- gw$ld_variant_id==FOCUS
    gw$neg_log10_p <- -log10(gw$p_lrt)
    log_stage("[2/4] read chr4 actual exon; definition overview and local window")
    e <- read_exons(); genes <- spans(e,"gene_id"); all_tx <- spans(e,"transcript_id")
    target <- genes[genes$id==CANDIDATE,]; partner <- genes[genes$id==PARTNER,]
    zoom <- range(c(gw$pos,target$start,target$end))+c(-3000,3000)
    tx <- all_tx[all_tx$start<=zoom[2] & all_tx$end>=zoom[1],]
    tx <- tx[order(tx$gene_id!=CANDIDATE,tx$gene_name,tx$id),]
    overview_genes <- genes[genes$id %in% c(CANDIDATE,PARTNER,tx$gene_id),]
    overview_genes <- overview_genes[order(overview_genes$start),]
    overview <- range(c(zoom,overview_genes$start,overview_genes$end))+c(-15000,15000)
    e <- e[e$gene_id %in% overview_genes$id,]
    e_union <- exon_union(e)
    local_genes <- genes[genes$id %in% tx$gene_id,]
    local_genes <- local_genes[order(local_genes$id!=CANDIDATE,local_genes$start),]

    candidate_exons <- e[e$gene_id==CANDIDATE,]
    gw$within_candidate_span <- gw$pos>=target$start & gw$pos<=target$end
    gw$within_candidate_exon <- vapply(gw$pos,function(p) any(p>=candidate_exons$start & p<=candidate_exons$end),logical(1))

    exon_overlap <- do.call(rbind,lapply(seq_len(nrow(gw)),function(i) {
      hit <- e[e$start<=gw$pos[i] & e$end>=gw$pos[i],,drop=FALSE]
      n <- nrow(hit)
      data.frame(variant_id=rep(gw$ld_variant_id[i],n),chrom=rep(CHROM,n),
        position=rep(gw$pos[i],n),pip=rep(gw$pip[i],n),is_cs=rep(gw$is_cs[i],n),
        gene_id=hit$gene_id,gene_name=hit$gene_name,transcript_id=hit$transcript_id,
        exon_start=hit$start,exon_end=hit$end,strand=hit$strand)
    }))
    tissue_map <- c(cerebellum="CB",cortex_frontal="FC",cortex_occipital="OC")
    sample_parts <- strsplit(unique(links$sample_tissue),"|",fixed=TRUE)
    sample_labels <- vapply(sample_parts,function(s) paste(sub("_[0-9]+$","",s[1]),tissue_map[s[2]],sep="/"),character(1))

    mechanisms <- unique(links$mechanism)
    if(length(mechanisms)==1L && !is.na(mechanisms) && mechanisms %in% names(LINK_COLOURS)) {
      link_annotation <- list(label=if(mechanisms=="repressing") "Predicted repressing LINK" else "Predicted enhancing LINK",
        colour=unname(LINK_COLOURS[mechanisms]))
    } else {
      link_annotation <- list(label=if(nrow(links)) "Reported LINK" else "No LINK record",
        colour="#777777")
    }
    log_stage("[3/4] export overview plus enlarge region figure")
    height <- HEIGHT_MM
    for(fmt in c("pdf","svg","png","tiff")) {
      path <- file.path(out,paste0("Figure5D_exploratory_positional_context.",fmt))
      if(fmt=="pdf") grDevices::cairo_pdf(path,width=WIDTH_MM/25.4,height=height/25.4,pointsize=7,family=FONT)
      if(fmt=="svg") svglite::svglite(path,width=WIDTH_MM/25.4,height=height/25.4,pointsize=7,bg="white")
      if(fmt=="png") grDevices::png(path,width=WIDTH_MM,height=height,units="mm",res=600,type="cairo",pointsize=7,bg="white")
      if(fmt=="tiff") grDevices::tiff(path,width=WIDTH_MM,height=height,units="mm",res=600,type="cairo",pointsize=7,compression="lzw",bg="white")
      tryCatch(draw_figure(gw,info,e_union,overview_genes,local_genes,overview,zoom,link_annotation),finally=grDevices::dev.off())
    }
    log_stage(paste("width mm=",WIDTH_MM,"height mm=",height))
    writeLines(c("width_mm=87; height_mm=100; intended_placement_width_mm=87; intended_scale=1",
      "font=Arial; panel_pt=10; track_title_pt=8; axis_title_pt=7.5; tick_pt=7; dense_labels_pt=6.5; internal_legend=none; shared_legend_font_pt=7",
      "PNG/TIFF=600ppi; TIFF=LZW; PDF_font_embedding=NOT_CHECKED; visual_layout=NOT_CHECKED",
      "assembled_Figure5_and_caption=NOT_CHECKED; overview and local coordinate scales differ",
      "layout_weights=20,20,23,37; total=100; first_two_tracks_aligned_with_current_C",
      "Use Figure5CD_shared_legend.pdf (178 x 8 mm) generated by plot_figure5C.R below the C/D pair; keep full width.",
      "Shared legend adds 8 mm plus spacing; this script neither reads nor modifies the C outputs."),
      file.path(out,"Figure5D_layout.txt"))
    log_stage("[4/4] save source table, figure legend, source and run record")
    tables <- list(gwas_variants=gw,locus_summary=info,credible_sets=cs,LINK_sources=links,
      overview_genes=overview_genes,local_genes=local_genes,local_transcripts=tx,exons=e,
      gene_exon_union=e_union,variant_exon_overlap=exon_overlap)
    for(n in names(tables)) write_tsv(tables[[n]],file.path(out,paste0("Figure5D_",n,".tsv")))
    write_tsv(data.frame(view=c("overview","local"),chrom=CHROM,start=c(overview[1],zoom[1]),
      end=c(overview[2],zoom[2]),coordinate_system="1-based inclusive"),file.path(out,"Figure5D_windows.tsv"))
    inputs <- c(gwas=file.path(GWAS_DIR,"variant_pip.tsv.gz"),summary=file.path(GWAS_DIR,"locus_summary.tsv"),
      cs=file.path(GWAS_DIR,"credible_sets.tsv"),links=LINK_FILE,gtf=GTF)
    write_tsv(data.frame(role=names(inputs),absolute_path=unname(inputs)),file.path(out,"Figure5D_input_paths.tsv"))
    f <- gw[gw$is_focus,]

    reference_hits <- unique(exon_overlap[exon_overlap$is_cs & exon_overlap$gene_id!=CANDIDATE,
      c("variant_id","gene_name","pip"),drop=FALSE])
    reference_note <- if(nrow(reference_hits)) paste0(
      "Other CS exon overlaps among the displayed genes: ",
      paste(sprintf("%s in %s (PIP %.6f)",reference_hits$variant_id,
        reference_hits$gene_name,reference_hits$pip),collapse="; "),
      ". The added candidate is not the only possible regional interpretation; exon overlap does not specify a coding consequence.") else character(0)
    cs_notes <- sprintf("CS %s: %d of %d modeled variants; upstream reference-LD absolute correlation min %.6f, mean %.6f. These are imported summaries, not LD recalculated by this plotting script.",
      cs$cs_id,cs$cs_size,nrow(gw),cs$min_abs_corr,cs$mean_abs_corr)
    writeLines(c("D. Exploratory positional overlap at DOG5_NEWGENE_00538 (CanFam3.1).",
      paste("Separation-related problems; modeled variants:",nrow(gw),"; GWAS CS members:",sum(gw$is_cs)),
      paste("Regional lead:",info$gwas_lead_variant_id,"; P:",info$gwas_lead_p_lrt,"; highest-PIP variant:",info$top_variant_id),
      paste("Candidate-overlapping variant:",FOCUS,"; position:",f$pos,"; P:",f$p_lrt,"; PIP:",f$pip,
        "; candidate exon overlap:",f$within_candidate_exon),
      "The bottom overview shows selected genes at true coordinates; the shaded interval corresponds to the enlarged association, PIP and gene tracks above. The overview uses a separate coordinate scale; the SNP guide line is confined to the local tracks. The five-digit label 00538 in both gene views denotes DOG5_NEWGENE_00538, with the leading zeros retained; this abbreviation follows panel C. Full identifiers remain in the source tables. Both annotation views show gene-level exon unions, not individual isoforms. Original transcript IDs and exon coordinates remain in the source tables; / denotes a clipped end. Annotation is not RNA coverage.",
      "This panel uses the separate C/D shared legend generated by the C plotting script. Selected variant denotes Affx-206206334 in D; Molecular CS applies only to C. No molecular posterior is displayed in D.",
      "Blue indicates CS membership, grey other modeled variants; circles mark the candidate-overlapping SNP and diamonds the regional GWAS lead. The P=1e-5 line is a reference only. The region remains exploratory; no threshold is relaxed by plotting it.",
      "P values are shown only for variants in this fine-mapping model, not all original GWAS variants. PIP is unchanged and not renormalized. No molecular PPA is used in this evidence path.",
      reference_note,cs_notes,
      "The highly correlated CS members are not clearly resolved into a unique causal variant. Variants outside the modeled input have not been excluded by this analysis.",
      paste("LINK source samples:",paste(sample_labels,collapse=", ")),
      paste("Bracket label:",link_annotation$label,
        "; classifications are predictions from Allelome.LINK, not experimentally established regulation."),
      "The undirected bracket connects the candidate and its coding partner, not the highlighted SNP. Enhancing (#478577) and repressing (#BC220D) follow the author colours also used in panel C.",
      "These are tissue observations, not independent dogs. LINK modes and scores are retained in the source table without assigning causal direction. No multigene flag in these records is not proof of functional regulation.",
      "Positional overlap, CS membership and gene-level LINK do not establish a causal variant, colocalization or a variant-to-behaviour mechanism."
    ),file.path(out,"Figure5D_caption.txt"))
    writeLines(c("Schwartzentruber et al. 2021, doi:10.1038/s41588-020-00776-w, Figure 4.",
      "Verified author master: https://github.com/Jeremy37/AD_finemap/blob/a15ca4e0d32f83faf92bb2cc1c918b67e3e1ac42/src/AD_finemap_master.sh",
      "Verified figure script: https://github.com/Jeremy37/AD_finemap/blob/a15ca4e0d32f83faf92bb2cc1c918b67e3e1ac42/src/finemapping_plots.R",
      "Read-only retrieval confirmed Figure 4 mapping, cowplot vertical alignment, probability stems and rescale_introns=F. Author code was not executed. This is an independent base-R implementation; no human coordinates, LD colours, ATAC or motifs are imported.",
      "Allelome.LINK mechanism colours: https://github.com/AndergassenLab/Allelome.LINK/blob/main/00_src/Allelome.LINK.R ; doi:10.1093/nar/gkaf912",
      "Official track configuration reference: https://pygenometracks.readthedocs.io/en/latest/content/examples.html#basic-examples"
    ),file.path(out,"Figure5D_references.txt"))
    capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
    writeLines(c("execution_status=success",paste0("output_dir=",out)),file.path(out,"execution_completion.marker"))
    log_stage("task completed / done; execution_status=success")
    log_stage("please check the figure layout separately, coordinates, input span with biological interpretation; execution success does not imply quality audit passed.")
  },error=function(e) {log_stage(paste("execution_status=failed;",conditionMessage(e)));stop(e)})
}
main()
