#!/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/Rscript


Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           NUMEXPR_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
ROOT <- "/path/to/project/analysis_20260817"
OUT_ROOT <- file.path(ROOT, "article/figure/figure4/D_candidate_molecular_finemapping")

WIDTH_MM <- 178
HEIGHT_MM <- 55
FONT <- "Arial"
BASE_PT <- 8
GENE <- "DOG5_NEWGENE_03157"
BLOCK <- "no_gwas__cortex_frontal__DOG5_NEWGENE_03157__chr30"
INPUT <- file.path(ROOT, "5-dogs/10.finemaping/formal/8.2_no_gwas_prepare_expression_ase_plasma_inputs/20260821_181104/plasma_inputs", BLOCK)
RESULT <- file.path(ROOT, "5-dogs/10.finemaping/formal/8.3_no_gwas_run_expression_ase_plasma_joint/20260822_233245/plasma_results", BLOCK)
GTF <- file.path(ROOT, "reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/GCF_000002285.3_CanFam3.1_genomic_chr.plus_five_only_new_gene.gtf")

read_tsv <- function(p) read.delim(p, check.names=FALSE, stringsAsFactors=FALSE)
attribute <- function(x, key) {
  pattern <- paste0('(?:^|;[ ]*)', key, ' "([^"]+)"')
  hit <- regexec(pattern, x, perl=TRUE)
  vapply(regmatches(x, hit), function(z) if (length(z)>1L) z[2L] else NA_character_, "")
}


read_region <- function(path, chromosome, limits) {
  con <- file(path, "rt"); on.exit(close(con))
  pieces <- list()
  repeat {
    lines <- readLines(con, n=50000L, warn=FALSE)
    if (!length(lines)) break
    lines <- lines[startsWith(lines, paste0(chromosome, "\t"))]
    if (!length(lines)) next
    z <- read.table(text=lines, sep="\t", header=FALSE, quote="", comment.char="",
                    stringsAsFactors=FALSE, col.names=c("chrom","source","feature","start","end","score","strand","frame","attributes"))
    z <- z[z$feature %in% c("transcript", "exon") & z$start<=limits[2] & z$end>=limits[1], ]
    if (nrow(z)) pieces[[length(pieces)+1L]] <- z
  }
  if (!length(pieces)) stop("GTF within window no structure annotation")
  z <- do.call(rbind, pieces)
  z$gene_id <- attribute(z$attributes, "gene_id")
  z$gene_name <- attribute(z$attributes, "gene_name")
  z
}


exon_union <- function(z) {
  z <- unique(z[z$feature=="exon", c("start","end")])
  if (!nrow(z)) return(z)
  z <- z[order(z$start,z$end), , drop=FALSE]
  out <- z[1L, , drop=FALSE]
  if (nrow(z)>1L) for (i in 2:nrow(z)) {
    j <- nrow(out)
    if (z$start[i]<=out$end[j]+1) out$end[j] <- max(out$end[j],z$end[i])
    else out <- rbind(out,z[i, , drop=FALSE])
  }
  out
}

main <- function() {
  stamp <- paste0(format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid())
  out <- file.path(OUT_ROOT,stamp)

  if (!dir.create(out,recursive=TRUE,showWarnings=FALSE)) {
    stop(paste0("this run output directory already exists or cannot create: ",out),call.=FALSE)
  }
  logfile <- file.path(out,"plot_figure4D.log")
  cat("detailed run log: ",logfile,"\n",sep="")
  con <- file(logfile,"wt")
  sink(con); sink(con,type="message")
  on.exit({sink(type="message"); sink(); close(con)},add=TRUE)
  log <- function(s) {cat(format(Sys.time(),"[%Y-%m-%d %H:%M:%S ]"),s,"\n"); flush(con)}
  write_tsv <- function(x,name) write.table(x,file.path(out,name),sep="\t",row.names=FALSE,quote=TRUE,na="NA")
  tryCatch({
    font_paths <- systemfonts::match_fonts(c(FONT, FONT), weight=c("normal", "bold"))$path
    expected_fonts <- c("/path/to/fonts/Arial/Arial.TTF",
                        "/path/to/fonts/Arial/Arialbd.TTF")
    if (!all(file.exists(expected_fonts)) ||
        !identical(normalizePath(font_paths), normalizePath(expected_fonts)))
      stop("Arial regular or bold font is unavailable; please install the specified font, do not silently substitute.")
    log("[1/4] read original block; concurrent_tasks=1 threads_per_task=1")
    cat("INPUT=",INPUT,"\nRESULT=",RESULT,"\nGTF=",GTF,"\nOUTPUT=",out,"\n",sep="")
    m <- read_tsv(file.path(INPUT,"markers.tsv"))
    s <- read_tsv(file.path(INPUT,"samples.tsv"))
    ppa <- scan(file.path(RESULT,"ppas.txt"),quiet=TRUE)
    cs <- scan(file.path(RESULT,"cset.txt"),quiet=TRUE)
    a <- as.matrix(read.table(file.path(INPUT,"hap_A.txt")))
    b <- as.matrix(read.table(file.path(INPUT,"hap_B.txt")))

    if (length(ppa)!=nrow(m) || length(cs)!=nrow(m) ||
        !identical(dim(a),dim(b)) || nrow(a)!=nrow(s) || ncol(a)!=nrow(m)) stop("Input matrix dimensions do not match the sample and marker ordering")
    if (anyDuplicated(m$marker_id) || any(!is.finite(ppa)) || any(ppa<0 | ppa>1) ||
        any(!cs %in% c(0,1)) || length(unique(m$chrom))!=1L) stop("variant identities or probability inputs are not suitable for plotting")
    if (any(!a %in% c(0,1)) || any(!b %in% c(0,1))) stop("haplotype matrix is not 0/1 encoding")
    m$ppa <- ppa; m$credible_set_member <- cs==1
    m$n_heterozygous_dogs <- colSums(a!=b)

    for (i in seq_len(nrow(s))) m[[paste0("GT_",s$sample[i])]] <- paste(a[i,],b[i,],sep="|")
    limits <- range(m$pos)
    log("[2/4] read local GTF and construct gene exon union")
    g <- read_region(GTF,unique(m$chrom),limits)
    if (!GENE %in% g$gene_id) stop("target candidate not in GTF within window")
    ids <- unique(g$gene_id)
    starts <- vapply(ids,function(id) min(g$start[g$gene_id==id]),0)
    ids <- c(GENE,ids[order(starts)][ids[order(starts)]!=GENE])
    ymax <- ceiling(max(m$ppa)*1.15*100)/100
    if (ymax<=0) stop("No positive PPA values are available to define the y-axis")
    height_mm <- HEIGHT_MM
    blue <- "#0072B2"; orange <- "#D55E00"; grey <- "#B6B6B6"


    models <- lapply(ids,function(id) {
      z <- g[g$gene_id==id,]
      labels <- unique(z$gene_name[!is.na(z$gene_name) & nzchar(z$gene_name)])
      strands <- unique(z$strand[z$strand %in% c("+","-")])
      list(id=id, label=if(length(labels)) labels[1] else id,
           start=min(z$start), end=max(z$end),
           strand=if(length(strands)==1L) strands else NA_character_, ex=exon_union(z))
    })
    ppa_height <- 5
    draw <- function() {

      par(mai=c(.50,.46,.16,.08),mgp=c(1.8,.5,0),tcl=-.18,cex=1,
          family=FONT,las=1,xaxs="i",yaxs="i")
      plot.new()
      plot.window(xlim=limits,ylim=c(0,1))

      gap <- diff(limits)*.018
      lane_ends <- numeric(); lane_load <- numeric(); lanes <- integer(length(models))
      label_x <- numeric(length(models))
      for(j in seq_along(models)) {
        z <- models[[j]]
        label <- if(j==1L) GENE else z$label
        w <- strwidth(label,cex=6.5/BASE_PT)
        lo <- max(z$start,limits[1]); hi <- min(z$end,limits[2])
        label_x[j] <- max(limits[1]+w/2,min((lo+hi)/2,limits[2]-w/2))
        left <- min(lo,label_x[j]-w/2); right <- max(hi,label_x[j]+w/2)
        if(j==1L) next
        available <- which(lane_ends+gap<left)

        lane <- if(length(available)) available[which.min(lane_load[available])] else length(lane_ends)+1L
        if(lane>length(lane_load)) lane_load[lane] <- 0
        lane_load[lane] <- lane_load[lane]+right-left
        lanes[j] <- lane; lane_ends[lane] <- right
      }
      nlanes <- length(lane_ends)
      gene_top <- 1.2*nlanes+.75

      ppa_bottom <- gene_top+1.2
      plot.window(xlim=limits,ylim=c(0,ppa_bottom+ppa_height+.5))
      ticks <- pretty(c(0,ymax),n=4); ticks <- ticks[ticks>=0 & ticks<=ymax]
      yy <- ppa_bottom+ticks/ymax*ppa_height
      segments(limits[1],yy,limits[2],yy,col="#F2F2F2",lwd=.4)
      axis(2,at=yy,labels=format(ticks,trim=TRUE),cex.axis=7/BASE_PT)
      mtext("PPA",side=2,line=1.9,at=ppa_bottom+ppa_height/2,las=0,cex=7.5/BASE_PT)
      points(m$pos,ppa_bottom+m$ppa/ymax*ppa_height,pch=16,cex=.38,col=grey)
      k <- m$credible_set_member
      points(m$pos[k],ppa_bottom+m$ppa[k]/ymax*ppa_height,pch=16,cex=.72,col=blue)
      segments(limits[1],ppa_bottom,limits[2],ppa_bottom,lwd=.65)
      for (j in seq_along(models)) {
        z <- models[[j]]; ex <- z$ex
        y <- if(j==1L) gene_top else gene_top-1.2*lanes[j]
        colour <- if(j==1L) orange else "#666666"
        lo <- max(z$start,limits[1]); hi <- min(z$end,limits[2])
        segments(lo,y,hi,y,col=colour,lwd=.45)

        if(!is.na(z$strand) && hi-lo>diff(limits)*.012) {
          mid <- (lo+hi)/2
          dx <- min((hi-lo)*.18,diff(limits)*.008)
          direction <- if(z$strand=="+") 1 else -1
          arrows(mid-direction*dx,y,mid+direction*dx,y,length=.035,
                 angle=25,col=colour,lwd=.5)
        }
        if (nrow(ex)) rect(pmax(ex$start,limits[1]),y-.19,pmin(ex$end,limits[2]),y+.19,col=colour,border=NA)

        dx <- diff(limits)*.0014
        if(z$start<limits[1]) {
          xx <- lo+c(.5,2)*dx
          segments(xx-dx/2,y-.15,xx+dx/2,y+.15,col=colour,lwd=.7)
        }
        if(z$end>limits[2]) {
          xx <- hi-c(.5,2)*dx
          segments(xx-dx/2,y-.15,xx+dx/2,y+.15,col=colour,lwd=.7)
        }
        text(label_x[j],y+.43,if(j==1L) GENE else z$label,
             adj=c(.5,0),cex=6.5/BASE_PT,col=colour)
      }
      xt <- pretty(limits,n=5); xt <- xt[xt>=limits[1] & xt<=limits[2]]
      axis(1,at=xt,labels=format(xt/1e6,nsmall=2,trim=TRUE),cex.axis=7/BASE_PT)
      mtext(paste0(unique(m$chrom)," position (Mb), CanFam3.1"),1,line=1.65,cex=7.5/BASE_PT)


      text(grconvertX(5/72,from="inches",to="user"),
           grconvertY(height_mm/25.4-5/72,from="inches",to="user"),"D",
           adj=c(0,1),font=2,cex=10/BASE_PT,xpd=NA)

      legend(limits[1],ppa_bottom+ppa_height+.45,xjust=0,yjust=1,
             legend=c("Other variants","95% CS members"),pch=16,
             col=c(grey,blue),bty="n",horiz=TRUE,cex=7/BASE_PT,pt.cex=.7)

      text(limits[1],grconvertY(3/25.4,from="inches",to="user"),
           "Original implementation output",adj=c(0,0),
           cex=6.5/BASE_PT,col="#777777",xpd=NA)
    }
    log("[3/4] export original result visualization draft and source table; do not run PLASMA")
    prefix <- file.path(out,"Figure4D_candidate_molecular_finemapping_draft")

    for (ext in c("pdf","svg","png","tiff")) {
      filename <- paste0(prefix,".",ext)
      if(ext=="pdf") cairo_pdf(filename,width=WIDTH_MM/25.4,height=height_mm/25.4,pointsize=BASE_PT)
      if(ext=="svg") svglite::svglite(filename,width=WIDTH_MM/25.4,height=height_mm/25.4,pointsize=BASE_PT)
      if(ext=="png") png(filename,width=WIDTH_MM,height=height_mm,units="mm",res=600,pointsize=BASE_PT,type="cairo")
      if(ext=="tiff") tiff(filename,width=WIDTH_MM,height=height_mm,units="mm",res=600,pointsize=BASE_PT,type="cairo",compression="lzw")
      tryCatch(draw(),finally=dev.off())
    }
    writeLines(c("width_mm=178; height_mm=55; font=Arial",
      "font_pt: panel=10; axis_title=7.5; ticks/legend=7; gene_labels/status=6.5",
      "placement=Figure4 full-width row beneath A/B/C; intended_scale=1",
      "margins_inches: bottom=0.50; left=0.46; top=0.16; right=0.08; ppa_title_line=1.9",
      "status_label_bottom_anchor_mm=3; status_label_vertical_alignment=bottom",
      "assembly_and_caption_page_fit=NOT_CHECKED",
      paste("font_paths",paste(font_paths,collapse=";"))),
      file.path(out,"Figure4D_layout.txt"))
    write_tsv(m,"Figure4D_variant_source.tsv")
    write_tsv(m[m$credible_set_member,],"Figure4D_credible_set_source.tsv")
    write_tsv(s,"Figure4D_sample_source.tsv")
    write_tsv(g,"Figure4D_annotation_source.tsv")
    unions <- do.call(rbind,lapply(ids,function(id) {
      e <- exon_union(g[g$gene_id==id,]); if(!nrow(e)) return(NULL)
      data.frame(gene_id=id,e)
    }))
    write_tsv(unions,"Figure4D_exon_union_source.tsv")
    writeLines(c("Figure 4D. Original PLASMA implementation output at DOG5_NEWGENE_03157 in frontal cortex.",
      sprintf("The frontal-cortex block includes %d dogs; orange denotes the candidate locus.",nrow(s)),
      sprintf("All %d modeled variants are shown. Blue denotes the %d members read directly from the upstream 95%% credible set; PPA values are unmodified.",nrow(m),sum(cs)),
      "The locus was selected for concentrated conditional posterior support, not as a typical block or based on GWAS/LINK status.",
      "The upstream model searches exactly one causal variant, omits the null configuration, and uses global overdispersion zero; these are not significance or causal-validation results.",
      "Gene tracks show exon unions across annotated transcripts, clipped to the observed marker range; unions are not single transcript models. Arrows indicate annotated strand where space permits; paired slashes indicate clipping. Reference tracks are packed for display. Nearby genes imply no regulatory link.",
      "No jitter, smoothing, LD colouring, unique lead selection, or PPA renormalization is applied. Credible-set membership is indicated by blue points without a duplicate position track. Nearby variants may overlap visually; in this block, the three members at chr30:38084279, 38084344 and 38084854 span 575 bp.",
      sprintf("The maximum PPA is %.6f and is shared by %d modeled variants; no unique lead is designated. The reported credible-set posterior sum is %.6f. The 95%% label is target cumulative model-conditional coverage, not the probability of each member.",
              max(ppa),sum(ppa==max(ppa)),sum(ppa[cs==1])),
      "Sample haplotype support is a sum of site-level allele counts, not independent RNA molecules; haplotype labels are individual-specific.",
      "Prior targeted diagnosis of this block identified allele-encoding dependence in total-expression residual calculation; a centered-residual single-factor comparison changed PPA and credible-set size. This figure preserves the original implementation output; variant priorities and mapping resolution are not finalized. Diagnosis is not rerun by this plotting script."),file.path(out,"Figure4D_caption.txt"))
    writeLines(c(paste0("input_dir=",INPUT),paste0("result_dir=",RESULT),paste0("gtf=",GTF),
      "displayed_result=original_PLASMA_output",
      "targeted_diagnosis=completed_for_this_block","encoding_dependence=identified",
      "diagnosis_status_source=prior_targeted_review_not_recomputed_by_plot",
      "diagnostic_artifacts=not_archived_with_this_figure",
      "final_result_version=not_finalized",
      "rendered_visual_QA=pending","alignment=single_shared_coordinate_system",
      "concurrent_tasks=1","threads_per_task=1",readLines(file.path(RESULT,"command.txt")),
      "Original stderr:",readLines(file.path(RESULT,"run.stderr.log"))),file.path(out,"Figure4D_provenance.txt"))
    capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
    log("[4/4] file write finished")
    writeLines("execution_status=success",file.path(out,"execution_completion.marker"))
    log("task completed / done; execution_status=success")
    log("please check separately output completeness, font, occlusion, figure interpretation and numerical diagnostics; this marker does not indicate that scientific or layout review has passed.")
  },error=function(e) {log(paste("execution_status=failed",conditionMessage(e))); stop(e)})
}
main()
