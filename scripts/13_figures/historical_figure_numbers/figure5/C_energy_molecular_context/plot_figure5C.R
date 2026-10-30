#!/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/Rscript


Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           NUMEXPR_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
ROOT <- "/path/to/project/analysis_20260817"
OUT_ROOT <- file.path(ROOT,"article/figure/figure5/C_energy_molecular_context")
BLOCK <- "no_gwas__cerebellum__DOG5_NEWGENE_01873__chr15"
MOL_INPUT <- file.path(ROOT,"5-dogs/10.finemaping/formal/8.2_no_gwas_prepare_expression_ase_plasma_inputs/20260821_181104/plasma_inputs",BLOCK)
MOL_RESULT <- file.path(ROOT,"5-dogs/10.finemaping/formal/8.3_no_gwas_run_expression_ase_plasma_joint/20260822_233245/plasma_results",BLOCK)
GWAS_DIR <- "/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.10_run_dog_behavior_susie_rss_finemapping/20260802_205352/main/primary/score_energy/DOG_GWAS_PRIMARY__score_energy__chr15_001"
GTF <- file.path(ROOT,"reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/GCF_000002285.3_CanFam3.1_genomic_chr.plus_five_only_new_gene.gtf")
LINK_INPUT <- file.path(ROOT,"article/figure/figure5/B_candidate_partner_support/preparation/20260929_231104/Figure5B_pair_sample_evidence.tsv")
CANDIDATE <- "DOG5_NEWGENE_01873"
PARTNER <- "gene-PAH"
FOCUS_ID <- "Affx-206394015"
FONT <- "Arial"
FOCUS_WIDTH_MM <- 87
LEFT_MARGIN_MM <- 13
RIGHT_MARGIN_MM <- 2
BLUE <- "#0072B2"
GREEN <- "#009E73"
GREY <- "#B5B5B5"


LINK_COLOURS <- c(enhancing="#478577", repressing="#BC220D")

read_tsv <- function(path) read.delim(path,quote="\"",comment.char="",encoding="UTF-8",
  stringsAsFactors=FALSE,check.names=FALSE)
write_tsv <- function(x,path) write.table(x,path,sep="\t",quote=TRUE,row.names=FALSE,
  na="NA",fileEncoding="UTF-8",qmethod="double")
log_stage <- function(x) {cat(format(Sys.time(),"[%Y-%m-%d %H:%M:%S ]"),x,"\n"); flush.console()}
attr_value <- function(x,key) {
  hit <- regmatches(x,regexec(paste0('(?:^|;)[[:space:]]*',key,' "([^"]+)"'),x,perl=TRUE))
  vapply(hit,function(a) if(length(a)>1) a[2] else NA_character_,character(1))
}

variant_key <- function(chrom,pos,a,b) paste(sub("^chr","",chrom),pos,
  pmin(toupper(a),toupper(b)),pmax(toupper(a),toupper(b)),sep=":")
read_chr15_exons <- function() {
  con <- file(GTF,"r"); on.exit(close(con))
  chunks <- list()
  repeat {
    lines <- readLines(con,n=50000,warn=FALSE)
    if(!length(lines)) break
    lines <- lines[grepl("^chr15\t[^\t]+\texon\t",lines)]
    if(!length(lines)) next
    d <- read.delim(text=paste(lines,collapse="\n"),header=FALSE,quote="",comment.char="")
    chunks[[length(chunks)+1]] <- data.frame(gene_id=attr_value(d$V9,"gene_id"),
      gene_name=attr_value(d$V9,"gene_name"),transcript_id=attr_value(d$V9,"transcript_id"),
      start=d$V4,end=d$V5,strand=d$V7)
  }
  unique(do.call(rbind,chunks))
}
spans <- function(e,key) {
  do.call(rbind,lapply(split(e,e[[key]]),function(z) data.frame(id=z[[key]][1],
    gene_id=z$gene_id[1],gene_name=z$gene_name[1],strand=z$strand[1],
    start=min(z$start),end=max(z$end))))
}
union_exons <- function(e) {
  do.call(rbind,lapply(split(e,e$gene_id),function(z) {
    z <- z[order(z$start,z$end),]; result <- list(); a <- z$start[1]; b <- z$end[1]
    if(nrow(z)>1) for(i in 2:nrow(z)) {
      if(z$start[i]<=b+1) b <- max(b,z$end[i]) else {
        result[[length(result)+1]] <- c(a,b); a <- z$start[i]; b <- z$end[i]
      }
    }
    result[[length(result)+1]] <- c(a,b); v <- do.call(rbind,result)
    data.frame(gene_id=z$gene_id[1],start=v[,1],end=v[,2],model="gene_exon_union")
  }))
}

pack_genes <- function(g,window,width_mm=178) {
  g <- g[order(g$start,g$end),]; g$row <- 0L; g$label_x <- NA_real_
  occupied_start <- occupied_end <- numeric(nrow(g))
  g$label <- sub("^DOG5_NEWGENE_","",g$gene_name)
  for(i in seq_len(nrow(g))) {
    center <- mean(c(max(g$start[i],window[1]),min(g$end[i],window[2])))

    weight <- if(g$id[i] %in% c(CANDIDATE,PARTNER,"gene-IGF1","gene-LOC111090065")) "bold" else "normal"
    label_mm <- systemfonts::string_width(g$label[i],family=FONT,weight=weight,
                                         size=6.5,res=72)/72*25.4
    pad <- (label_mm/2+0.8)/(width_mm-LEFT_MARGIN_MM-RIGHT_MARGIN_MM)*diff(window)
    center <- max(window[1]+pad,min(window[2]-pad,center))
    g$label_x[i] <- center
    occupied_start[i] <- min(max(g$start[i],window[1]),center-pad)
    occupied_end[i] <- max(min(g$end[i],window[2]),center+pad)+diff(window)*0.006
  }


  first_row <- which(g$id %in% c(CANDIDATE,PARTNER))
  g$row[first_row] <- 1L
  rows <- list(first_row)
  remaining <- which(!g$id %in% c(CANDIDATE,PARTNER))

  priority <- match(g$id[remaining],c("gene-IGF1","gene-LOC111090065"))
  priority[is.na(priority)] <- 3L
  remaining <- remaining[order(priority,occupied_start[remaining],occupied_end[remaining])]
  for(i in remaining) {
    chosen <- length(rows)+1L
    if(length(rows)>=2L) for(r in 2:length(rows)) {
      jj <- rows[[r]]
      collision <- occupied_start[i] <= occupied_end[jj] &
                   occupied_end[i] >= occupied_start[jj]
      if(!any(collision)) { chosen <- r; break }
    }
    if(chosen>length(rows)) rows[[chosen]] <- integer(0)
    rows[[chosen]] <- c(rows[[chosen]],i)
    g$row[i] <- chosen
  }
  g
}
canvas <- function(window,ylim) plot(NA,xlim=window/1e6,ylim=ylim,xaxs="i",yaxs="i",
  axes=FALSE,xlab="",ylab="",bty="n")
coord_axis <- function(window) {
  ticks <- pretty(window/1e6,n=if(par("fin")[1]*25.4<100) 3 else 5); ticks <- ticks[ticks>=window[1]/1e6 & ticks<=window[2]/1e6]
  segments(window[1]/1e6,par("usr")[3],window[2]/1e6,par("usr")[3],lwd=0.6,col="#444444")
  axis(1,at=ticks,labels=format(ticks,trim=TRUE),tcl=-0.2,lwd=0,lwd.ticks=0.5)
  mtext("Chr15 position (Mb)",1,line=1.7,cex=7.5/7)
}
prob_track <- function(d,value,cs,window,label,colour,focus,lead_id=NULL,
                       y_label, y_max=NULL, subtitle=NULL) {
  ymax <- if(is.null(y_max)) max(d[[value]],na.rm=TRUE)*1.22 else y_max
  canvas(window,c(0,ymax)); abline(h=0,col="#DDDDDD")

  if(value=="neg_log10_p") {
    abline(h=5,col="#AAAAAA",lty=2,lwd=0.6)

    text((window[2]-diff(window)*0.008)/1e6,5+ymax*0.14,
      labels="P = 1e-5",adj=c(1,0.5),cex=6.5/7,col="#666666")
  }
  abline(v=focus/1e6,col="#777777",lty=3,lwd=0.6)
  points(d$pos/1e6,d[[value]],pch=16,cex=0.45,col=GREY)
  ii <- which(d[[cs]])
  points(d$pos[ii]/1e6,d[[value]][ii],pch=16,cex=0.65,col=colour)
  ff <- which(d$is_focus)
  points(d$pos[ff]/1e6,d[[value]][ff],pch=21,cex=1.05,lwd=0.9,bg=NA,col="#222222")
  if(!is.null(lead_id)) {
    tt <- which(d$ld_variant_id==lead_id)
    points(d$pos[tt]/1e6,d[[value]][tt],pch=18,cex=0.9,col="#222222")
  }

  ticks <- if(is.null(y_max)) pretty(c(0,ymax),n=3) else seq(0,ymax,by=0.05)
  ticks <- ticks[ticks>=0 & ticks<=ymax]
  segments(window[1]/1e6,0,window[1]/1e6,ymax,col="#444444",lwd=0.6)
  axis(2,at=ticks,las=1,tcl=-0.2,lwd=0,lwd.ticks=0.5)
  mtext(y_label,2,line=2.6,las=0,cex=7.5/7)

  mtext(label,3,line=0.4,adj=0,cex=8/7,font=1)
  if(!is.null(subtitle)) mtext(subtitle,3,line=0.4,adj=1,cex=6.5/7,font=1,col="#555555")
}
draw_structure <- function(g,e,window,key,y,colour) {
  a <- max(g$start,window[1]); b <- min(g$end,window[2])
  segments(a/1e6,y,b/1e6,y,col=colour,lwd=0.65)
  ee <- e[e[[key]]==g$id & e$end>=window[1] & e$start<=window[2],,drop=FALSE]
  rect(pmax(ee$start,window[1])/1e6,y-0.11,pmin(ee$end,window[2])/1e6,y+0.11,
    col=colour,border=NA)

  if(b-a>diff(window)*0.015) {
    mid <- mean(c(a,b)); dx <- min((b-a)*0.25,diff(window)*0.015)
    if(g$strand=="-") dx <- -dx

    arrow_col <- if(g$gene_id==CANDIDATE) "white" else colour
    arrows((mid-dx)/1e6,y,(mid+dx)/1e6,y,length=0.045,col=arrow_col,lwd=0.8)
  }
  if(g$start<window[1]) text(window[1]/1e6,y,"/",cex=6.5/7,col=colour)
  if(g$end>window[2]) text(window[2]/1e6,y,"/",cex=6.5/7,col=colour)
}
draw_region <- function(gw,mol,genes,exons,window,focus,lead_id,n_samples,link_annotation,focused=FALSE) {
  nr <- max(genes$row)

  layout(matrix(1:4,ncol=1),heights=if(focused) c(20,20,20,40) else c(20,20,20,nr*5+16))

  par(family=FONT,ps=7,cex=1,cex.axis=1,omi=c(if(focused) 0 else 10,0,4,0)/25.4,mgp=c(1.7,0.45,0),tcl=-0.2)
  par(mai=c(1.5,LEFT_MARGIN_MM,5.5,RIGHT_MARGIN_MM)/25.4)
  prob_track(gw,"neg_log10_p","is_cs",window,"Energy GWAS",BLUE,focus,lead_id,
    y_label="-log10(P)")
  prob_track(gw,"pip","is_cs",window,"GWAS fine-mapping",BLUE,focus,lead_id,
    y_label="PIP",y_max=0.20,subtitle=paste0(sum(gw$is_cs)," CS variants"))


  f <- gw[gw$is_focus,,drop=FALSE]

  focus_label <- sprintf("%s\nPIP = %.3f",f$ld_variant_id,f$pip)
  label_width <- max(strwidth(strsplit(focus_label,"\n")[[1]],cex=6.5/7))
  label_x <- min((f$pos+diff(window)*0.015)/1e6,window[2]/1e6-label_width)
  text(label_x,0.185,focus_label,adj=c(0,1),cex=6.5/7,col="#444444")
  prob_track(mol,"ppa","is_cs",window,"Molecular fine-mapping",GREEN,focus,
    y_label="PPA",y_max=0.20,subtitle=paste0("01873 / CB; n = ",n_samples," dogs"))

  par(mai=c(9,LEFT_MARGIN_MM,0.5,RIGHT_MARGIN_MM)/25.4)
  canvas(window,c(0,nr+2.2))
  abline(v=focus/1e6,col="#777777",lty=3,lwd=0.6)
  candidate_gene <- genes[genes$id==CANDIDATE,]
  partner_gene <- genes[genes$id==PARTNER,]
  link_x <- c(mean(c(candidate_gene$start,candidate_gene$end)),
              mean(c(partner_gene$start,partner_gene$end)))/1e6
  lines(c(link_x[1],link_x[1],link_x[2],link_x[2]),
        c(nr+0.8,nr+1.05,nr+1.05,nr+0.8),col=link_annotation$colour,lwd=0.6)
  text(mean(link_x),nr+1.55,paste0("Predicted ",tolower(link_annotation$label)," LINK"),
       cex=6.5/7,col=link_annotation$colour)
  for(i in seq_len(nrow(genes))) {
    z <- genes[i,]; y <- nr-z$row+0.5
    key_gene <- z$id %in% c(CANDIDATE,PARTNER,"gene-IGF1","gene-LOC111090065")
    label_font <- if(key_gene) 2 else 1
    col <- if(z$id==CANDIDATE) GREEN else if(key_gene) "#333333" else "#808080"
    draw_structure(z,exons,window,"gene_id",y,col)

    half_width <- strwidth(z$label,units="user",cex=6.5/7,font=label_font)/2
    inset <- diff(window)/1e6*0.004
    label_x <- max(window[1]/1e6+half_width+inset,
      min(window[2]/1e6-half_width-inset,z$label_x/1e6))
    text(label_x,y+0.34,z$label,cex=6.5/7,col=col,font=label_font)
  }
  coord_axis(window)


  if (focused) grid::grid.text("C",x=grid::unit(2,"mm"),
    y=grid::unit(1,"npc")-grid::unit(1.5,"mm"),just=c("left","top"),
    gp=grid::gpar(fontfamily=FONT,fontsize=10,fontface="bold",col="#252525"))
  if (!focused) {
  mtext("Blue: GWAS CS   Green: molecular CS   Grey: other modeled variants",1,outer=TRUE,line=0.3,cex=6.5/7)
  mtext("Circle: selected shared variant   Diamond: GWAS lead",1,outer=TRUE,line=1.3,cex=6.5/7)
  mtext("Original molecular output; model review pending (see caption).",1,outer=TRUE,line=2.3,cex=6.5/7)
  }
}
draw_local <- function(e,tx,window) {
  par(family=FONT,ps=7,cex=1,cex.axis=1,mar=c(3.5,1.5,5,1.5),mgp=c(2,0.4,0))
  canvas(window,c(0,nrow(tx)+0.8))
  for(i in seq_len(nrow(tx))) {
    z <- tx[i,]; y <- nrow(tx)-i+0.5; col <- if(z$gene_id==CANDIDATE) GREEN else "#444444"
    draw_structure(z,e,window,"transcript_id",y,col)
    text(window[1]/1e6,y+0.35,paste(z$gene_name,z$id,paste0("(",z$strand,")"),sep=" | "),adj=0,cex=1,col=col)
  }
  coord_axis(window)
  mtext("01873 local transcript structure",3,line=3,font=2,adj=0,cex=8/7)
  mtext("Annotated exons and strands; no RNA coverage or allelic phase inference",3,line=1.6,adj=0,cex=1)
  mtext("Slash: transcript continues outside this local window",3,line=0.4,adj=0,cex=6.5/7)
}
export_plot <- function(out,stem,height,draw,width=178) {
  for(fmt in c("pdf","svg","png","tiff")) {
    f <- file.path(out,paste0(stem,".",fmt))
    if(fmt=="pdf") grDevices::cairo_pdf(f,width=width/25.4,height=height/25.4,pointsize=7,family=FONT)
    if(fmt=="svg") svglite::svglite(f,width=width/25.4,height=height/25.4,pointsize=7,bg="white")
    if(fmt=="png") grDevices::png(f,width=width,height=height,units="mm",res=600,type="cairo",pointsize=7,bg="white")
    if(fmt=="tiff") grDevices::tiff(f,width=width,height=height,units="mm",res=600,type="cairo",pointsize=7,compression="lzw",bg="white")
    tryCatch(draw(),finally=grDevices::dev.off())
  }
  writeLines(c(paste0("width_mm=",width,"; height_mm=",height,"; intended_scale=1"),
    "font=Arial; panel_pt=10; track_heading_pt=8; axis_title_pt=7.5; tick_pt=7; dense_label_pt=6.5",
    "panel_inset_mm=left:2,top:1.5; molecular_sample_note_pt=6.5; background_gene_colour=#808080",
    "raster_ppi=600; TIFF=LZW; PDF_font_embedding=NOT_CHECKED; visual_layout=NOT_CHECKED",
    paste0("assembled_Figure5_and_caption=NOT_CHECKED; intended_width_mm=",width,"; no reduction")),
    file.path(out,paste0(stem,"_layout.txt")))
  log_stage(paste(stem,"width mm=",width,"height mm=",height))
}

draw_shared_legend <- function() {
  grid::grid.newpage()
  labels <- c("GWAS CS","Molecular CS","Other modeled variant","Selected variant","GWAS lead")
  colours <- c(BLUE,GREEN,GREY,"#222222","#222222")
  shapes <- c(16,16,16,1,18)
  widths <- vapply(labels,function(label) grid::convertWidth(grid::grobWidth(
    grid::textGrob(label,gp=grid::gpar(fontfamily=FONT,fontsize=7))),"mm",valueOnly=TRUE),numeric(1))

  item_widths <- widths+9
  left <- (178-sum(item_widths))/2
  for(i in seq_along(labels)) {
    grid::grid.points(x=grid::unit(left+1,"mm"),y=grid::unit(4,"mm"),pch=shapes[i],
      size=grid::unit(2,"mm"),gp=grid::gpar(col=colours[i],lwd=0.8))
    grid::grid.text(labels[i],x=grid::unit(left+4,"mm"),y=grid::unit(4,"mm"),just="left",
      gp=grid::gpar(fontfamily=FONT,fontsize=7,col="#252525"))
    left <- left+item_widths[i]
  }
}
main <- function() {
  out <- file.path(OUT_ROOT,format(Sys.time(),"%Y%m%d_%H%M%S"))
  if(!dir.create(out,showWarnings=FALSE)) stop("cannot create output directory: ",out)
  con <- file(file.path(out,"plot_figure5C.log"),"wt")
  sink(con); sink(con,type="message")
  on.exit({sink(type="message");sink();close(con)},add=TRUE)
  tryCatch({
    log_stage("[1/4] read original fine-mapping inputs; concurrent_tasks=1 threads_per_task=1")
    font_paths <- systemfonts::match_fonts(c(FONT,FONT),weight=c("normal","bold"))$path
    expected_fonts <- c("/path/to/fonts/Arial/Arial.TTF",
                        "/path/to/fonts/Arial/Arialbd.TTF")
    if(!all(file.exists(expected_fonts)) ||
       !identical(normalizePath(font_paths),normalizePath(expected_fonts)))
      stop("Arial regular or bold font is unavailable; do not silently substitute fonts.")
    log_stage(paste("Arial font",paste(font_paths,collapse=";")))
    gw <- read_tsv(file.path(GWAS_DIR,"variant_pip.tsv.gz"))

    gw$neg_log10_p <- -log10(gw$p_lrt)
    info <- read_tsv(file.path(GWAS_DIR,"locus_summary.tsv"))
    mol <- read_tsv(file.path(MOL_INPUT,"markers.tsv"))
    samples <- read_tsv(file.path(MOL_INPUT,"samples.tsv"))
    mol$ppa <- scan(file.path(MOL_RESULT,"ppas.txt"),quiet=TRUE)
    mol$is_cs <- scan(file.path(MOL_RESULT,"cset.txt"),quiet=TRUE)==1
    mol$input_marker_index <- seq_len(nrow(mol))
    gw$is_cs <- gw$in_credible_set=="yes"
    gw$key <- variant_key(gw$chrom,gw$pos,gw$ld_a1,gw$ld_a2)
    mol$key <- variant_key(mol$chrom,mol$pos,mol$ref,mol$alt)
    shared <- intersect(gw$key[gw$is_cs],mol$key[mol$is_cs])
    gw$shared_cs <- gw$is_cs & gw$key %in% shared
    mol$shared_cs <- mol$is_cs & mol$key %in% shared
    gw$is_focus <- gw$ld_variant_id==FOCUS_ID
    focus_key <- gw$key[gw$is_focus]; focus <- gw$pos[gw$is_focus]
    mol$is_focus <- mol$key %in% focus_key
    joined <- merge(gw,mol,by="key",suffixes=c("_gwas","_molecular"),sort=FALSE)


    coverage <- mol[mol$is_cs,,drop=FALSE]
    coverage$gwas_position_present <- paste(sub("^chr","",coverage$chrom),coverage$pos) %in%
      paste(sub("^chr","",gw$chrom),gw$pos)
    coverage$gwas_exact_match <- coverage$key %in% gw$key
    gwas_lookup <- gw[,c("key","ld_variant_id","p_lrt","pip","is_cs")]
    names(gwas_lookup) <- c("key","gwas_variant_id","gwas_p_lrt","gwas_pip","gwas_is_cs")
    coverage <- merge(coverage,gwas_lookup,by="key",all.x=TRUE,sort=FALSE)
    coverage <- coverage[order(coverage$input_marker_index),]
    links <- read_tsv(LINK_INPUT)
    links <- links[links$new_gene_id==CANDIDATE & links$partner_gene_id==PARTNER,,drop=FALSE]

    link_scope <- links[links$tissue=="cerebellum" & links$pair_state=="detected",,drop=FALSE]
    mechanisms <- unique(na.omit(link_scope$reported_mechanisms))
    if(length(mechanisms)!=1L || !mechanisms %in% names(LINK_COLOURS))
      stop("cerebellum LINK mechanism missing or is not unique, cannot assign a single enhancing label to the link/repressing color")
    mechanism <- mechanisms[1]
    dogs <- paste(sort(unique(sub("_.*", "", link_scope$animal_id))),collapse=", ")
    link_annotation <- list(colour=unname(LINK_COLOURS[mechanism]),
      label=if(mechanism=="enhancing") "Enhancing" else "Repressing")
    log_stage(paste("GWAS variant",nrow(gw),"molecular variant",nrow(mol),"molecular sample",nrow(samples),"complete CS intersection",length(shared)))
    log_stage("[2/4] read original GTF; region use exon union, local retain transcript_id")
    e <- read_chr15_exons(); g <- spans(e,"gene_id")
    target <- g[g$id==CANDIDATE,]; partner <- g[g$id==PARTNER,]
    window <- range(c(gw$pos,mol$pos,target$start,target$end,partner$start,partner$end))+c(-10000,10000)
    genes <- g[g$start<=window[2] & g$end>=window[1],]
    region_e <- e[e$gene_id %in% genes$id,]; unions <- union_exons(region_e)
    genes <- pack_genes(genes,window)

    focused_window <- range(c(gw$pos[gw$is_cs],mol$pos,target$start,target$end,
      partner$start,partner$end,gw$pos[gw$ld_variant_id==info$gwas_lead_variant_id]))+c(-10000,10000)

    focused_genes_all <- g[g$start<=focused_window[2] & g$end>=focused_window[1],]
    focus_gene_ids <- c(CANDIDATE,PARTNER,"gene-IGF1","gene-LOC111090065")
    focused_genes_all$displayed_in_focus <- focused_genes_all$id %in% focus_gene_ids
    focused_genes <- pack_genes(focused_genes_all[focused_genes_all$displayed_in_focus,],
                               focused_window,width_mm=FOCUS_WIDTH_MM)
    log_stage(paste("focus window all gene",nrow(focused_genes_all),
                    "main text display gene",nrow(focused_genes),"track count",max(focused_genes$row)))
    local_window <- c(target$start-2000,target$end+2000)
    tx <- spans(region_e,"transcript_id")
    tx <- tx[tx$start<=local_window[2] & tx$end>=local_window[1],]
    tx <- tx[order(tx$gene_id==CANDIDATE,tx$gene_name,tx$id),]
    local_e <- region_e[region_e$transcript_id %in% tx$id,]

    selected_tx <- tx[tx$gene_id==CANDIDATE | tx$id=="rna-XR_002614240.1",]
    log_stage("[3/4] export regional integration plots and separate transcript-structure detail plots")
    export_plot(out,"Figure5C_energy_molecular_context",103+max(genes$row)*7,
      function() draw_region(gw,mol,genes,unions,window,focus,info$gwas_lead_variant_id,nrow(samples),link_annotation))
    export_plot(out,"Figure5C_energy_molecular_context_focus",max(100,80+max(focused_genes$row)*5),
      function() draw_region(gw,mol,focused_genes,unions,focused_window,focus,info$gwas_lead_variant_id,nrow(samples),link_annotation,focused=TRUE),width=FOCUS_WIDTH_MM)
    export_plot(out,"Figure5C_01873_transcript_detail",35+nrow(tx)*10,
      function() draw_local(local_e,tx,local_window))
    export_plot(out,"Figure5C_01873_transcript_focus",35+nrow(selected_tx)*10,
      function() draw_local(local_e,selected_tx,local_window))
    export_plot(out,"Figure5CD_shared_legend",8,draw_shared_legend,width=178)
    writeLines(c("Standalone C/D shared legend: intended width 178 mm, height 8 mm, Arial 7 pt.",
      "Selected variant: C, the matched variant Affx-206394015; D, candidate-overlapping Affx-206206334.",
      "Blue/green denote GWAS/molecular CS membership; grey denotes other modeled variants; diamonds denote GWAS leads.",
      "When using this shared legend below C/D, omit the duplicate D legend during an explicitly coordinated assembly; this script does not modify D.",
      "Adding this strip increases assembled height by 8 mm plus any gap; verify the full figure and caption together."),
      file.path(out,"Figure5CD_shared_legend_notes.txt"))
    log_stage("[4/4] save source table, source, figure legend with run status")
    tables <- list(gwas_variants=gw,molecular_variants=mol,exact_overlap=joined,
      molecular_CS_GWAS_coverage=coverage,
      shared_cs=joined[joined$is_cs_gwas & joined$is_cs_molecular,,drop=FALSE],
      molecular_samples=samples,region_genes=genes,region_exon_union=unions,
      focused_genes=focused_genes,focused_genes_all=focused_genes_all,local_transcripts=tx,local_transcripts_selected=selected_tx,
      local_exons=local_e,LINK_context=links)
    for(n in names(tables)) write_tsv(tables[[n]],file.path(out,paste0("Figure5C_",n,".tsv")))
    write_tsv(data.frame(view=c("region","focused_region","local_transcripts"),chrom="chr15",
      start=c(window[1],focused_window[1],local_window[1]),end=c(window[2],focused_window[2],local_window[2]),coordinate_system="1-based inclusive"),file.path(out,"Figure5C_windows.tsv"))
    inputs <- c(gwas=file.path(GWAS_DIR,"variant_pip.tsv.gz"),gwas_summary=file.path(GWAS_DIR,"locus_summary.tsv"),
      markers=file.path(MOL_INPUT,"markers.tsv"),samples=file.path(MOL_INPUT,"samples.tsv"),
      ppa=file.path(MOL_RESULT,"ppas.txt"),cs=file.path(MOL_RESULT,"cset.txt"),gtf=GTF,links=LINK_INPUT)
    write_tsv(data.frame(role=names(inputs),absolute_path=unname(inputs)),file.path(out,"Figure5C_input_paths.tsv"))

    writeLines(c("C. Genetic and molecular fine-mapping at the Energy-associated region (CanFam3.1).",
      paste0("The mechanism bracket uses cerebellum LINK records from ",dogs,"; reported classification: ",mechanism,"."),
      paste0("Tracks show GWAS -log10(P), GWAS PIP and molecular PPA for DOG5_NEWGENE_01873 in cerebellum (",nrow(samples),
        " dogs). The focused view displays ",sum(gw$pos>=focused_window[1] & gw$pos<=focused_window[2])," of ",nrow(gw),
        " GWAS-modeled variants and ",sum(mol$pos>=focused_window[1] & mol$pos<=focused_window[2])," molecular variants."),
      paste0("Blue and green indicate GWAS and molecular credible-set (CS) membership (n=",sum(gw$is_cs)," and ",sum(mol$is_cs),
        "); grey indicates other modeled variants. Diamonds mark the GWAS lead; circles and vertical dotted lines mark ",FOCUS_ID,
        ". The horizontal dashed line denotes P=1e-5."),
      paste0("Molecular CS coverage in the GWAS input: ",sum(coverage$gwas_exact_match),"/",sum(mol$is_cs),
        " exact coordinate-and-allele-set matches; shared CS membership: ",length(shared),
        ". Missing GWAS records do not indicate lack of association."),
      "The focused panel displays selected genes only: DOG5_NEWGENE_01873, PAH, IGF1 and LOC111090065. Other genes are omitted from this panel for readability, not because they lack biological relevance. The full regional view retains all overlapping annotated genes; focused_genes_all records every gene in the focused window and its display status.",
      "Five-digit gene labels use the DOG5_NEWGENE_ prefix; full IDs are retained in the source tables. The molecular sample count is the model input size, not the number of dogs supporting LINK in panel A or B.",
      "Gene tracks show exon unions; 01873 and PAH occupy the first gene row, with other genes packed below for readability. Vertical row placement has no biological meaning. The undirected bracket denotes the predicted 01873-PAH LINK mechanism from detected cerebellum records (see B); the label explicitly identifies enhancing or repressing as a predicted LINK classification, not an experimentally validated effect. Mechanism colours follow the author Allelome.LINK BEDPE implementation: enhancing #478577 and repressing #BC220D (https://github.com/AndergassenLab/Allelome.LINK/blob/main/00_src/Allelome.LINK.R; doi:10.1093/nar/gkaf912). Other-tissue calls are not pooled into this mechanism label. This classification is not experimental proof of enhancement or repression, nor does it establish which gene regulates the other. LINK attribution and endpoint phase comparability remain unresolved. Molecular PPA retains the original PLASMA output and awaits model review; this plotting script does not assess numerical diagnostics or implementation stability. CS overlap does not establish colocalization or a causal regulatory mechanism."
    ),file.path(out,"Figure5C_caption.txt"))
    writeLines(c("Figure 5C: detailed methods, provenance and interpretation notes.",
      paste("GWAS modeled variants:",nrow(gw),"; molecular modeled variants:",nrow(mol),"; cerebellum dogs:",nrow(samples)),
      paste("Focused window: GWAS variants shown",sum(gw$pos>=focused_window[1] & gw$pos<=focused_window[2]),"of",nrow(gw),
        "; molecular variants shown",sum(mol$pos>=focused_window[1] & mol$pos<=focused_window[2]),"of",nrow(mol),
        "; exact variant keys shared between models",length(intersect(gw$key,mol$key))),
      "All positions use the canine CanFam3.1 assembly. The association track shows -log10(p_lrt) for variants included in the GWAS fine-mapping model, not all variants tested in the original regional GWAS. Blue indicates GWAS CS membership, not LD or a P-value threshold.",
      "The dashed horizontal line at -log10(P)=5 marks the single-variant P=1e-5 reference. It does not filter plotted variants, change CS membership or replace the regional-lead inclusion criterion.",
      "The focused gene track shows only DOG5_NEWGENE_01873, PAH, IGF1 and LOC111090065. All regional genes remain in the full view and source tables; focused_genes_all retains all focused-window genes with displayed_in_focus. Vertical gene placement prioritizes the candidate-partner pair and prevents label overlap. Gene display selection does not filter variants or change the coordinate window.",
      "The molecular expression object is DOG5_NEWGENE_01873, not PAH. Probabilities and CS flags are read unchanged in original marker order. Missing model records are not assigned zero.",
      paste("GWAS CS members:",sum(gw$is_cs),"; molecular CS members:",sum(mol$is_cs),"; exact shared CS members:",length(shared)),
      paste("Molecular CS coverage counts: total",sum(mol$is_cs),"; represented GWAS coordinates",sum(coverage$gwas_position_present),
        "; exact allele-set matches",sum(coverage$gwas_exact_match),
        "; members of both CSs",sum(coverage$gwas_is_cs,na.rm=TRUE)),
      "Absence from the GWAS fine-mapping input is not evidence against a molecular CS variant. The reason for absence is not established here; probabilities are not renormalized over shared variants. Some molecular CS points overlap at this display scale.",
      "CS membership is encoded by blue/green points; separate CS position tracks are omitted. Open circles highlight only the selected integration variant in the association and posterior tracks. All shared members are retained in the shared_cs source table. Allele-set matching does not harmonize effect direction or perform formal colocalization.",
      "The regional window includes all modeled variants and PAH. Gene structures are exon unions. Local detail retains individual transcript IDs, with clipped ends marked; it contains annotation only, not RNA coverage.",
      "The focused regional view includes all GWAS CS members, all molecular variants, the GWAS lead and PAH, with 10 kb padding; other GWAS variants outside this window are visible in the full view and retained in the source table. Probability scales are unchanged between views.",
      "The compact local view displays the candidate transcript and its annotation reference rna-XR_002614240.1; the full local detail includes all overlapping transcripts. White arrows on the candidate exon indicate annotated strand only.",
      "The undirected LINK connects gene spans, not a tested SNP-to-gene mechanism; sample support is shown in panel B. GTF coordinates are 1-based inclusive; exported LINK endpoint coordinates retain their original convention.",
      "In the inspected LINK input, all four detected records carry nonoverlap_multigene_review; all three candidate information sites in each record carry multigene attribution flags, and endpoint phase comparability is not_assessed. These flags do not by themselves establish read multimapping or invalidate the reported links.",
      "Molecular PPA uses the original PLASMA implementation and awaits model review. This script does not read diagnostic or sensitivity-test outputs and cannot establish block-specific numerical warnings or the effect of implementation dependence. These require separate upstream verification; no corrected molecular result is inferred here.",
      "Multigene attribution and phase comparability remain upstream review items. ATAC, motifs, LD colouring and RNA coverage are not inferred by this script."
    ),file.path(out,"Figure5C_methods_notes.txt"))
    writeLines(c("Layout reference: Schwartzentruber et al., Nature Genetics (2021), doi:10.1038/s41588-020-00776-w; user-supplied regional figure.",
      "Author code locator: https://github.com/Jeremy37/AD_finemap/blob/master/src/finemapping_plots.R ; exact source body and final figure mapping could not be retrieved in this session. No code copied.",
      "Verified official track example and configuration: https://pygenometracks.readthedocs.io/en/latest/content/examples.html#basic-examples (genes, coordinate axis and vlines).",
      "Biological reference: Hasenbein et al., NAR (2025), doi:10.1093/nar/gkaf912; https://github.com/AndergassenLab/Allelome.LINK ; not treated as a verified figure-reproduction script.",
      "Implementation: independent R 4.5.1 base graphics; local help for read.table, scan, plot.default, segments, arrows, layout, axis, png and cairo_pdf reviewed."
    ),file.path(out,"Figure5C_references.txt"))
    capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
    writeLines(c("execution_status=success",paste0("output_dir=",out)),file.path(out,"execution_completion.marker"))
    log_stage("task completed / done; execution_status=success")
    log_stage("execution completed does not imply quality audit passed; please check separately layout, CS matching, coordinates with biological interpretation.")
  },error=function(e) {log_stage(paste("execution_status=failed;",conditionMessage(e)));stop(e)})
}
main()
