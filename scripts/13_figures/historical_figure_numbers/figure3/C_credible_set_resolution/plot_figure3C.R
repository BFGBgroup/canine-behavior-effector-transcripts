#!/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/Rscript


Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
           NUMEXPR_NUM_THREADS="1",VECLIB_MAXIMUM_THREADS="1")
ROOT <- "/path/to/repository/scripts/13_figures/historical_figure_numbers/figure3/C_credible_set_resolution"
WORK <- "/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work"
MANIFEST <- file.path(WORK,"1.10_run_dog_behavior_susie_rss_finemapping/20260802_205352/locus_run_manifest.tsv")
REGIONS <- file.path(WORK,"1.7_define_dog_behavior_ld_loci/20260802_171634/all_locus_regions.tsv")
METADATA <- file.path(WORK,"1.2_build_pnas_dog_behavior_locus_manifest/20260801_202130/trait_metadata.tsv")
PRODUCER <- file.path(WORK,"scripts/1.10_run_dog_behavior_susie_rss_finemapping/1.10_run_dog_behavior_susie_rss_finemapping.R")
REGION_PRODUCER <- file.path(WORK,"scripts/1.7_define_dog_behavior_ld_loci/1.7_define_dog_behavior_ld_loci.py")
META_PRODUCER <- file.path(WORK,"scripts/1.2_build_pnas_dog_behavior_locus_manifest/1.2_build_pnas_dog_behavior_locus_manifest.py")
OUT <- file.path(ROOT,paste0(format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid()))
if(!dir.create(OUT,recursive=TRUE,showWarnings=FALSE)) stop("cannot create new output directory")
con <- file(file.path(OUT,"plot_figure3C.log"),"wt")
sink(con); sink(con,type="message")
log_stage <- function(...) {cat(format(Sys.time(),"[%Y-%m-%d %H:%M:%S] "),...,"\n",sep="");flush(con)}
read_tsv <- function(path) read.delim(path,quote="",comment.char="",check.names=FALSE)
write_tsv <- function(x,name) write.table(x,file.path(OUT,name),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
main <- function() {

  .libPaths(c("/path/to/software/R/figure3-library",.libPaths()))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggbeeswarm))

  font_paths <- systemfonts::match_fonts(c("Arial","Arial"),weight=c("normal","bold"))$path
  expected_fonts <- c("/path/to/fonts/Arial/Arial.TTF",
                      "/path/to/fonts/Arial/Arialbd.TTF")
  if (!all(file.exists(expected_fonts)) ||
      !identical(normalizePath(font_paths),normalizePath(expected_fonts)))
    stop("Arial regular or bold font was not resolved correctly, please check font installation")
  log_stage("[1/3] read original credible sets; concurrent_tasks=1; threads_per_task=1")
  log_stage("MANIFEST=",MANIFEST,"\nREGIONS=",REGIONS,"\nMETADATA=",METADATA,"\nOUT=",OUT)
  m <- read_tsv(MANIFEST)
  m <- m[m$branch=="main" & m$analysis_scope %in% c("primary","exploratory"),]
  regions <- read_tsv(REGIONS)
  meta <- read_tsv(METADATA)
  meta <- meta[order(meta$trait),c("trait","trait_label")]

  if(anyNA(m$status) || any(m$status!="PASS")) stop("Some main tasks did not succeed; these must not be treated as having no credible set")
  idx <- match(m$locus_id,regions$locus_id)
  m$region_lead_p_lrt <- regions$lead_p_lrt[idx]
  m$trait_label <- meta$trait_label[match(m$trait,meta$trait)]
  rows <- vector("list",nrow(m))
  m$n_cs <- integer(nrow(m))
  for(i in seq_len(nrow(m))) {
    cs <- read_tsv(m$credible_sets[i])
    m$n_cs[i] <- nrow(cs)
    cs$source_file <- rep(m$credible_sets[i],nrow(cs))
    cs$region_lead_p_lrt <- rep(m$region_lead_p_lrt[i],nrow(cs))
    rows[[i]] <- cs
  }
  d <- do.call(rbind,rows)

  if(any(!is.finite(d$cs_size) | d$cs_size<1 | d$cs_size!=floor(d$cs_size)))
    stop("credible-set member counts contain non-positive integers or missing values")
  tiers <- c("Genome-wide","Suggestive","Exploratory")
  d$tier <- factor(ifelse(d$region_lead_p_lrt<2.97e-6,tiers[1],
                   ifelse(d$region_lead_p_lrt<1e-5,tiers[2],tiers[3])),levels=tiers)
  d$trait_label <- meta$trait_label[match(d$trait,meta$trait)]
  d <- d[order(d$trait,d$cs_size,d$tier,d$locus_id,d$cs_id),]
  d$y <- nrow(meta)+1-match(d$trait,meta$trait)

  focus <- d[d$locus_id=="DOG_GWAS_PRIMARY__score_energy__chr15_001" &
             d$branch=="main" & d$cs_id=="L1" & d$cs_size==30,]
  if(nrow(focus)!=1L) stop("not found unique of Figure 3D corresponding Energy credible set")
  focus$label_x <- focus$cs_size*1.33


  summary <- do.call(rbind,lapply(seq_len(nrow(meta)),function(i) {
    z <- d[d$trait==meta$trait[i],]; r <- m[m$trait==meta$trait[i],]
    data.frame(trait=meta$trait[i],trait_label=meta$trait_label[i],
      n_regions=nrow(r),regions_with_cs=sum(r$n_cs>0),regions_without_cs=sum(r$n_cs==0),
      n_credible_sets=nrow(z),median_cs_size=if(nrow(z)) median(z$cs_size) else NA_real_,
      min_cs_size=if(nrow(z)) min(z$cs_size) else NA_real_,
      max_cs_size=if(nrow(z)) max(z$cs_size) else NA_real_)
  }))
  log_stage("[2/3] region=",nrow(m),"; credible set=",nrow(d),"; regions without retained credible sets=",sum(m$n_cs==0))
  write_tsv(d,"Figure3C_credible_sets.tsv")
  write_tsv(summary,"Figure3C_trait_summary.tsv")
  write_tsv(m,"Figure3C_region_status.tsv")
  write_tsv(m[m$n_cs==0,],"Figure3C_regions_without_cs.tsv")
  write_tsv(data.frame(role=c("run_manifest","regions","metadata",rep("credible_sets",nrow(m))),
    path=c(MANIFEST,REGIONS,METADATA,m$credible_sets),
    producer=c(PRODUCER,REGION_PRODUCER,META_PRODUCER,rep(PRODUCER,nrow(m)))),"input_manifest.tsv")

  colors <- setNames(c("#0F4D92","#3775BA","#CFCECE"),tiers)
  outlines <- setNames(c("#0F4D92","#3775BA","#767676"),tiers)
  xmax <- max(d$cs_size)
  xupper <- max(55,ceiling((xmax+1)/5)*5)
  xb <- c(1,2,5,10,20,50)
  p <- ggplot(d,aes(cs_size,y,colour=tier,fill=tier,group=trait))+
    geom_hline(yintercept=seq_len(nrow(meta)),colour="#EEEEEE",linewidth=.22)+

    ggbeeswarm::geom_beeswarm(orientation="y",method="swarm",priority="ascending",
      preserve.data.axis=TRUE,dodge.width=NULL,corral="none",cex=1.2,
      shape=21,size=1.3,stroke=.18,alpha=1)+

    geom_point(data=focus,aes(x=cs_size,y=y),inherit.aes=FALSE,
      shape=21,fill=NA,colour="#222222",size=2.4,stroke=.35)+
    geom_text(data=focus,aes(x=label_x,y=y,label="D"),inherit.aes=FALSE,
      hjust=0,size=6.5,size.unit="pt",family="Arial",colour="#222222")+
    scale_fill_manual(values=colors,drop=FALSE)+
    scale_colour_manual(values=outlines,drop=FALSE)+

    scale_x_log10(breaks=xb,labels=as.character(xb),minor_breaks=NULL,
      limits=c(1,xupper),expand=expansion(mult=c(.035,.035)))+
    scale_y_continuous(breaks=seq_len(nrow(meta)),labels=rev(meta$trait_label),
      limits=c(.45,nrow(meta)+.55),expand=expansion(mult=0))+
    labs(x="Variants per CS (log scale)",y=NULL,title="Credible-set resolution",tag="C")+
    theme_classic(base_size=7,base_family="Arial")+
    theme(axis.text=element_text(size=6.5,colour="#222222"),

      axis.text.y=element_blank(),
      axis.title.x=element_text(size=7.5,margin=margin(t=5)),
      axis.line.y=element_blank(),axis.ticks.y=element_blank(),
      axis.line.x=element_line(linewidth=.3),axis.ticks.x=element_line(linewidth=.3),
      legend.position="none",plot.title=element_text(size=8,face="bold",margin=margin(b=5)),
      plot.tag=element_text(size=10,face="bold"),plot.tag.position="topleft",
      plot.margin=margin(5,5,5,5))
  log_stage("[3/3] export plot and description")
  prefix <- file.path(OUT,"Figure3C_credible_set_resolution")
  ggsave(paste0(prefix,".pdf"),p,device=grDevices::cairo_pdf,width=50,height=82,units="mm",bg="white")
  ggsave(paste0(prefix,".png"),p,device=ragg::agg_png,width=50,height=82,units="mm",dpi=600,bg="white")
  ggsave(paste0(prefix,".svg"),p,device=svglite::svglite,width=50,height=82,units="mm",bg="white")
  ggsave(paste0(prefix,".tiff"),p,device=ragg::agg_tiff,width=50,height=82,units="mm",dpi=600,bg="white",compression="lzw")
  saveRDS(p,paste0(prefix,".rds"))
  writeLines(c("panel_width_mm=50", "panel_height_mm=82",
    "font=Arial; panel_tag=10; title=8; axis_title=7.5; axis_text/D_label=6.5 pt",
    "trait_labels=shared_with_panel_A; order=ascending_trait_id",
    "y_limits=0.45 to n_traits+0.55; plot_margins=5 pt; x_axis_title=single_line",
    "ABC planned widths=75/45/50 mm; gaps=4 mm",
    "A/B/C title and axis typography/margins matched; rendered alignment=NOT_CHECKED",
    "beeswarm offsets are visual only; inspect final point separation and D ring after rendering",
    "PDF/SVG vector; PNG/TIFF 600 dpi; TIFF LZW; white background",
    paste("font_paths",paste(font_paths,collapse=";"))),
    file.path(OUT,"Figure3C.layout_parameters.txt"))
  writeLines(c("Figure 3C. Credible-set resolution by behavioral trait.",
    sprintf("Each point represents one retained 95%% credible set in the main branch: %d sets from %d of %d trait-region units.",nrow(d),sum(m$n_cs>0),nrow(m)),
    "The upstream susieR run requested nominal 95% coverage and used min_abs_corr=0.5 when identifying credible sets; this panel applies no additional purity filter.",
    "CS denotes credible set. Trait names are shared with panel A and omitted from this assembly panel. The horizontal axis shows the actual number of variants per credible set on a log10 scale, with ticks at 1, 2, 5, 10, 20 and 50. Equal distances represent equal ratios, not equal count differences. Rows follow panels A and B. ggbeeswarm arranges all colors within each trait together using deterministic vertical offsets while preserving the data axis. Vertical offsets have no quantitative meaning; exported y values are trait baselines, not rendered offsets. Final-size rendering still requires visual inspection.",
    "All points are equally sized circles. Dark blue: region lead LRT P < 2.97e-6; medium blue: 2.97e-6 <= P < 1e-5; light gray with gray outlines: 1e-5 <= P < 5e-5. Use the shared legend exported with panel A. Categories are based on region lead P, not member PIP.",
    "The black ring and D label identify only the Energy chr15_001 main-branch L1 credible set containing 30 variants, which is the regional example in panel D. This ring does not denote the specific integration-linked SNP highlighted within panel D.",
    sprintf("%d regions have no retained credible set; these are reported separately and not plotted as zeros. The plotted distribution is conditional on retaining a credible set.",sum(m$n_cs==0)),
    "All locus-CS records are retained without collapsing across traits or genomic overlaps. No PIP filter, median-only reduction, density estimation or between-trait significance test is applied. Credible sets are not biological replicates or experimentally confirmed causal signals.",
    "50 x 82 mm; Arial; vector PDF/SVG, 600 dpi PNG/LZW TIFF and ggplot RDS. Final-size layout, point overlap and font embedding require separate visual inspection."),file.path(OUT,"Figure3C_caption.txt"))
  writeLines(c(
    "Plotting implementation: ggbeeswarm 0.7.3, geom_beeswarm; horizontal orientation, swarm method, ascending priority, preserve.data.axis=TRUE, no dodging and no point omission.",
    "Source: https://github.com/eclarke/ggbeeswarm/blob/main/R/geom-beeswarm.R",
    "Examples: https://github.com/eclarke/ggbeeswarm#examples",
    "Local library: /path/to/software/R/figure3-library",
    "Original project composition; not a reproduction of a published figure."),file.path(OUT,"Figure3C_plotting_reference.txt"))
  writeLines(capture.output(sessionInfo()),file.path(OUT,"sessionInfo.txt"))
  writeLines("execution_status=success",file.path(OUT,"execution_completion.marker"))
  log_stage("task completed / done; execution_status=success")
  log_stage("please check separately plot, output completeness and biological interpretation; completion marker only indicates execution finished.")
}
status <- tryCatch({main();0L},error=function(e){log_stage("execution_status=failed; ",conditionMessage(e));1L})
sink(type="message");sink();close(con)
quit(save="no",status=status)
