#!/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/Rscript


Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
ROOT <- "/path/to/repository/scripts/13_figures/historical_figure_numbers/figure3/A_gwas_candidate_variant_counts"
UPSTREAM <- "/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.2_build_pnas_dog_behavior_locus_manifest/20260801_202130"
INPUT <- file.path(UPSTREAM, "gwas_p_lrt_candidate_variants.tsv")
METADATA <- file.path(UPSTREAM, "trait_metadata.tsv")
PRODUCER <- "/path/to/repository/scripts/11_gwas_finemapping/1.2_build_pnas_dog_behavior_locus_manifest/1.2_build_pnas_dog_behavior_locus_manifest.py"
OUT <- file.path(ROOT, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
if (!dir.create(OUT, recursive=TRUE, showWarnings=FALSE)) stop("cannot create new output directory")
con <- file(file.path(OUT, "plot_figure3A_candidate_variants.log"), "wt")
sink(con, split=FALSE); sink(con, type="message")
log_stage <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep="")
  flush(con)
}
write_tsv <- function(x, name) write.table(x, file.path(OUT,name), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
main <- function() {
  suppressPackageStartupMessages(library(ggplot2))

  font_paths <- systemfonts::match_fonts(c("Arial","Arial"),weight=c("normal","bold"))$path
  expected_fonts <- c("/path/to/fonts/Arial/Arial.TTF",
                      "/path/to/fonts/Arial/Arialbd.TTF")
  if (!all(file.exists(expected_fonts)) ||
      !identical(normalizePath(font_paths),normalizePath(expected_fonts)))
    stop("Arial regular or bold font was not resolved correctly, please check font installation")
  log_stage("[1/3] read upstream candidate; concurrent_tasks=1; threads_per_task=1")
  log_stage("INPUT=",INPUT,"\nMETADATA=",METADATA,"\nOUT=",OUT)
  d <- read.delim(INPUT, quote="", comment.char="", check.names=FALSE)
  meta <- read.delim(METADATA, quote="", comment.char="", check.names=FALSE)
  meta <- meta[order(meta$trait), c("trait","trait_label")]
  input_rows <- nrow(d)

  d <- d[is.finite(d$p_lrt) & d$p_lrt >= 0 & d$p_lrt < 5e-5, ]
  d$allele_set <- vapply(seq_len(nrow(d)), function(i)
    paste(sort(unique(toupper(c(d$allele1[i],d$allele0[i])))),collapse="/"), character(1))
  d$variant_key <- paste(d$chrom,d$pos,d$allele_set,sep=":")

  d <- d[order(d$trait,d$p_lrt,d$variant_key), ]
  d <- d[!duplicated(d[,c("trait","variant_key")]), ]
  tiers <- c("Genome-wide", "Suggestive", "Exploratory")
  d$tier <- ifelse(d$p_lrt < 2.97e-6,tiers[1],ifelse(d$p_lrt < 1e-5,tiers[2],tiers[3]))
  counts <- as.data.frame(table(
    trait=factor(d$trait,levels=meta$trait),
    tier=factor(d$tier,levels=tiers)),stringsAsFactors=FALSE)
  names(counts)[3] <- "n_variants"
  counts$trait <- as.character(counts$trait)
  counts$tier <- factor(counts$tier,levels=tiers)
  counts$trait_label <- meta$trait_label[match(counts$trait,meta$trait)]
  counts$trait_display <- factor(counts$trait_label,levels=rev(meta$trait_label))
  totals <- aggregate(n_variants ~ trait, counts, sum)
  totals$trait_label <- meta$trait_label[match(totals$trait,meta$trait)]
  totals$trait_display <- factor(totals$trait_label,levels=rev(meta$trait_label))
  log_stage("[2/3] input row=",input_rows,"; deduplicate candidate=",nrow(d),"; trait=",nrow(meta))
  write_tsv(d,"Figure3A_candidate_membership.tsv")
  write_tsv(counts[,c("trait","trait_label","tier","n_variants")],"Figure3A_candidate_counts.tsv")
  write_tsv(totals[,c("trait","trait_label","n_variants")],"Figure3A_candidate_totals.tsv")
  write_tsv(data.frame(role=c("candidate_table","trait_metadata"),path=c(INPUT,METADATA),producer=PRODUCER),"input_manifest.tsv")

  labels <- c("Genome-wide: P_LRT < 2.97e-6",
              "Suggestive: 2.97e-6 <= P_LRT < 1e-5",
              "Exploratory: 1e-5 <= P_LRT < 5e-5")

  colors <- setNames(c("#0F4D92","#3775BA","#CFCECE"),tiers)
  outlines <- setNames(c("#0F4D92","#3775BA","#767676"),tiers)

  counts$y <- nrow(meta)+1-match(counts$trait,meta$trait)
  totals$y <- nrow(meta)+1-match(totals$trait,meta$trait)
  xmax <- max(totals$n_variants)
  xbreaks <- pretty(c(0,xmax),n=4)
  xbreaks <- xbreaks[xbreaks>=0 & xbreaks<=xmax & xbreaks==floor(xbreaks)]

  label_gap <- max(1,.025*xmax)
  p <- ggplot(counts,aes(x=n_variants,y=y,fill=tier,group=tier)) +
    geom_col(orientation="y",width=.58,position=position_stack(reverse=TRUE)) +
    geom_text(data=totals,aes(x=n_variants+label_gap,y=y,label=n_variants),
              inherit.aes=FALSE,hjust=0,size=6.5,size.unit="pt",
              family="Arial",colour="#222222") +
    scale_fill_manual(values=colors,breaks=tiers,labels=labels,drop=FALSE) +
    scale_x_continuous(breaks=xbreaks,limits=c(0,xmax*1.16),
                       expand=expansion(mult=0)) +
    scale_y_continuous(breaks=seq_len(nrow(meta)),labels=rev(meta$trait_label),
                       limits=c(.45,nrow(meta)+.55),expand=expansion(mult=0)) +
    labs(x="Candidate variants",y=NULL,fill=NULL,tag="A",
         title="GWAS candidates") +
    theme_classic(base_size=7,base_family="Arial") +
    theme(axis.text=element_text(size=6.5,colour="#222222"),
          axis.text.y=element_text(lineheight=.95,margin=margin(r=5)),
          axis.title.x=element_text(size=7.5,margin=margin(t=5)),
          axis.line.y=element_blank(),axis.ticks.y=element_blank(),
          axis.line.x=element_line(linewidth=.3),axis.ticks.x=element_line(linewidth=.3),
          legend.position="none",
          plot.title=element_text(size=8,face="bold",margin=margin(b=5)),
          plot.tag=element_text(size=10,face="bold"),plot.tag.position="topleft",
          plot.margin=margin(5,5,5,5))
  log_stage("[3/3] export plot with figure legend")
  prefix <- file.path(OUT,"Figure3A_candidate_variant_counts")
  ggsave(paste0(prefix,".pdf"),p,device=grDevices::cairo_pdf,width=75,height=82,units="mm",bg="white")
  ggsave(paste0(prefix,".png"),p,device=ragg::agg_png,width=75,height=82,units="mm",dpi=600,bg="white")
  ggsave(paste0(prefix,".svg"),p,device=svglite::svglite,width=75,height=82,units="mm",bg="white")
  ggsave(paste0(prefix,".tiff"),p,device=ragg::agg_tiff,width=75,height=82,units="mm",dpi=600,bg="white",compression="lzw")
  saveRDS(p,file.path(OUT,"Figure3A_candidate_variant_counts.rds"))

  legend_data <- data.frame(tier=factor(tiers,levels=tiers),y=3:1,label=labels)
  legend_plot <- ggplot(legend_data,aes(y=y)) +

    geom_point(aes(x=0.055,fill=tier,colour=tier),shape=21,size=1.6,stroke=.18) +

    geom_text(aes(x=0.090,label=label),hjust=0,size=7,size.unit="pt",parse=FALSE,
              family="Arial",colour="#222222") +
    scale_fill_manual(values=colors,drop=FALSE) +
    scale_colour_manual(values=outlines,drop=FALSE) +
    scale_x_continuous(limits=c(0,1),expand=expansion(mult=0)) +
    scale_y_continuous(limits=c(.5,3.5),expand=expansion(mult=0)) +
    theme_void(base_size=7,base_family="Arial") +
    theme(legend.position="none",plot.margin=margin(2,2,2,2))
  ggsave(file.path(OUT,"Figure3ABC_shared_legend.pdf"),legend_plot,
         device=grDevices::cairo_pdf,width=178,height=12,units="mm",bg="white")
  ggsave(file.path(OUT,"Figure3ABC_shared_legend.png"),legend_plot,
         device=ragg::agg_png,width=178,height=12,units="mm",dpi=600,bg="white")
  ggsave(file.path(OUT,"Figure3ABC_shared_legend.svg"),legend_plot,
         device=svglite::svglite,width=178,height=12,units="mm",bg="white")
  ggsave(file.path(OUT,"Figure3ABC_shared_legend.tiff"),legend_plot,
         device=ragg::agg_tiff,width=178,height=12,units="mm",dpi=600,bg="white",compression="lzw")
  saveRDS(legend_plot,file.path(OUT,"Figure3ABC_shared_legend.rds"))
  writeLines(c("panel_width_mm=75", "panel_height_mm=82",
    "shared_legend_width_mm=178", "shared_legend_height_mm=12",
    "font=Arial; panel_tag=10; title=8; axis_title=7.5; axis_text/count=6.5; legend=7 pt",
    "ABC planned widths=75/45/50 mm; gaps=4 mm; A retains full trait labels",
    "ABC row alignment=NOT_CHECKED; B/C must be aligned to the final rendered A data panel",
    "PDF/SVG vector; PNG/TIFF 600 dpi; TIFF LZW; white background",
    paste("font_paths",paste(font_paths,collapse=";"))),
    file.path(OUT,"Figure3A.layout_parameters.txt"))
  writeLines(c(
    "Figure 3A. GWAS candidate variants preceding fine-mapping region definition.",
    sprintf("The plot includes %d trait-variant records across %d behavioral traits selected from complete GWAS summaries using variant LRT P < 5e-5.",nrow(d),nrow(meta)),
    "Horizontal stacked bars show absolute candidate counts per trait. Dark-blue, medium-blue and light-gray segments, from left to right, denote mutually exclusive categories based on each candidate variant's own LRT P value: genome-wide P_LRT < 2.97e-6; suggestive 2.97e-6 <= P_LRT < 1e-5; exploratory 1e-5 <= P_LRT < 5e-5. Labels at bar ends give totals; segment counts are available in the source table. Bars are not normalized to percentages.",
    "Counts precede LD compatibility filtering and include candidates absent from the final LD reference. They are not counts of all tested GWAS variants, independent signals, or LD-defined regions.",
    "A variant is deduplicated within a trait by chromosome, position and unordered allele set; associations with different traits remain separate records.",
    "The shared legend denotes LRT P-value thresholds. Panel A classifies each variant by its own LRT P value; panels B and C classify results by the region lead LRT P value. Variant and region counts have different units and must not be interpreted as a retention percentage.",
    "Colors and trait ordering match the regional-count panel. The separate shared legend uses equally sized colored circles for the association tiers; the same tier colors fill the bars in panels A and B and the points in panel C. All metadata traits and zero tier counts are retained. No significance tests or uncertainty intervals are applied to these descriptive counts.",
    "Panel: 75 x 82 mm; separate shared legend: 178 x 12 mm. Arial. The plot and legend are saved as vector PDF/SVG, 600 dpi PNG/LZW TIFF, and RDS. The shared legend uses plain scientific notation at 7 pt to avoid undersized superscripts. Final assembly alignment and PDF text/collision checks remain separate post-render reviews."),file.path(OUT,"Figure3A_candidate_caption.txt"))
  writeLines(capture.output(sessionInfo()),file.path(OUT,"sessionInfo.txt"))
  writeLines("execution_status=success",file.path(OUT,"execution_completion.marker"))
  log_stage("task completed / done; execution_status=success")
  log_stage("please check separately actual plot, source table and biological interpretation; completion marker only indicates execution finished.")
}
status <- tryCatch({main();0L},error=function(e){log_stage("execution_status=failed; ",conditionMessage(e));1L})
sink(type="message");sink();close(con)
quit(save="no",status=status)
