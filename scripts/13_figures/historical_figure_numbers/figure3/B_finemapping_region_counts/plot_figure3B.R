#!/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/Rscript


Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
ROOT <- "/path/to/repository/scripts/13_figures/historical_figure_numbers/figure3/B_finemapping_region_counts"
INPUT <- "/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.7_define_dog_behavior_ld_loci/20260802_171634/all_locus_regions.tsv"
METADATA <- "/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.2_build_pnas_dog_behavior_locus_manifest/20260801_202130/trait_metadata.tsv"
PRODUCER <- "/path/to/repository/scripts/11_gwas_finemapping/1.7_define_dog_behavior_ld_loci/1.7_define_dog_behavior_ld_loci.py"
META_PRODUCER <- "/path/to/repository/scripts/11_gwas_finemapping/1.2_build_pnas_dog_behavior_locus_manifest/1.2_build_pnas_dog_behavior_locus_manifest.py"
OUT <- file.path(ROOT, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
if (!dir.create(OUT, recursive=TRUE, showWarnings=FALSE)) stop("cannot create new output directory")
con <- file(file.path(OUT, "plot_figure3B_regions.log"), "wt")
sink(con, split=FALSE); sink(con, type="message")
log_stage <- function(...) {cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S ]"), ..., "\n", sep=""); flush(con)}
write_tsv <- function(x, name) write.table(x, file.path(OUT,name), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
main <- function() {
  suppressPackageStartupMessages(library(ggplot2))

  font_paths <- systemfonts::match_fonts(c("Arial","Arial"),weight=c("normal","bold"))$path
  expected_fonts <- c("/path/to/fonts/Arial/Arial.TTF",
                      "/path/to/fonts/Arial/Arialbd.TTF")
  if (!all(file.exists(expected_fonts)) ||
      !identical(normalizePath(font_paths),normalizePath(expected_fonts)))
    stop("Arial regular or bold font was not resolved correctly, please check font installation")
  log_stage("[1/3] read final region list; concurrent_tasks=1; threads_per_task=1")
  log_stage("INPUT=",INPUT,"\nMETADATA=",METADATA,"\nOUT=",OUT)
  d <- read.delim(INPUT, quote="", comment.char="", check.names=FALSE)
  meta <- read.delim(METADATA, quote="", comment.char="", check.names=FALSE)
  meta <- meta[order(meta$trait), c("trait","trait_label")]


  d <- d[!duplicated(d[,c("trait","locus_id")]), ]
  tiers <- c("Genome-wide", "Suggestive", "Exploratory")
  d$tier <- ifelse(d$lead_p_lrt < 2.97e-6,tiers[1],
                   ifelse(d$lead_p_lrt < 1e-5,tiers[2],tiers[3]))
  counts <- as.data.frame(table(
    trait=factor(d$trait,levels=meta$trait),
    tier=factor(d$tier,levels=tiers)),stringsAsFactors=FALSE)
  names(counts)[3] <- "n_regions"
  counts$trait <- as.character(counts$trait)
  counts$tier <- factor(counts$tier,levels=tiers)
  counts$trait_label <- meta$trait_label[match(counts$trait,meta$trait)]
  counts$trait_display <- factor(counts$trait_label,levels=rev(meta$trait_label))
  totals <- aggregate(n_regions ~ trait, counts, sum)
  totals$trait_label <- meta$trait_label[match(totals$trait,meta$trait)]
  totals$trait_display <- factor(totals$trait_label,levels=rev(meta$trait_label))
  log_stage("[2/3] final trait region count=",nrow(d),"; trait=",nrow(meta))
  write_tsv(d,"Figure3B_region_membership.tsv")
  write_tsv(counts[,c("trait","trait_label","tier","n_regions")],"Figure3B_region_counts.tsv")
  write_tsv(totals[,c("trait","trait_label","n_regions")],"Figure3B_region_totals.tsv")
  write_tsv(data.frame(role=c("final_region_table","trait_metadata"),path=c(INPUT,METADATA),producer=c(PRODUCER,META_PRODUCER)),"input_manifest.tsv")

  labels <- tiers

  colors <- setNames(c("#0F4D92","#3775BA","#CFCECE"),tiers)

  counts$y <- nrow(meta)+1-match(counts$trait,meta$trait)
  totals$y <- nrow(meta)+1-match(totals$trait,meta$trait)
  xmax <- max(totals$n_regions)
  xbreaks <- pretty(c(0,xmax),n=4)

  xbreaks <- xbreaks[xbreaks>=0 & xbreaks<=xmax*1.18 & xbreaks==floor(xbreaks)]

  label_gap <- max(.25,.035*xmax)
  p <- ggplot(counts,aes(x=n_regions,y=y,fill=tier,group=tier)) +
    geom_col(orientation="y",width=.58,position=position_stack(reverse=TRUE)) +
    geom_text(data=totals,aes(x=n_regions+label_gap,y=y,label=n_regions),
              inherit.aes=FALSE,hjust=0,size=6.5,size.unit="pt",
              family="Arial",colour="#222222") +
    scale_fill_manual(values=colors,breaks=tiers,labels=labels,drop=FALSE) +
    scale_x_continuous(breaks=xbreaks,limits=c(0,xmax*1.18),
                       expand=expansion(mult=0)) +
    scale_y_continuous(breaks=seq_len(nrow(meta)),labels=rev(meta$trait_label),
                       limits=c(.45,nrow(meta)+.55),expand=expansion(mult=0)) +
    labs(x="Regions",y=NULL,fill=NULL,tag="B",
         title="Fine-mapping input regions") +
    theme_classic(base_size=7,base_family="Arial") +
    theme(axis.text=element_text(size=6.5,colour="#222222"),

          axis.text.y=element_blank(),
          axis.title.x=element_text(size=7.5,margin=margin(t=5)),
          axis.line.y=element_blank(),axis.ticks.y=element_blank(),
          axis.line.x=element_line(linewidth=.3),axis.ticks.x=element_line(linewidth=.3),
          legend.position="none",
          plot.title=element_text(size=8,face="bold",margin=margin(b=5)),
          plot.tag=element_text(size=10,face="bold"),plot.tag.position="topleft",
          plot.margin=margin(5,5,5,5))
  log_stage("[3/3] export plot with figure legend")
  prefix <- file.path(OUT,"Figure3B_finemapping_region_counts")
  ggsave(paste0(prefix,".pdf"),p,device=grDevices::cairo_pdf,width=45,height=82,units="mm",bg="white")
  ggsave(paste0(prefix,".png"),p,device=ragg::agg_png,width=45,height=82,units="mm",dpi=600,bg="white")
  ggsave(paste0(prefix,".svg"),p,device=svglite::svglite,width=45,height=82,units="mm",bg="white")
  ggsave(paste0(prefix,".tiff"),p,device=ragg::agg_tiff,width=45,height=82,units="mm",dpi=600,bg="white",compression="lzw")
  saveRDS(p,file.path(OUT,"Figure3B_finemapping_region_counts.rds"))
  writeLines(c("panel_width_mm=45", "panel_height_mm=82",
    "font=Arial; panel_tag=10; title=8; axis_title=7.5; axis_text/count=6.5 pt",
    "trait_labels=shared_with_panel_A; order=ascending_trait_id",
    "y_limits=0.45 to n_traits+0.55; bar_width=0.58; plot_margins=5 pt",
    "ABC planned widths=75/45/50 mm; gaps=4 mm",
    "A/B title and x-axis typography/margins matched; rendered alignment=NOT_CHECKED",
    "shared_legend=supplied_by_A_at_assembly; no_A_output_used_as_data_input",
    "PDF/SVG vector; PNG/TIFF 600 dpi; TIFF LZW; white background",
    paste("font_paths",paste(font_paths,collapse=";"))),
    file.path(OUT,"Figure3B.layout_parameters.txt"))
  writeLines(c(
    "Figure 3B. Trait-specific regions submitted to fine-mapping.",
    sprintf("Counts include %d trait-region input units across %d behavioral traits from the final 1.7 region manifest.",nrow(d),nrow(meta)),
    "Each trait-locus_id is counted once, regardless of credible-set output. Main and strict model branches are not counted separately. Genomic intervals shared across traits remain distinct trait-region units.",
    "Horizontal stacked bars show absolute region counts per trait. Dark-blue, medium-blue and light-gray segments, from left to right, denote region lead P_LRT < 2.97e-6, 2.97e-6 <= P_LRT < 1e-5, and 1e-5 <= P_LRT < 5e-5, respectively. Labels at bar ends give totals; segment counts are available in the source table. Bars are not normalized to percentages.",
    "Panel A classifies each candidate variant by its own P_LRT value, whereas this panel classifies each region by its lead P_LRT value. Variant and region counts are different units and do not form a retention percentage.",
    "All 14 metadata traits and zero category counts are retained. Rows follow panel A in ascending trait ID order. Trait names are shared with panel A and omitted from this narrow assembly panel; row identities remain in the source tables. The shared LRT P-value legend is supplied by panel A during assembly.",
    "These are descriptive counts, not counts of independent causal signals. No significance test or error bar is applied.",

    sprintf("The source manifest includes %d standalone no-proxy fallback windows and %d merged regions containing at least one fallback component; region definitions are inherited rather than recomputed here.",
            sum(d$region_source == "FALLBACK_LEAD_PLUS_MINUS_WINDOW_NO_R2_PROXY"),
            sum(d$region_source == "MERGED_COMPONENT_REGIONS" & d$fallback_component_count > 0)),
    "Panel size: 45 x 82 mm; Arial. Vector PDF/SVG, 600 dpi PNG/LZW TIFF, and ggplot RDS are exported. Final panel alignment, text legibility and font embedding require separate rendered-figure inspection."),
    file.path(OUT,"Figure3B_region_caption.txt"))
  writeLines(capture.output(sessionInfo()),file.path(OUT,"sessionInfo.txt"))
  writeLines("execution_status=success",file.path(OUT,"execution_completion.marker"))
  log_stage("task completed / done; execution_status=success")
  log_stage("please check separately actual plot, source table and biological interpretation; completion marker only indicates execution finished.")
}
status <- tryCatch({main();0L},error=function(e){log_stage("execution_status=failed; ",conditionMessage(e));1L})
sink(type="message");sink();close(con)
quit(save="no",status=status)
