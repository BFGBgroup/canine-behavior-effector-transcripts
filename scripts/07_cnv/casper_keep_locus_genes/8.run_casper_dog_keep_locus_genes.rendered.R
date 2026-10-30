#!/usr/bin/env Rscript


suppressPackageStartupMessages({
  library(CaSpER)
  library(GenomicRanges)
  library(IRanges)
  library(signal)
  library(HMMcopy)
  library(mclust)
  library(pheatmap)
  library(RColorBrewer)
})

options(stringsAsFactors = FALSE)


args <- commandArgs(trailingOnly = TRUE)

tissue <- "cerebellum"
cnv.scale <- 3L
loh.scale <- 3L
window.length <- 50L
length.iterations <- 50L
baf.mode <- "hybrid"
min.shift.threshold <- 0.10
high.rate.threshold <- 0.50
high.rate.fallback.threshold <- 0.15

i <- 1L
while (i <= length(args)) {
  if (args[i] == "--tissue") {
    tissue <- args[i + 1L]
    i <- i + 2L
  } else if (args[i] == "--cnv-scale") {
    cnv.scale <- as.integer(args[i + 1L])
    i <- i + 2L
  } else if (args[i] == "--loh-scale") {
    loh.scale <- as.integer(args[i + 1L])
    i <- i + 2L
  } else if (args[i] == "--window-length") {
    window.length <- as.integer(args[i + 1L])
    i <- i + 2L
  } else if (args[i] == "--length-iterations") {
    length.iterations <- as.integer(args[i + 1L])
    i <- i + 2L
  } else if (args[i] == "--baf-mode") {
    baf.mode <- tolower(args[i + 1L])
    i <- i + 2L
  } else if (args[i] == "--min-shift-threshold") {
    min.shift.threshold <- as.numeric(args[i + 1L])
    i <- i + 2L
  } else if (args[i] == "--high-rate-threshold") {
    high.rate.threshold <- as.numeric(args[i + 1L])
    i <- i + 2L
  } else if (args[i] == "--high-rate-fallback-threshold") {
    high.rate.fallback.threshold <- as.numeric(args[i + 1L])
    i <- i + 2L
  } else {
    stop("Unknown argument: ", args[i])
  }
}

baf.mode <- switch(
  baf.mode,
  "hybrid" = "hybrid",
  "off" = "off",
  "none" = "off",
  "expression-only" = "off",
  "expression_only" = "off",
  "auto" = "auto",
  stop("Unsupported --baf-mode: ", baf.mode, " (use hybrid, off, or auto)")
)


expr_file <- "/path/to/project/legacy/five-dogs/5.cnv/5.casper_input/bulk_rna.cnv_expression.keep_locus_genes.all_genes.log2_norm.primary_chr.casper_ordered.tsv"
annot_file <- "/path/to/project/legacy/five-dogs/5.cnv/5.casper_input/CanFam3.1.keep_locus_genes.all_genes.gene_annotation.primary_chr.for_casper.tsv"
baf_dir <- "/path/to/project/legacy/five-dogs/5.cnv/6.baf_signal/keep_locus_genes.all_genes/primary/per_sample_baf"
fai_file <- "/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna.fai"
out_root <- "/path/to/project/analysis_20260609/5-dogs/4.cnv-rna/8.casper_dog_keep_locus_genes-v2"

out_tissue <- file.path(out_root, tissue)
dir.create(out_tissue, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_tissue, paste0("7.run_casper_dog.", tissue, ".log"))
log_con <- file(log_file, open = "w")

log_msg <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  x <- paste0("[", timestamp, "] ", msg)
  cat(x, "\n")
  cat(x, "\n", file = log_con)
  flush(log_con)
}

on.exit({
  try(close(log_con), silent = TRUE)
}, add = TRUE)


is_dog_autosome <- function(chr) {
  chr %in% paste0("chr", 1:38)
}

chr_order_key <- function(chr) {
  x <- sub("^chr", "", as.character(chr))
  suppressWarnings(num <- as.integer(x))
  ifelse(!is.na(num), num, 999L)
}

map_raw_state <- function(state) {
  map <- c("del", "del", "neutral", "amp", "amp")
  map[pmin(pmax(as.integer(state), 1L), 5L)]
}

map_expression_only_state <- function(state) {
  out <- rep("neutral", length(state))
  out[as.integer(state) == 1L] <- "del"
  out[as.integer(state) == 5L] <- "amp"
  out
}

prepareAnnotationDog <- function(annot_df) {
  log_msg("Preparing dog annotation: keeping chr1-chr38 only")

  required <- c("Gene", "GeneSymbol", "Chr", "start", "end")
  missing <- setdiff(required, colnames(annot_df))
  if (length(missing) > 0L) {
    stop("Annotation missing columns: ", paste(missing, collapse = ", "))
  }

  annot_df <- annot_df[is_dog_autosome(annot_df$Chr), , drop = FALSE]
  annot_df$start <- as.numeric(annot_df$start)
  annot_df$end <- as.numeric(annot_df$end)
  annot_df$Position <- (annot_df$start + annot_df$end) / 2

  annot_df$chr_key <- chr_order_key(annot_df$Chr)
  annot_df <- annot_df[order(annot_df$chr_key, annot_df$start, annot_df$end, annot_df$Gene), , drop = FALSE]
  annot_df$chr_key <- NULL

  annot_df$cytoband <- annot_df$Chr
  annot_df$isCentromer <- "no"
  annot_df$new_positions <- ave(seq_len(nrow(annot_df)), annot_df$Chr, FUN = seq_along)
  rownames(annot_df) <- annot_df$Gene

  log_msg(paste0("Annotation autosome genes: ", nrow(annot_df)))
  annot_df
}

getDogCytoband <- function(fai_file) {
  if (!file.exists(fai_file)) {
    stop("CanFam3.1 FASTA index not found: ", fai_file)
  }

  fai <- read.delim(fai_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(fai)[1:2] <- c("chr", "size")
  fai <- fai[is_dog_autosome(fai$chr), c("chr", "size"), drop = FALSE]
  fai$chr_key <- chr_order_key(fai$chr)
  fai <- fai[order(fai$chr_key), , drop = FALSE]

  chrs <- paste0("chr", 1:38)
  missing <- setdiff(chrs, fai$chr)
  if (length(missing) > 0L) {
    stop("FASTA index missing dog autosome(s): ", paste(missing, collapse = ", "))
  }

  fai <- fai[match(chrs, fai$chr), , drop = FALSE]
  chr_sizes <- as.numeric(fai$size)

  data.frame(
    V1 = chrs,
    V2 = chrs,
    V3 = chr_sizes,
    stringsAsFactors = FALSE
  )
}

smoothBafListDog <- function(loh_list, loh.scale, n = 50L, scale.iteration = 50L) {
  out <- list()

  for (sample_name in names(loh_list)) {
    maf_temp <- loh_list[[sample_name]]

    for (scale_i in seq_len(loh.scale)) {
      window <- n + (scale_i - 1L) * scale.iteration
      maf_2 <- NULL

      for (m in 1:38) {
        chr_str <- paste0("chr", m)
        chrBAF <- maf_temp[as.character(maf_temp$chr) == chr_str, , drop = FALSE]
        if (nrow(chrBAF) == 0L) next

        chrBAF <- chrBAF[order(as.numeric(chrBAF$pos)), , drop = FALSE]
        data_smoothed <- as.numeric(chrBAF$dev)

        if (length(data_smoothed) >= (window + 1L)) {
          data_smoothed <- round(
            signal::filter(signal::MedianFilter(window + 1L), data_smoothed),
            digits = 4
          )
        }

        chrBAF$dev <- as.numeric(data_smoothed)
        maf_2 <- rbind(maf_2, chrBAF)
      }

      if (is.null(maf_2) || nrow(maf_2) == 0L) {
        maf_2 <- maf_temp
      }

      maf_temp <- maf_2
    }

    out[[sample_name]] <- maf_temp
  }

  out
}

callBafShiftGMMDog <- function(
  loh_smoothed,
  min_snps_for_gmm = 30L,
  fallback_quantile = 0.75,
  min_shift_threshold = 0.10,
  high_rate_threshold = 0.50,
  high_rate_fallback_threshold = 0.15
) {
  out <- list()
  qc <- data.frame()

  for (sample_name in names(loh_smoothed)) {
    d <- loh_smoothed[[sample_name]]
    if (is.null(d) || nrow(d) == 0L) {
      d <- data.frame(
        chr = character(),
        pos = integer(),
        dev = numeric(),
        loh_call = integer(),
        loh_threshold = numeric()
      )
      out[[sample_name]] <- d
      next
    }

    d$loh_call <- 0L
    d$loh_threshold <- NA_real_

    for (m in 1:38) {
      chr_str <- paste0("chr", m)
      idx <- which(as.character(d$chr) == chr_str & is.finite(d$dev))
      if (length(idx) == 0L) next

      vals <- as.numeric(d$dev[idx])
      vals <- vals[is.finite(vals)]

      threshold_raw <- NA_real_
      threshold <- NA_real_
      method <- "fallback_quantile"

      if (length(vals) >= min_snps_for_gmm && length(unique(vals)) >= 3L) {
        fit <- tryCatch(
          mclust::Mclust(vals, G = 2, verbose = FALSE),
          error = function(e) NULL
        )

        if (!is.null(fit) && !is.null(fit$classification)) {
          means <- fit$parameters$mean
          shift_class <- which.max(means)
          shifted_vals <- vals[fit$classification == shift_class]

          if (length(shifted_vals) > 0L) {
            threshold_raw <- suppressWarnings(min(shifted_vals, na.rm = TRUE))
            method <- "GMM"
          }
        }
      }

      if (!is.finite(threshold_raw)) {
        threshold_raw <- as.numeric(stats::quantile(vals, fallback_quantile, na.rm = TRUE, names = FALSE))
        method <- "fallback_quantile"
      }

      threshold <- max(threshold_raw, min_shift_threshold)
      tmp_call <- as.integer(as.numeric(d$dev[idx]) >= threshold)
      tmp_rate <- mean(tmp_call == 1L, na.rm = TRUE)

      if (is.finite(tmp_rate) && tmp_rate > high_rate_threshold) {
        threshold <- high_rate_fallback_threshold
        method <- paste0(method, "_highRateFallback")
        tmp_call <- as.integer(as.numeric(d$dev[idx]) >= threshold)
        tmp_rate <- mean(tmp_call == 1L, na.rm = TRUE)
      } else if (threshold > threshold_raw) {
        method <- paste0(method, "_floor")
      }

      d$loh_threshold[idx] <- threshold
      d$loh_call[idx] <- tmp_call

      qc <- rbind(
        qc,
        data.frame(
          sample_id = sample_name,
          chr = chr_str,
          n_snps = length(idx),
          threshold_raw = threshold_raw,
          threshold = threshold,
          method = method,
          called_shift_sites = sum(d$loh_call[idx] == 1L, na.rm = TRUE),
          shift_rate = mean(d$loh_call[idx] == 1L, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      )
    }

    out[[sample_name]] <- d
  }

  list(loh = out, qc = qc)
}

PerformSegmentationWithHMMDog <- function(object, cnv.scale) {
  log_msg(paste0("Performing HMM segmentation: cnv.scale=", cnv.scale))

  object <- generateParam(object, cnv.scale = cnv.scale)

  data <- object@control.normalized[[cnv.scale]]
  annotation <- object@annotation.filt[]
  annotation$cytoband <- annotation$Chr

  keep <- annotation$Chr %in% names(which(table(annotation$Chr) > 1))
  annotation <- annotation[keep, , drop = FALSE]
  data <- data[keep, , drop = FALSE]

  log_msg(paste0("HMM genes: ", nrow(annotation), "; samples: ", ncol(data)))

  segments <- NULL

  for (i in seq_len(ncol(data))) {
    sample_id <- colnames(data)[i]

    rdata <- GRanges(
      ranges = IRanges(start = annotation$start, end = annotation$end),
      seqnames = annotation$Chr,
      copy = data[, i],
      chr = annotation$Chr,
      space = annotation$Chr
    )

    rdata <- data.frame(rdata)
    rdata$chr <- as.factor(rdata$chr)

    hmm.segments <- HMMcopy::HMMsegment(correctOut = rdata, param = object@hmmparam, verbose = FALSE)
    segs_i <- hmm.segments$segs

    segs_i$state <- as.integer(segs_i$state)
    segs_i$rawCNVState <- map_raw_state(segs_i$state)
    segs_i$exprOnlyCNVState <- map_expression_only_state(segs_i$state)
    segs_i$states2 <- segs_i$rawCNVState
    segs_i$finalCNVState <- segs_i$rawCNVState
    segs_i$baf_decision <- "not_checked"
    segs_i$baf_num_snps <- 0L
    segs_i$baf_support_frac <- NA_real_
    segs_i$baf_median_shift <- NA_real_

    segments <- rbind(segments, data.frame(ID = sample_id, segs_i))
  }

  if (is.null(segments) || nrow(segments) == 0L) {
    stop("No HMM segments generated")
  }

  cytoband <- object@cytoband
  arms <- cytoband$V1
  arm_sizes <- cytoband$V3

  object@segments <- segments
  object@segments$event_scale <- ""
  object@segments$size <- as.numeric(object@segments$end) - as.numeric(object@segments$start)
  object@segments$num.marks <- 0L
  object@segments$arm.size <- NA_real_
  object@segments$segment.ratio <- NA_real_

  for (i in seq_len(nrow(object@segments))) {
    chr_i <- as.character(object@segments$chr[i])

    ind <- which(
      annotation$Chr == chr_i &
        annotation$Position >= as.numeric(object@segments$start[i]) &
        annotation$Position <= as.numeric(object@segments$end[i])
    )

    object@segments$num.marks[i] <- length(ind)

    pair_arm_size <- arm_sizes[match(chr_i, arms)]
    object@segments$arm.size[i] <- pair_arm_size
    object@segments$segment.ratio[i] <- object@segments$size[i] / pair_arm_size

    if (is.finite(object@segments$segment.ratio[i]) &&
        object@segments$segment.ratio[i] >= 1/3) {
      object@segments$event_scale[i] <- "large_scale"
    } else if (is.finite(object@segments$segment.ratio[i]) &&
               object@segments$segment.ratio[i] >= 1/10) {
      object@segments$event_scale[i] <- "focal"
    } else {
      object@segments$event_scale[i] <- "small"
    }
  }

  object
}

applyExpressionOnlyMode <- function(segments) {
  segments$states2 <- segments$exprOnlyCNVState
  segments$finalCNVState <- segments$exprOnlyCNVState
  segments$baf_decision <- "baf_disabled_expression_only"
  segments
}

correctSegmentsByBafDog <- function(segments, loh_sample) {
  if (is.null(segments) || nrow(segments) == 0L) return(segments)
  if (is.null(loh_sample) || nrow(loh_sample) == 0L) {
    segments$states2 <- segments$exprOnlyCNVState
    segments$finalCNVState <- segments$exprOnlyCNVState
    segments$baf_decision[segments$state %in% c(2L, 4L)] <- "missing_baf_expression_only"
    segments$baf_decision[segments$state == 3L] <- "missing_baf_cnloh_unchecked"
    segments$baf_decision[segments$state %in% c(1L, 5L)] <- "no_baf_needed"
    return(segments)
  }

  if (!"loh_call" %in% colnames(loh_sample)) {
    warning("loh_sample has no loh_call column; falling back to expression-only interpretation")
    return(applyExpressionOnlyMode(segments))
  }

  segments$states2 <- segments$rawCNVState
  segments$finalCNVState <- segments$rawCNVState
  segments$baf_decision <- ifelse(segments$state %in% c(2L, 3L, 4L), "pending_baf_check", "no_baf_needed")

  for (i in seq_len(nrow(segments))) {
    chr_str <- as.character(segments$chr[i])

    in_seg <- as.character(loh_sample$chr) == chr_str &
      as.numeric(loh_sample$pos) >= as.numeric(segments$start[i]) &
      as.numeric(loh_sample$pos) <= as.numeric(segments$end[i])

    n_snp <- sum(in_seg, na.rm = TRUE)
    segments$baf_num_snps[i] <- n_snp

    if (n_snp > 0L) {
      segments$baf_median_shift[i] <- median(as.numeric(loh_sample$dev[in_seg]), na.rm = TRUE)
      segments$baf_support_frac[i] <- mean(as.numeric(loh_sample$loh_call[in_seg]) == 1L, na.rm = TRUE)
    }

    if (as.integer(segments$state[i]) == 3L) {
      loh_shift_thr <- suppressWarnings(median(as.numeric(loh_sample$loh_threshold[in_seg]), na.rm = TRUE))
      if (!is.finite(loh_shift_thr)) {
        segments$baf_decision[i] <- "insufficient_baf_sites_cnloh_neutral"
        next
      }

      if (is.finite(segments$baf_median_shift[i]) &&
          segments$baf_median_shift[i] > loh_shift_thr) {
        segments$states2[i] <- "cnloh"
        segments$finalCNVState[i] <- "cnloh"
        segments$baf_decision[i] <- "cnloh_supported_by_baf"
      } else {
        segments$baf_decision[i] <- "no_cnloh_baf_support"
      }
      next
    }

    if (!(as.integer(segments$state[i]) %in% c(2L, 4L))) next

    loh_shift_thr <- suppressWarnings(median(as.numeric(loh_sample$loh_threshold[in_seg]), na.rm = TRUE))
    if (!is.finite(loh_shift_thr)) {
      segments$states2[i] <- segments$exprOnlyCNVState[i]
      segments$finalCNVState[i] <- segments$exprOnlyCNVState[i]
      segments$baf_decision[i] <- "insufficient_baf_sites_expression_only"
      next
    }

    if (is.finite(segments$baf_median_shift[i]) &&
        segments$baf_median_shift[i] > loh_shift_thr) {
      segments$baf_decision[i] <- "supported_by_baf"
    } else {
      segments$states2[i] <- segments$exprOnlyCNVState[i]
      segments$finalCNVState[i] <- segments$exprOnlyCNVState[i]
      segments$baf_decision[i] <- "unsupported_by_baf_expression_only"
    }
  }

  segments
}

extractEventsDog <- function(segments, cytoband, state_col = "finalCNVState", type = c("amp", "del")) {
  type <- match.arg(type)

  sample_ids <- unique(as.character(segments$ID))
  arms <- as.character(cytoband$V1)
  arm_sizes <- as.numeric(cytoband$V3)

  mat_large <- matrix(0L, nrow = length(sample_ids), ncol = length(arms), dimnames = list(sample_ids, arms))
  mat_focal <- matrix(0L, nrow = length(sample_ids), ncol = length(arms), dimnames = list(sample_ids, arms))

  for (sample_id in sample_ids) {
    sub_segments <- segments[as.character(segments$ID) == sample_id, , drop = FALSE]

    for (chr_i in unique(as.character(sub_segments$chr))) {
      idx <- which(as.character(sub_segments$chr) == chr_i &
                     as.character(sub_segments[[state_col]]) == type)
      if (length(idx) == 0L) next

      chr_size <- arm_sizes[match(chr_i, arms)]
      if (!is.finite(chr_size)) next

      arm_size <- sum(as.numeric(sub_segments$size[idx]), na.rm = TRUE)

      if (arm_size >= chr_size * (1/3)) {
        mat_large[sample_id, chr_i] <- 1L
      } else if (arm_size >= chr_size * (1/10)) {
        mat_focal[sample_id, chr_i] <- 1L
      }
    }
  }

  list(large = mat_large, focal = mat_focal)
}

generateLargeScaleEventsDog <- function(object, state_col = "finalCNVState") {
  amp <- extractEventsDog(object@segments, object@cytoband, state_col = state_col, type = "amp")
  del <- extractEventsDog(object@segments, object@cytoband, state_col = state_col, type = "del")

  sample_ids <- rownames(amp$large)
  final <- data.frame(
    LargeScaleAmp = "",
    FocalAmp = "",
    LargeScaleAmpNum = 0L,
    FocalAmpNum = 0L,
    LargeScaleDel = "",
    FocalDel = "",
    LargeScaleDelNum = 0L,
    FocalDelNum = 0L,
    stringsAsFactors = FALSE
  )
  final <- final[rep(1, length(sample_ids)), , drop = FALSE]
  rownames(final) <- sample_ids

  for (sample_id in sample_ids) {
    final[sample_id, "LargeScaleAmp"] <- paste(colnames(amp$large)[amp$large[sample_id, ] == 1L], collapse = " ")
    final[sample_id, "FocalAmp"] <- paste(colnames(amp$focal)[amp$focal[sample_id, ] == 1L], collapse = " ")
    final[sample_id, "LargeScaleAmpNum"] <- sum(amp$large[sample_id, ] == 1L)
    final[sample_id, "FocalAmpNum"] <- sum(amp$focal[sample_id, ] == 1L)

    final[sample_id, "LargeScaleDel"] <- paste(colnames(del$large)[del$large[sample_id, ] == 1L], collapse = " ")
    final[sample_id, "FocalDel"] <- paste(colnames(del$focal)[del$focal[sample_id, ] == 1L], collapse = " ")
    final[sample_id, "LargeScaleDelNum"] <- sum(del$large[sample_id, ] == 1L)
    final[sample_id, "FocalDelNum"] <- sum(del$focal[sample_id, ] == 1L)
  }

  object@large.scale.cnv.events <- final
  object
}

extractLargeScaleEventsDog <- function(final.objects, thr = 0.75) {
  sample_names <- rownames(final.objects[[1]]@large.scale.cnv.events)
  chrs <- paste0("chr", 1:38)

  amp_count <- matrix(0L, nrow = length(sample_names), ncol = length(chrs), dimnames = list(sample_names, chrs))
  del_count <- matrix(0L, nrow = length(sample_names), ncol = length(chrs), dimnames = list(sample_names, chrs))

  for (object in final.objects) {
    ev <- object@large.scale.cnv.events

    for (sample_id in sample_names) {
      amp_chr <- unlist(strsplit(as.character(ev[sample_id, "LargeScaleAmp"]), "\\s+"))
      del_chr <- unlist(strsplit(as.character(ev[sample_id, "LargeScaleDel"]), "\\s+"))

      amp_chr <- intersect(amp_chr[amp_chr != ""], chrs)
      del_chr <- intersect(del_chr[del_chr != ""], chrs)

      if (length(amp_chr) > 0L) amp_count[sample_id, amp_chr] <- amp_count[sample_id, amp_chr] + 1L
      if (length(del_chr) > 0L) del_count[sample_id, del_chr] <- del_count[sample_id, del_chr] + 1L
    }
  }

  min_agree <- ceiling(thr * length(final.objects))

  mat <- matrix(0L, nrow = length(sample_names), ncol = length(chrs), dimnames = list(sample_names, chrs))
  amp_call <- amp_count >= min_agree
  del_call <- del_count >= min_agree
  conflict <- amp_call & del_call

  mat[amp_call] <- 1L
  mat[del_call] <- -1L
  mat[conflict] <- 0L

  mat
}

collapse_state_calls <- function(calls) {
  states <- c("amp", "del", "cnloh")
  called <- states[calls[states]]
  if (length(called) == 0L) "neutral" else paste(called, collapse = "+")
}

extractLargeScaleStateMatrixDog <- function(final.objects, state_col = "finalCNVState", thr = 0.75) {
  sample_names <- unique(as.character(final.objects[[1]]@segments$ID))
  chrs <- paste0("chr", 1:38)
  states <- c("amp", "del", "cnloh")
  min_agree <- ceiling(thr * length(final.objects))

  counts <- lapply(
    states,
    function(x) matrix(0L, nrow = length(sample_names), ncol = length(chrs), dimnames = list(sample_names, chrs))
  )
  names(counts) <- states

  cytoband <- final.objects[[1]]@cytoband
  chr_sizes <- setNames(as.numeric(cytoband$V3), as.character(cytoband$V1))

  for (object in final.objects) {
    segs <- object@segments

    for (sample_id in sample_names) {
      ss <- segs[as.character(segs$ID) == sample_id, , drop = FALSE]
      if (nrow(ss) == 0L) next

      for (chr_i in intersect(unique(as.character(ss$chr)), chrs)) {
        chr_size <- chr_sizes[[chr_i]]
        if (!is.finite(chr_size)) next

        for (state_i in states) {
          idx <- which(as.character(ss$chr) == chr_i & as.character(ss[[state_col]]) == state_i)
          if (length(idx) == 0L) next

          state_size <- sum(as.numeric(ss$size[idx]), na.rm = TRUE)
          if (state_size >= chr_size * (1/3)) {
            counts[[state_i]][sample_id, chr_i] <- counts[[state_i]][sample_id, chr_i] + 1L
          }
        }
      }
    }
  }

  out <- matrix("neutral", nrow = length(sample_names), ncol = length(chrs), dimnames = list(sample_names, chrs))
  for (sample_id in sample_names) {
    for (chr_i in chrs) {
      calls <- vapply(states, function(state_i) counts[[state_i]][sample_id, chr_i] >= min_agree, logical(1))
      out[sample_id, chr_i] <- collapse_state_calls(calls)
    }
  }

  out
}

geneLevelCNVDog <- function(final.objects, annot_data, state_col = "finalCNVState", thr = 0.75) {
  sample_names <- unique(unlist(lapply(final.objects, function(o) unique(o@segments$ID))))
  min_agree <- ceiling(thr * length(final.objects))

  amp_count <- matrix(0L, nrow = nrow(annot_data), ncol = length(sample_names), dimnames = list(annot_data$Gene, sample_names))
  del_count <- matrix(0L, nrow = nrow(annot_data), ncol = length(sample_names), dimnames = list(annot_data$Gene, sample_names))

  for (object in final.objects) {
    segs <- object@segments

    for (sample_id in sample_names) {
      ss <- segs[as.character(segs$ID) == sample_id, , drop = FALSE]
      if (nrow(ss) == 0L) next

      for (i in seq_len(nrow(ss))) {
        chr_i <- as.character(ss$chr[i])
        idx <- which(
          annot_data$Chr == chr_i &
            annot_data$Position >= as.numeric(ss$start[i]) &
            annot_data$Position <= as.numeric(ss$end[i])
        )

        if (length(idx) == 0L) next

        state_i <- as.character(ss[[state_col]][i])
        if (state_i == "amp") {
          amp_count[idx, sample_id] <- amp_count[idx, sample_id] + 1L
        } else if (state_i == "del") {
          del_count[idx, sample_id] <- del_count[idx, sample_id] + 1L
        }
      }
    }
  }

  mat <- matrix(0L, nrow = nrow(annot_data), ncol = length(sample_names), dimnames = list(annot_data$Gene, sample_names))
  amp_call <- amp_count >= min_agree
  del_call <- del_count >= min_agree
  conflict <- amp_call & del_call

  mat[amp_call] <- 1L
  mat[del_call] <- -1L
  mat[conflict] <- 0L

  mat
}

geneLevelCNVStateDog <- function(final.objects, annot_data, state_col = "finalCNVState", thr = 0.75) {
  sample_names <- unique(unlist(lapply(final.objects, function(o) unique(o@segments$ID))))
  min_agree <- ceiling(thr * length(final.objects))
  states <- c("amp", "del", "cnloh")

  counts <- lapply(
    states,
    function(x) matrix(0L, nrow = nrow(annot_data), ncol = length(sample_names), dimnames = list(annot_data$Gene, sample_names))
  )
  names(counts) <- states

  for (object in final.objects) {
    segs <- object@segments

    for (sample_id in sample_names) {
      ss <- segs[as.character(segs$ID) == sample_id, , drop = FALSE]
      if (nrow(ss) == 0L) next

      for (i in seq_len(nrow(ss))) {
        chr_i <- as.character(ss$chr[i])
        idx <- which(
          annot_data$Chr == chr_i &
            annot_data$Position >= as.numeric(ss$start[i]) &
            annot_data$Position <= as.numeric(ss$end[i])
        )

        if (length(idx) == 0L) next

        state_i <- as.character(ss[[state_col]][i])
        if (state_i %in% states) {
          counts[[state_i]][idx, sample_id] <- counts[[state_i]][idx, sample_id] + 1L
        }
      }
    }
  }

  out <- matrix("neutral", nrow = nrow(annot_data), ncol = length(sample_names), dimnames = list(annot_data$Gene, sample_names))
  for (sample_id in sample_names) {
    for (gene_i in annot_data$Gene) {
      calls <- vapply(states, function(state_i) counts[[state_i]][gene_i, sample_id] >= min_agree, logical(1))
      out[gene_i, sample_id] <- collapse_state_calls(calls)
    }
  }

  out
}

to_output_df <- function(mat, first_col = "sample_id") {
  out <- data.frame(sample_id = rownames(mat), mat, check.names = FALSE)
  colnames(out)[1] <- first_col
  out
}

count_state_calls <- function(mat, state) {
  sum(grepl(paste0("(^|\\+)", state, "(\\+|$)"), as.character(mat)))
}

empty_baf_qc <- function() {
  data.frame(
    sample_id = character(),
    chr = character(),
    n_snps = integer(),
    threshold_raw = numeric(),
    threshold = numeric(),
    method = character(),
    called_shift_sites = integer(),
    shift_rate = numeric(),
    scale_pair = character(),
    stringsAsFactors = FALSE
  )
}


tryCatch({

  log_msg(paste0("Starting final Step 7 for tissue: ", tissue))
  log_msg(paste0(
    "cnv.scale=", cnv.scale,
    ", loh.scale=", loh.scale,
    ", window.length=", window.length,
    ", length.iterations=", length.iterations,
    ", baf.mode=", baf.mode
  ))

  log_msg("Loading expression matrix")
  expr_df <- read.delim(expr_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"gene_id" %in% colnames(expr_df)) {
    stop("Expression matrix must contain first column named gene_id")
  }
  rownames(expr_df) <- expr_df$gene_id
  expr_df$gene_id <- NULL
  expr_data_all <- as.matrix(expr_df)
  mode(expr_data_all) <- "numeric"

  tissue_suffix <- paste0("_", tissue)
  tissue_samples <- colnames(expr_data_all)[endsWith(colnames(expr_data_all), tissue_suffix)]
  if (length(tissue_samples) == 0L) {
    tissue_samples <- grep(paste0("_", tissue, "$"), colnames(expr_data_all), value = TRUE)
  }
  if (length(tissue_samples) == 0L) {
    stop("No expression samples found for tissue: ", tissue)
  }
  log_msg(paste0("Found expression samples: ", length(tissue_samples)))
  log_msg(paste(tissue_samples, collapse = ", "))

  expr_data <- expr_data_all[, tissue_samples, drop = FALSE]

  log_msg("Loading annotation")
  annot_raw <- read.delim(annot_file, check.names = FALSE, stringsAsFactors = FALSE)
  annot_data <- prepareAnnotationDog(annot_raw)

  common_genes <- intersect(rownames(expr_data), annot_data$Gene)
  annot_data <- annot_data[match(common_genes, annot_data$Gene), , drop = FALSE]
  expr_data <- expr_data[match(common_genes, rownames(expr_data)), , drop = FALSE]

  stopifnot(identical(rownames(expr_data), annot_data$Gene))

  log_msg(paste0("Final autosome expression genes: ", nrow(expr_data)))
  log_msg(paste0("Final samples: ", ncol(expr_data)))

  effective_baf_mode <- baf.mode
  loh_list <- list()

  if (baf.mode != "off") {
    log_msg("Loading BAF files")

    baf_problem_messages <- character()
    usable_baf_samples <- 0L

    for (sample_name in tissue_samples) {
      baf_file <- file.path(baf_dir, paste0(sample_name, ".baf.primary.tsv"))

      if (!file.exists(baf_file)) {
        baf_problem_messages <- c(baf_problem_messages, paste0("Missing BAF file for ", sample_name))
        loh_list[[sample_name]] <- data.frame(chr = character(), pos = integer(), dev = numeric())
        next
      }

      baf_data <- read.delim(baf_file, check.names = FALSE, stringsAsFactors = FALSE)
      required_baf <- c("chr", "position", "total_count", "baf")
      missing_baf <- setdiff(required_baf, colnames(baf_data))

      if (length(missing_baf) > 0L) {
        baf_problem_messages <- c(
          baf_problem_messages,
          paste0("BAF file missing columns in ", sample_name, ": ", paste(missing_baf, collapse = ", "))
        )
        loh_list[[sample_name]] <- data.frame(chr = character(), pos = integer(), dev = numeric())
        next
      }

      baf_data <- baf_data[
        is_dog_autosome(baf_data$chr) &
          as.numeric(baf_data$total_count) >= 20 &
          as.numeric(baf_data$baf) > 0.2 &
          as.numeric(baf_data$baf) < 0.8,
        ,
        drop = FALSE
      ]

      if (nrow(baf_data) == 0L) {
        baf_problem_messages <- c(baf_problem_messages, paste0("No usable BAF sites after filtering for ", sample_name))
        loh_list[[sample_name]] <- data.frame(chr = character(), pos = integer(), dev = numeric())
        next
      }

      baf_data$chr <- as.character(baf_data$chr)
      baf_data$pos <- as.integer(baf_data$position)
      baf_data$dev <- abs(as.numeric(baf_data$baf) - 0.5)

      loh_list[[sample_name]] <- baf_data[, c("chr", "pos", "dev"), drop = FALSE]
      usable_baf_samples <- usable_baf_samples + 1L
      log_msg(paste0("Loaded BAF ", sample_name, ": ", nrow(loh_list[[sample_name]]), " usable sites"))
    }

    if (length(baf_problem_messages) > 0L) {
      for (msg in baf_problem_messages) {
        log_msg(paste0("BAF sample fallback reason: ", msg))
      }
    }

    if (usable_baf_samples == 0L) {
      if (baf.mode == "auto") {
        effective_baf_mode <- "off"
        loh_list <- list()
        log_msg("BAF auto mode fallback triggered; no samples have usable BAF, switching to expression-only mode")
      } else {
        stop("No samples have usable BAF sites after filtering; use --baf-mode auto/off or fix BAF inputs")
      }
    } else if (!identical(names(loh_list), colnames(expr_data))) {
      stop("BAF sample names and expression column names are not identical")
    } else if (length(baf_problem_messages) > 0L) {
      effective_baf_mode <- paste0(baf.mode, "_sample_fallback")
      log_msg(paste0("Continuing with expression+BAF for usable samples; sample-level fallback count: ",
                     length(tissue_samples) - usable_baf_samples))
    }
  }

  cytoband_dog <- getDogCytoband(fai_file)

  log_msg(paste0("Creating CaSpER object (effective_baf_mode=", effective_baf_mode, ")"))

  if (effective_baf_mode == "off") {
    casper_obj <- CreateCasperObject(
      raw.data = expr_data,
      annotation = annot_data,
      control.sample.ids = colnames(expr_data),
      cytoband = cytoband_dog,
      loh.name.mapping = NULL,
      cnv.scale = cnv.scale,
      loh.scale = loh.scale,
      method = "iterative",
      loh = NULL,
  project = paste0("dog_", tissue, "_final_expr_only"),
      matrix.type = "normalized",
      sequencing.type = "bulk",
      expr.cutoff = 4.5,
      log.transformed = TRUE,
      window.length = window.length,
      length.iterations = length.iterations,
      genomeVersion = "canFam3",
      filter = "median"
    )
  } else {
    loh.name.mapping <- data.frame(
      loh.name = tissue_samples,
      sample.name = tissue_samples,
      stringsAsFactors = FALSE
    )

    casper_obj <- CreateCasperObject(
      raw.data = expr_data,
      annotation = annot_data,
      control.sample.ids = colnames(expr_data),
      cytoband = cytoband_dog,
      loh.name.mapping = loh.name.mapping,
      cnv.scale = cnv.scale,
      loh.scale = loh.scale,
      method = "iterative",
      loh = loh_list,
      project = paste0("dog_", tissue, "_final_hybrid"),
      matrix.type = "normalized",
      sequencing.type = "bulk",
      expr.cutoff = 4.5,
      log.transformed = TRUE,
      window.length = window.length,
      length.iterations = length.iterations,
      genomeVersion = "canFam3",
      filter = "median"
    )
  }

  casper_obj@cytoband <- cytoband_dog
  if (effective_baf_mode != "off") {
    casper_obj@loh <- loh_list
  }

  final_objects <- list()
  raw_hmm_objects <- list()
  baf_qc_all <- empty_baf_qc()

  if (effective_baf_mode == "off") {
    log_msg("Running expression-only mode from issue #22 guidance")

    for (cs in seq_len(cnv.scale)) {
      key <- paste0("cs", cs, "_exprOnly")
      log_msg(paste0("Processing ", key))

      obj_raw <- PerformSegmentationWithHMMDog(casper_obj, cnv.scale = cs)
      obj_final <- obj_raw
      obj_final@segments <- applyExpressionOnlyMode(obj_final@segments)

      raw_hmm_objects[[key]] <- generateLargeScaleEventsDog(obj_raw, state_col = "rawCNVState")
      final_objects[[key]] <- generateLargeScaleEventsDog(obj_final, state_col = "finalCNVState")
    }
  } else {
    log_msg(paste0("Running ", cnv.scale, " x ", loh.scale, " CNV/BAF scale pairs"))

    for (cs in seq_len(cnv.scale)) {
      obj_base <- PerformSegmentationWithHMMDog(casper_obj, cnv.scale = cs)

      for (ls in seq_len(loh.scale)) {
        key <- paste0("cs", cs, "_ls", ls)
        log_msg(paste0("Processing ", key))

        loh_smoothed <- smoothBafListDog(
          loh_list,
          loh.scale = ls,
          n = window.length,
          scale.iteration = length.iterations
        )

        baf_called <- callBafShiftGMMDog(
          loh_smoothed,
          min_shift_threshold = min.shift.threshold,
          high_rate_threshold = high.rate.threshold,
          high_rate_fallback_threshold = high.rate.fallback.threshold
        )

        baf_qc <- baf_called$qc
        baf_qc$scale_pair <- key
        baf_qc_all <- rbind(baf_qc_all, baf_qc)

        obj_raw <- obj_base
        obj_final <- obj_base

        for (sample_name in unique(as.character(obj_final@segments$ID))) {
          rows <- as.character(obj_final@segments$ID) == sample_name
          obj_final@segments[rows, ] <- correctSegmentsByBafDog(
          obj_final@segments[rows, , drop = FALSE],
            baf_called$loh[[sample_name]]
          )
        }

        raw_hmm_objects[[key]] <- generateLargeScaleEventsDog(obj_raw, state_col = "rawCNVState")
        final_objects[[key]] <- generateLargeScaleEventsDog(obj_final, state_col = "finalCNVState")
      }
    }
  }

  log_msg("Merging large-scale CNV calls")
  final_chr_mat <- extractLargeScaleEventsDog(final_objects, thr = 0.75)
  raw_chr_mat <- extractLargeScaleEventsDog(raw_hmm_objects, thr = 0.75)
  final_chr_state_mat <- extractLargeScaleStateMatrixDog(final_objects, state_col = "finalCNVState", thr = 0.75)
  raw_chr_state_mat <- extractLargeScaleStateMatrixDog(raw_hmm_objects, state_col = "rawCNVState", thr = 0.75)

  log_msg("Generating gene-level CNV matrices")
  gene_cnv_mat <- geneLevelCNVDog(final_objects, annot_data, state_col = "finalCNVState", thr = 0.75)
  gene_cnv_raw_mat <- geneLevelCNVDog(raw_hmm_objects, annot_data, state_col = "rawCNVState", thr = 0.75)
  gene_cnv_state_mat <- geneLevelCNVStateDog(final_objects, annot_data, state_col = "finalCNVState", thr = 0.75)
  gene_cnv_raw_state_mat <- geneLevelCNVStateDog(raw_hmm_objects, annot_data, state_col = "rawCNVState", thr = 0.75)

  log_msg("Saving outputs")

  save(
    final_objects,
    raw_hmm_objects,
    casper_obj,
    expr_data,
    annot_data,
    loh_list,
    file = file.path(out_tissue, "casper_objects.rda")
  )

  segments_all <- NULL
  for (nm in names(final_objects)) {
    tmp <- final_objects[[nm]]@segments
    tmp$scale_pair <- nm
    segments_all <- rbind(segments_all, tmp)
  }

  write.table(
    segments_all,
    file = file.path(out_tissue, "segments.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  write.table(
    to_output_df(final_chr_state_mat, "sample_id"),
    file = file.path(out_tissue, "large_scale_events.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  write.table(
    to_output_df(raw_chr_state_mat, "sample_id"),
    file = file.path(out_tissue, "large_scale_events_raw_hmm.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  write.table(
    data.frame(Gene = rownames(gene_cnv_state_mat), gene_cnv_state_mat, check.names = FALSE),
    file = file.path(out_tissue, "gene_level_cnv.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  write.table(
    data.frame(Gene = rownames(gene_cnv_raw_state_mat), gene_cnv_raw_state_mat, check.names = FALSE),
    file = file.path(out_tissue, "gene_level_cnv_raw_hmm.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  write.table(
    baf_qc_all,
    file = file.path(out_tissue, "baf_shift_qc.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  run_config <- data.frame(
    key = c(
      "tissue",
      "cnv_scale",
      "loh_scale",
      "requested_baf_mode",
      "effective_baf_mode",
      "window_length",
      "length_iterations",
      "min_shift_threshold",
      "high_rate_threshold",
      "high_rate_fallback_threshold"
    ),
    value = c(
      tissue,
      cnv.scale,
      loh.scale,
      baf.mode,
      effective_baf_mode,
      window.length,
      length.iterations,
      min.shift.threshold,
      high.rate.threshold,
      high.rate.fallback.threshold
    ),
    stringsAsFactors = FALSE
  )

  write.table(
    run_config,
    file = file.path(out_tissue, "run_config.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  log_msg("Generating heatmaps")

  png(file.path(out_tissue, "large_scale_events_heatmap.png"), width = 1600, height = 900, res = 150)
  pheatmap(
    final_chr_mat,
    color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(11),
    breaks = seq(-1, 1, length.out = 12),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    main = paste0("Dog autosome CNV final: ", tissue, " (", effective_baf_mode, ")")
  )
  dev.off()

  png(file.path(out_tissue, "large_scale_events_raw_hmm_heatmap.png"), width = 1600, height = 900, res = 150)
  pheatmap(
    raw_chr_mat,
    color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(11),
    breaks = seq(-1, 1, length.out = 12),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    main = paste0("Dog autosome CNV raw HMM: ", tissue)
  )
  dev.off()

  log_msg(paste0("Saved segments: ", nrow(segments_all)))
  log_msg(paste0(
    "Final large-scale calls: amp=", count_state_calls(final_chr_state_mat, "amp"),
    ", del=", count_state_calls(final_chr_state_mat, "del"),
    ", cnloh=", count_state_calls(final_chr_state_mat, "cnloh")
  ))
  log_msg(paste0(
    "Raw HMM large-scale calls: amp=", count_state_calls(raw_chr_state_mat, "amp"),
    ", del=", count_state_calls(raw_chr_state_mat, "del"),
    ", cnloh=", count_state_calls(raw_chr_state_mat, "cnloh")
  ))
  log_msg(paste0(
    "Final gene-level calls: amp=", count_state_calls(gene_cnv_state_mat, "amp"),
    ", del=", count_state_calls(gene_cnv_state_mat, "del"),
    ", cnloh=", count_state_calls(gene_cnv_state_mat, "cnloh")
  ))
  log_msg(paste0(
    "Raw HMM gene-level calls: amp=", count_state_calls(gene_cnv_raw_state_mat, "amp"),
    ", del=", count_state_calls(gene_cnv_raw_state_mat, "del"),
    ", cnloh=", count_state_calls(gene_cnv_raw_state_mat, "cnloh")
  ))
  log_msg("Final Step 7 completed successfully")

}, error = function(e) {
  log_msg(paste0("ERROR: ", conditionMessage(e)))
  quit(save = "no", status = 1)
})
