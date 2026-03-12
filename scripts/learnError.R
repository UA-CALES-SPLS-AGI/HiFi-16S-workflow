library(dada2)

# Read from command line arguments
args <- commandArgs(trailingOnly = TRUE)
fnFs <- args[1]
cpu <- args[2]
error_model_mode <- ifelse(length(args) >= 3, args[3], "auto")

# Cap CPU to 16. Higher doesn't really help
if (cpu > 16) {
  cpu <- 16
}

### Error model selection and auto-detection ###

# Detect quality score distribution from FASTQ file(s)
detect_qscore_bins <- function(fastq_files, n_reads = 5000) {
  fq <- fastq_files[1]
  if (grepl("\\.gz$", fq)) {
    con <- gzfile(fq, "r")
  } else {
    con <- file(fq, "r")
  }
  on.exit(close(con))
  lines <- readLines(con, n = n_reads * 4)
  qual_lines <- lines[seq(4, length(lines), 4)]
  all_quals <- unlist(lapply(qual_lines, function(q) {
    as.integer(charToRaw(q)) - 33L
  }))
  sort(unique(all_quals))
}

# Polyfill for makeBinnedQualErrfun when dada2 < 1.32
binned_err_polyfill <- function(bin_quals) {
  function(trans) {
    qq <- as.numeric(colnames(trans))
    est <- matrix(0, nrow = 0, ncol = length(qq))
    for (nti in c("A", "C", "G", "T")) {
      for (ntj in c("A", "C", "G", "T")) {
        if (nti != ntj) {
          errs <- trans[paste0(nti, "2", ntj), ]
          tot <- colSums(trans[paste0(nti, "2", c("A", "C", "G", "T")), ])
          rlogp <- log10((errs + 1) / tot)
          rlogp[is.infinite(rlogp)] <- NA
          bin_idx <- which(qq %in% bin_quals)
          if (length(bin_idx) >= 2) {
            pred <- approx(qq[bin_idx], rlogp[bin_idx], xout = qq, rule = 2)$y
          } else {
            pred <- rlogp
          }
          pred[is.na(pred)] <- min(rlogp, na.rm = TRUE)
          est <- rbind(est, 10^pred)
        }
      }
    }
    MAX_ERROR_RATE <- 0.25
    MIN_ERROR_RATE <- 1e-7
    est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE
    est[est < MIN_ERROR_RATE] <- MIN_ERROR_RATE
    err <- rbind(1 - colSums(est[1:3, ]), est[1:3, ],
                 est[4, ], 1 - colSums(est[4:6, ]), est[5:6, ],
                 est[7:8, ], 1 - colSums(est[7:9, ]), est[9, ],
                 est[10:12, ], 1 - colSums(est[10:12, ]))
    rownames(err) <- paste0(rep(c("A", "C", "G", "T"), each = 4), "2", c("A", "C", "G", "T"))
    colnames(err) <- colnames(trans)
    return(err)
  }
}

# Get the binned error function, using native dada2 if available, polyfill otherwise
get_binned_errfun <- function(scores) {
  if (exists("makeBinnedQualErrfun", where = asNamespace("dada2"), inherits = FALSE)) {
    cat("Using dada2::makeBinnedQualErrfun (native)\n")
    return(dada2::makeBinnedQualErrfun(scores))
  } else {
    cat("dada2::makeBinnedQualErrfun not available (requires dada2 >= 1.32).",
        "Using polyfill with piecewise linear interpolation.\n")
    return(binned_err_polyfill(scores))
  }
}

# Determine the appropriate error estimation function
get_error_function <- function(mode, fastq_files) {
  if (mode == "pacbio") {
    cat("Error model: PacBioErrfun (legacy PacBio CCS with continuous Q-scores)\n")
    return(dada2:::PacBioErrfun)
  } else if (mode == "loess") {
    cat("Error model: loessErrfun (standard Illumina)\n")
    return(dada2:::loessErrfun)
  }

  scores <- detect_qscore_bins(fastq_files)
  cat("Detected", length(scores), "unique quality scores, range Q",
      min(scores), "-Q", max(scores), "\n")

  if (mode == "binned") {
    cat("Error model: Binned with scores:", paste(scores, collapse = ", "), "\n")
    return(get_binned_errfun(scores))
  }

  # Auto mode
  if (max(scores) >= 93) {
    cat("Auto: Q93 present — using PacBioErrfun (Sequel/Sequel II continuous Q-scores)\n")
    return(dada2:::PacBioErrfun)
  } else if (length(scores) <= 15) {
    cat("Auto: Binned quality scores detected (", length(scores), " unique values, max Q",
        max(scores), ") — using makeBinnedQualErrfun\n")
    cat("  Bin centers:", paste(scores, collapse = ", "), "\n")
    return(get_binned_errfun(scores))
  } else {
    cat("Auto: Continuous quality scores detected — using loessErrfun\n")
    return(dada2:::loessErrfun)
  }
}

# Get the error estimation function
errfun <- get_error_function(error_model_mode, fnFs)

# Learn the error rates
errF <- learnErrors(fnFs, multithread = cpu, errorEstimationFunction = errfun)
err_plot <- plotErrors(errF)
pdf("plot_error_model.pdf", width = 12, height = 8, useDingbats = FALSE)
print(err_plot)
dev.off()

# Save as RDS
saveRDS(errF, file = "errorfun.rds")
