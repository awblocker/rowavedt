#!/usr/bin/env Rscript

# Load libraries
if(!suppressMessages(require(optparse, quietly=TRUE))) {
  cat("Error - optparse library is required.\n", file=stderr())
  q(status=1)
}


# Set constants
kDefaults <- list(alpha=0.001,
                  method='BH',
                  sep=' ',
                  help=FALSE)

kInputSep <- ' '
kColumns <- list(id=1,
                 dim.full=3,
                 dim.low=4,
                 llr=6)

kOptList <- list(
  make_option(c('--alpha', '-a'), type='numeric', default=0.001,
              help=paste('False discovery rate to use for detection.',
                         '\n\t\tDefaults to %default.', sep='')),
  make_option(c('--method', '-m'), default='BH',
              help=paste('Method to use for FDR-controlling detection.',
                         '\nMust be in p.adjust.methods.',
                         '\n\t\tDefaults to %default.', sep='')),
  make_option(c('--sep', '-s'), default=' ',
              help=paste('Separator for output.',
                         '\n\t\tDefaults to %default.', sep=''))
)

kUsage <- '
Usage: Rscript screen_time_series.R [options] INPUT DETECTIONS_PATH STATS_PATH'
kEpilogue <- '
Runs screening procedure described in Blocker and Protopapas (2012) on output
from rowavedt. Outputs a list of detected time series to DETECTIONS_PATH and a
set of detection statistics to STATS_PATH.

The list of detections is SEP-separated and has 3 columns:
  1 - IDs of detected time series
  3 - P-values for detections
  2 - Approximate q-values for detections

The statistics file is SEP-separated and contains two columns:
  1 - Labels of statistics
  2 - Values of statistics

Requires optparse R package.

'

kParser <- OptionParser(usage=kUsage, option_list=kOptList,
                        epilogue=kEpilogue)


# Function definitions

#' Read single column via pipe
#'
#' Uses awk for very, very fast column extract (compared to scan).
#'
#' @param path Path to data file
#' @param column Index of column to extract (starting from 1)
#' @param ... Additional arguments to scan
#'
#' @returns A vector containing the values of the requested column
#'
ReadColumnViaPipe <- function(path, column, ...) {
  kPipePattern <- "awk '{print $%d}' %s" # Pattern for awk field extraction
  conn <- pipe(sprintf(kPipePattern, column, path))
  
  dat <- NULL
  tryCatch(dat <- scan(conn, ...), finally=close(conn))
  return(dat)
}


# Script


# Parse options
optargs <- parse_args(kParser, commandArgs(TRUE), positional_arguments=TRUE) 
opts <- optargs$options
args <- optargs$args

if (length(args) < 3) {
  cat('Error - need INPUT_PATH, DETECTIONS_PATH, and STATS_PATH.\n',
      file=stderr())
  q(status=1)
}

input.path <- args[1]
detections.path <- args[2]
stats.path <- args[3]


# Load columns from input
id.vec <- ReadColumnViaPipe(input.path, kColumns$id, what='')
dim.full.vec <- ReadColumnViaPipe(input.path, kColumns$dim.full,
                                  what=integer(0))
dim.low.vec <- ReadColumnViaPipe(input.path, kColumns$dim.low,
                                 what=integer(0))
llr.vec <- ReadColumnViaPipe(input.path, kColumns$llr, what=numeric(0))

# Compute p-values using chisq approximation
p.values <- pchisq(llr.vec, df=dim.full.vec - dim.low.vec, lower.tail=FALSE)

# Adjust p-values to approximate q-values
q.values.approx <- p.adjust(p.values, method=opts$method)

# Build data.frame of screened IDs and q-values
detections <- data.frame(id=id.vec[q.values.approx <= opts$alpha],
                         p=p.values[q.values.approx <= opts$alpha],
                         q=q.values.approx[q.values.approx <= opts$alpha])

# Write detections to DETECTIONS_PATH
write.table(detections, file=detections.path, quote=FALSE, sep=opts$sep,
            row.names=FALSE, col.names=FALSE)

# Compute detection statistics
detection.stats <- data.frame(n.detected=nrow(detections),
                              pct.detected=nrow(detections) / length(id.vec),
                              avg.q.value=mean(q.values.approx[q.values.approx <
                                               opts$alpha]))

# Write detection statistics to STATS_PATH
write.table(t(detection.stats), file=stats.path, quote=FALSE, sep=opts$sep,
            row.names=TRUE, col.names=FALSE)

