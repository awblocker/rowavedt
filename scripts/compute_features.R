#!/usr/bin/env Rscript

# Load libraries
if(!suppressMessages(require(optparse, quietly=TRUE))) {
  cat("Error - optparse library is required.\n", file=stderr())
  q(status=1)
}


# Set constants

kInputSep <- ' '
kCoefColOffset <- 8
kInputFmt <- list(id=character(0),
                  n=integer(0),
                  p=integer(0),
                  p.low=integer(0),
                  nu=numeric(0),
                  llr=numeric(0),
                  lpr=numeric(0),
                  tau=numeric(0)) 

kOptList <- list(
  make_option(c('--detections'), type='character', default='',
              help='Optional detections file to filter input.'),
  make_option(c('--sep', '-s'), default=' ',
              help=paste('Separator for output.',
                         '\n\t\tDefaults to %default.', sep=''))
)

kUsage <- '
Usage: Rscript compute_features.R [options] INPUT BASIS [> output.txt]'
kEpilogue <- '
Constructs features defined in Blocker and Protopapas (2012) on output
from rowavedt. These are engineered for the separation of distinct,
temporally-isolated events from other types of variability. BASIS must
correspond to the full basis matrix used when running rowavedt.

Optionally filters INPUT using detections file as output by
screen_time_series.R. This requires that both INPUT and DETECTIONS are sorted by
their first columns (ID) before running this script. This allows for the join
command to be used for filtering, greatly reducing the amount of data that must
be read into memory.

Outputs a list a SEP-delimited list to stdout with 3 columns:
  1 - ID from INPUT
  2 - CUSUM feature as defined in Blocker and Protopapas (2012)
  3 - Directed variation (DV) feature as defined in Blocker and Protopapas
      (2012)

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

#' Compute features from coefficients and basis
#'
#' @param b Vector of estimated coefficients. Can be complete or subsetted.
#' @param X Basis matrix with ncol(X)==length(b). Can be complete or subsetted.
#' @param subset.cols Indices by which to subset b and columns of X. Useful in
#'  restricting features to time scales of interest.
#'
#' @returns A named numeric vector of length 2 containing the cusum and dv
#'  features as defined in Blocker and Protopapas 2012.
#'
ComputeFeatures <- function(b, X, subset.cols) {
  # Compute standardized fitted time series at time scales of interest
  y.hat <- X[,subset.cols] %*% b[subset.cols]
  z <- scale(y.hat)
  n <- length(z)

  # Compute CUSUM feature
  S <- cumsum(z^2 - 1)
  cusum <- log(1 + (max(S) - min(S)) / sqrt(n))

  # Compute directed variation feature (DV)
  z.med <- median(z)
  dv <- mean(z[z >= z.med]^2) - mean(z[z < z.med]^2)

  return(c(cusum=cusum, dv=dv))
}

# Parse options
optargs <- parse_args(kParser, commandArgs(TRUE), positional_arguments=TRUE) 
opts <- optargs$options
args <- optargs$args

if (length(args) < 2) {
  cat('Error - need INPUT and BASIS.\n',
      file=stderr())
  q(status=1)
}

input.path <- args[1]
basis.path <- args[2]

# Load basis
basis.conn <- file(basis.path)
basis.line.1 <- readLines(basis.conn, n=1)
basis.n.col <- length(strsplit(basis.line.1, '\\s', perl=TRUE)[[1]])
basis <- matrix(scan(basis.conn, 0.), ncol=basis.n.col, byrow=TRUE)
close(basis.conn)

# Load input
coef.cols <- kCoefColOffset + 1:basis.n.col
input.fmt <- c(kInputFmt, rep(0., basis.n.col))

if (opts$detections != '') {
  # Filter input using detections; simple command-line method
  join.cmd <- 'join %s %s'
  
  # Open connection and read data
  input.fmt <- c(input.fmt, list(p=0., q=0.))
  conn <- pipe(sprintf(join.cmd, input.path, opts$detections))
  input <- scan(conn, input.fmt, sep=kInputSep, flush=TRUE)
  close(conn)
} else {
  # Load all of INPUT
  input <- scan(input.path, input.fmt, sep=kInputSep)
}

input <- data.frame(input)

# Compute features on each observation
n.obs <- nrow(input)
features <- sapply(1:n.obs, function(ii, input, coef.cols)
                   ComputeFeatures(b=as.numeric(input[ii, coef.cols]),
                                   X=basis,
                                   subset.cols=(input$p.low[ii]+1):basis.n.col),
                   input=input)
features <- t(features)

# Build output
output <- data.frame(id=input$id, features)

# Write output to stdout
write.table(output, file=stdout(), row.names=FALSE, col.names=FALSE,
            quote=FALSE)

