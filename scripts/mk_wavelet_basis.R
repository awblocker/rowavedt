#!/usr/bin/env Rscript

# Load libraries
if (!suppressMessages(require(wavethresh, quietly=TRUE))) {
  cat("Error - wavethresh library is required.\n", file=stderr())
  q(status=1)
}

if(!suppressMessages(require(optparse, quietly=TRUE))) {
  cat("Error - optparse library is required.\n", file=stderr())
  q(status=1)
}


# Set constants

kOptList <- list(
  make_option(c('--filter.family', '-f'), default='DaubLeAsymm',
              help=paste('A filter family recognized by wavethresh.',
                         '\n\t\tDefaults to %default.', sep='')),
  make_option(c('--filter.number', '-i'), type='integer', default=4,
              help=paste('Smoothness of wavelet to use for basis.',
                         '\n\t\tDefaults to %default.', sep='')),
  make_option(c('--boundary', '-b'), default='periodic',
              help=paste('Boundary handling setting.',
                         '\n\t\tDefaults to %default.', sep='')),
  make_option(c('--B', '-B'), type='integer', default=2048,
              help=paste('Size of basis to construct.',
                         '\n\t\tDefaults to %default.', sep='')),
  make_option(c('--c', '-c'), type='integer', default=0,
              help=paste('Number of columns of basis to output (starting from ',
                         'lowest\n\t\tfrequencies). Defaults to B.', sep='')),
  make_option(c('--sep', '-s'), default=' ',
              help=paste('Separator for output.',
                         '\n\t\tDefaults to %default.', sep=''))
)

kUsage <- 'Usage: Rscript mk_wavelet_basis.R [options] [> basis.txt]'
kEpilogue <- '
Constructs basis matrix corresponding to given wavelet filter, outputting the
first n columns to stdout as SEP-separated text.
Requires the wavethresh and optparse R packages (available on CRAN).

'

kParser <- OptionParser(usage=kUsage, option_list=kOptList,
                        epilogue=kEpilogue)


# Function definitions

# Single basis vector reconstruction
basis.vector <- function(k, level, wt) {
    # Make basis vector in wavelet space
    v <- rep(0, 2^level)
    v[k] <- 1
    
    # Reconstruct and return basis vector
    tmp <- putD(wt, level, v)
    return( wr(tmp) )
}

# Single level basis reconstruction
basis.level <- function(level, wt) {
    # Run basis.vector for given level across coefficients
    sapply( seq(1, 2^level), basis.vector,
        level=level, wt=wt )
}

# Full basis reconstruction
basis.reconstruct <- function(wt) {    
    # Extract size information
    nlevels <- wt$nlevels
    n <- 2^nlevels
    
    # Run basis.level across levels
    X <- sapply( seq(0, wt$nlevels-1), basis.level,
        wt=wt )
    # Make into a matrix
    X <- matrix(unlist(X), nrow=n)
    
    # Add intercept column
    X <- cbind( rep( 1/sqrt(n), n), X)
}


# Script

# Parse options
optargs <- parse_args(kParser, commandArgs(TRUE), positional_arguments=TRUE) 
opts <- optargs$options

# Handle c separately
if (opts$c < 1)
  opts$c <- opts$B


# Build zero sequence and base wt
y <- rep(0, opts$B)
wt <- wd(y, filter.number=opts$filter.number, family=opts$filter.family,
         bc=opts$boundary)

# Reconstruct wavelet basis
X <- basis.reconstruct(wt)


# Write basis to stdout
write.table(X[, 1:opts$c], file=stdout(), sep=opts$sep,
            row.names=FALSE, col.names=FALSE)

