#!/usr/bin/env Rscript

# Load libraries
if (!suppressMessages(require(wavethresh, quietly=TRUE))) {
  cat("Error - wavethresh library is required.\n", file=stderr())
  q(status=1)
}

if(!suppressMessages(require(getopt, quietly=TRUE))) {
  cat("Error - getopt library is required.\n", file=stderr())
  q(status=1)
}


# Set constants

kDefaults <- list(filter.family='DaubLeAsymm',
                  filter.number=4,
                  bc='symmetric',
                  B=2048,
                  sep=' ',
                  help=FALSE)

kOptSpec <- matrix(c(
  'filter.family', 'f', 1, 'character',
  'filter.number', 'i', 1, 'integer',
  'bc', 'b', 1, 'character',
  'B', 'B', 1, 'integer',
  'n', 'n', 1, 'integer',
  'sep', 's', 1, 'character',
  'help', 'h', 0, 'logical')
  ncol=5, byrow=TRUE)

kHelp <- '
Usage: Rscript mk_wavelet_basis.R [options] [> basis.txt]
  -f|--filter.family    A filter family recognized by wavethresh.
                        Defaults to DaubLeAsymm.
  -i|--filter.number    Smoothness of wavelet to use for basis.
                        Defaults to 4.
  -b|--bc               Boundary handling.
                        Defaults to symmetric.
  -B                    Size of basis to construct. Must be a power of 2.
                        Defaults to 2048.
  -n                    Number of columns of basis to output (starting from
                        lowest frequencies). Defaults to B.
  -s|--sep              Separator for output. Defaults to " ".
  -h|--help             Display this help message and exit.

Constructs basis matrix corresponding to given wavelet filter, outputting the
first n columns to stdout as SEP-separated text.

Requires the wavethresh and getopt R packages (available on CRAN).

'


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
opts <- getopt(kOptSpec)

# Set defaults, if needed
for (oo in names(kDefaults)) {
  if (!(oo %in% names(opts)))
    opts[[oo]] <- kDefaults[[oo]]
}

# Handle n separately
if (!('n' %in% names(opts)))
  opts$n <- opts$B

# Display help message and exit, if requested
if (opts$help) {
  cat(kHelp)
  q(status=1)
}


# Build zero sequence and base wt
y <- rep(0, opts$B)
wt <- wd(y, filter.number=opts$filter.number, family=opts$filter.family,
         bc=opts$bc)

# Reconstruct wavelet basis
X <- basis.reconstruct(wt)


# Write basis to stdout
write.table(X[, 1:opts$n], file=stdout(), sep=opts$sep,
            row.names=FALSE, col.names=FALSE)

