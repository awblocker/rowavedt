rowavedt : Robust (Semi-parametric) Wavelet-based Event Detection
================================================================================ 

This package provides the core routines from Blocker and Protopapas 2012
([arXiv:1301.3027 [stat.AP]](http://arxiv.org/abs/1301.3027)). It includes the
PX-EM algorithm for estimation of the semi-parametric model presented in that
paper, as well as R scripts for the construction of the wavelet basis,
frequency-based screening for time series with interesting variation, and
construction of features for subsequent classification.

Code for the training and use of the classifier is not included. The former used
existing techniques from the [arm](http://cran.r-project.org/web/packages/arm/)
package for R in combination with scientific simulations (which are not easily
generalized to other uses). The latter is straightforward, once the features are
constructed and the classifier is trained.

Requirements
--------------------------------------------------------------------------------

This package is written for *nix systems, although it may work on Mac as well
(untested). The core estimation code (written in C) requires a BLAS library with
[LAPACK](http://www.netlib.org/lapack/) configured and the [GNU Science
Library](https://www.gnu.org/software/gsl/). The
[ATLAS](http://math-atlas.sourceforge.net/) and
[Intel MKL](http://software.intel.com/en-us/intel-mkl/) BLAS libraries have been
tested successfully, and compilation options for both are included in the
Makefile.

All R scripts require the
[optparse](http://cran.r-project.org/web/packages/optparse/) package. The script
to construct wavelet basis matrices (`mk_wavelet_basis.R`) requires the
[wavethresh](http://cran.r-project.org/web/packages/wavethresh/) package as
well.


Installation
--------------------------------------------------------------------------------

To compile and install the core estimation routine (written in C), first edit
the included Makefile as needed.


