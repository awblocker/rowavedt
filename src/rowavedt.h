/*
 * waveletLMt.h
 *
 *  Created on: Jun 17, 2010
 *      Author: awblocker
 */

#ifndef WAVELETLMT_H_
#define WAVELETLMT_H_

#define GSL_C99_INLINE

// C libraries
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <unistd.h>

// GSL
#include <gsl_cdf.h>
#include <gsl_math.h>
#include <gsl_statistics_double.h>

// BLAS-LAPACK interface (direct to Fortran)
#include "interfaceBLAS-LAPACK.h"

// utils.c
void checkPtr(const void * ptr, const char * msg);
int arrayMinMax(double * X, int n, double * min, double * max);
int allocateMatrix(double*** matrix, int nRows, int nCols);
void freeMatrix(double*** matrix, int nRows, int nCols);
double arraymean(double x[], int start, int stop);
double arrayvar(double x[], int start, int stop);
double bisect(double x, double data[], int start, int end);
double quantile(double x, double data[], int n);
double quantile_int(double x, double data[], int n);
int compare_dbl(const void * a, const void * b);
void readToDoubleMatrix (const char * fname, int nRows, int nCols, double *X);
int readToDoubleVector (const char * fname, int nRows, int col, double* X);
int readToDoubleVectorDynamic(const char * fname, int startRows, int col,
    double** X);

// wls.c
int wls(double* X, int n, int k,
    double* y, double* w,
    double* XTX, double *sqw, double* sqwX, double* sqwy,
    double* Xy);
int calcFitted(double* X, int n, int k,
    double* y,
    double* coef,
    double* fitted);
int calcResid(double* X, int n, int k,
    double* y,
    double* coef,
    double* resid);

// dist.c
double dnorm_log(double * x, int n, double location, double scale);
double dt_log(double * x, int n, double df, double location, double scale);

// lmT.c
int lmT(double * basisMat, int basisRows, int basisCols,
    double * valueVec, int n,
    double * timeVec,
    double * priorVec,
    double nu, int k,
    int maxIter, double tol,
    double * logPosterior,
    double * logLikelihood,
    double * coef,
    double * tau);

#endif /* WAVELETLMT_H_ */
