/*
 * wls.c
 *
 *  Created on: Jun 18, 2010
 *      Author: awblocker
 */

#include "rowavedt.h"

int wls(double* X, int n, int k,
		double* y, double* w,
		double* XTX, double *sqw, double* sqwX, double* sqwy,
		double* coef) {
	// Assuming column-major order
	// Initializations
	int i, j;

	// Compute sqw
	for (i=0; i<n; i++)
	{
		sqw[i] = sqrt(w[i]);

		// Compute sqwX and sqwy
		for (j=0; j<k; j++) {
			sqwX[i + j*n] = X[i + j*n] * sqw[i];
		}
		sqwy[i] = sqw[i] * y[i];
	}

	// Compute XTX
	dsyrk('u', 't', k, n, 1, sqwX, n, 0, XTX, k);

	// Compute Xy
	dgemv('t', n, k, 1, sqwX, n, sqwy, 1, 0, coef, 1);

	// Obtain least-squares coefficients by solving normal equations
	int info;
	info = dposv('u', k, 1, XTX, k, coef, k);

//	// Obtain least-squares coefficients via LAPACK DGELS driver routine
//	int info;
//	info = dgels('n', n, k, 1, sqwX, n, sqwy, n);
//
//	// Copy solutions to coef vector
//	dcopy(k, sqwy, 1, coef, 1);

	return info;
}

int calcFitted(double* X, int n, int k,
		double* y,
		double* coef,
		double* fitted)
{
	// X \beta = \hat{X}
	dgemv('n', n, k, 1, X, n, coef, 1, 0, fitted, 1);

	return 0;
}

int calcResid(double* X, int n, int k,
		double* y,
		double* coef,
		double* resid)
{
	// Calculate fitted values; store temporarily in resid
	dgemv('n', n, k, 1, X, n, coef, 1, 0, resid, 1);

	// Calculate residuals
	dscal(n, -1, resid, 1);
	daxpy(n, 1, y, 1, resid, 1);

	return 0;
}
