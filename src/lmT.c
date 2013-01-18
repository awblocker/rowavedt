
#include "rowavedt.h"

/*
 * Function to run wavelet model for irregularly sampled data
 * Returns number of iterations run
 * Log-posterior, log-likelihood, and coefficients are returned by reference
 */

int lmT(double * basisMat, int basisRows, int basisCols,
		double * yVec, int n,
		double * timeVec,
		double * priorVec,
		double nu, int k,
		int maxIter, double tol,
		double * logPosterior,
		double * logLikelihood,
		double * coef,
		double * tau)
{
	// Initialize workspace variables
	int iter, i, j, m, tme;
	double logPosterior_tm1, logPrior, delta;
	m = n + k - 1;

	// Rescale times
	double minTime, maxTime;
	arrayMinMax(timeVec, n, &minTime, &maxTime);

	for (i=0; i<n; i++)
	{
		timeVec[i] = (timeVec[i] - minTime) / (maxTime - minTime);
		// timeVec[i] = timeVec[i] * (n-1);
	    timeVec[i] = timeVec[i] * (basisRows-1);
	}

	/*
	 * Setup dmat and dvec
	 */

	double * dMat, * dVec;
	dMat = calloc( m * k, sizeof(double));
	if (dMat==NULL)
	{
		fprintf(stderr, "Error -- out of memory\n");
		exit(1);
	}

	dVec = calloc( m, sizeof(double) );
	if (dVec==NULL)
	{
		fprintf(stderr, "Error -- out of memory\n");
		exit(1);
	}

	// Copy observations to dVec; leave remaining entries 0 for prior
	dcopy(n, yVec, 1, dVec, 1);

	// Copy correct rows from basisMat to dMat
	for (i=0; i<n; i++)
	{
		tme = (int) floor( timeVec[i] + 0.5 );
		for (j=0; j<k; j++)
		{
			dMat[i + j*m] = basisMat[tme + j*basisRows];
		}
	}

	// Initialize remaining rows of dMat to identity (w/o first row) for prior
	for (i=n; i<m; i++)
	{
		j = i - n + 1;
		dMat[i + j*m] = 1;
	}

	/*
	 * Allocate workspace arrays
	 */

	// Allocate workspace matrices
	double * sqwX, * XTX;
	sqwX = malloc(m * k * sizeof(double));
	if (sqwX==NULL)
	{
		fprintf(stderr, "Error -- out of memory\n");
		exit(1);
	}

	XTX = malloc(k * k * sizeof(double));
	if (XTX==NULL)
	{
		fprintf(stderr, "Error -- out of memory\n");
		exit(1);
	}

	// Allocate workspace vectors
	double * w, * sqw, * sqwy, * resid;

	w = calloc(m, sizeof(double));
	if (w==NULL)
	{
		fprintf(stderr, "Error -- out of memory\n");
		exit(1);
	}

	sqw = calloc(m, sizeof(double));
	if (sqw==NULL)
	{
		fprintf(stderr, "Error -- out of memory\n");
		exit(1);
	}

	sqwy = calloc(m, sizeof(double));
	if (sqwy==NULL)
	{
		fprintf(stderr, "Error -- out of memory\n");
		exit(1);
	}

	resid = calloc(m, sizeof(double));
	if (resid==NULL)
	{
		fprintf(stderr, "Error -- out of memory\n");
		exit(1);
	}

	/*
	 * Initialize quantities before EM iterations
	 */

	// Initialize weights for observations
	for (i=0; i<n; i++)
	{
		w[i] = 1;
	}

	// Copy prior information to weights
	dcopy(k-1, priorVec, 1, &w[n], 1);

	/*
	 * First iteration
	 */

	// Run regression to obtain initial coefficients
	wls(dMat, m, k, dVec, w, XTX, sqw, sqwX, sqwy, coef);

	// Calculate residuals
	calcResid(dMat, m, k, dVec, coef, resid);

	// Calculate tau
	(*tau) = 0;
	for (i=0; i<m; i++)
	{
		(*tau) += resid[i] * resid[i] * w[i];
	}
	(*tau) /= dasum(m, w, 1);

	// Calculate initial log-posterior
	(*logLikelihood) = dt_log(resid, n, nu, 0, sqrt(*tau));
	logPrior = dnorm_log(&resid[n], k-1, 0, sqrt((*tau)/priorVec[0]));
	(*logPosterior) = (*logLikelihood) + logPrior;

	logPosterior_tm1 = (*logPosterior);

	/*
	 * EM loop
	 */
	for (iter=0; iter<maxIter; iter++)
	{
		// E step: Update u | beta, tau
		for (i=0; i<n; i++)
		{
			w[i] = (nu + 1) / (nu + resid[i]*resid[i]/(*tau));
		}

		// M step: Update beta, tau | u

		// Run regression to obtain coefficients
		wls(dMat, m, k, dVec, w, XTX, sqw, sqwX, sqwy, coef);

		// Calculate residuals
		calcResid(dMat, m, k, dVec, coef, resid);

		// Calculate tau
		(*tau) = 0;
		for (i=0; i<m; i++)
		{
			(*tau) += resid[i] * resid[i] * w[i];
		}
		(*tau) /= dasum(m, w, 1);

		// Calculate log-posterior
		(*logLikelihood) = dt_log(resid, n, nu, 0, sqrt(*tau));
		logPrior = dnorm_log(&resid[n], k-1, 0, sqrt((*tau)/priorVec[0]));
		(*logPosterior) = (*logLikelihood) + logPrior;

		// Check convergence
		delta = ((*logPosterior) - logPosterior_tm1) /
				fabs((*logPosterior) + logPosterior_tm1) * 2;

		if (delta < tol)
		{
			logPosterior_tm1 = (*logPosterior);
			break;
		}
		logPosterior_tm1 = (*logPosterior);
	}

	/*
	 * Free allocated memory
	 */

	// Free dMat and dVec
	free(dMat);
	dMat = NULL;

	free(dVec);
	dVec = NULL;

	// Free workspace matrices
	free(sqwX);
	sqwX = NULL;

	free(XTX);
	XTX = NULL;

	// Free workspace vectors
	free(w);
	w=NULL;

	free(sqw);
	sqw=NULL;

	free(sqwy);
	sqw=NULL;

	free(resid);
	resid=NULL;

	return iter;
}
