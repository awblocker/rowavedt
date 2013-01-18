/*
 * dist.c
 *
 *  Created on: Jun 22, 2010
 *      Author: awblocker
 */

#include "rowavedt.h"

// Compute sum of log-density for normal
// Omits normalizing constants
double dnorm_log(double * x, int n, double location, double scale)
{
	int i;
	double logDensity = 0, z, logScale;

	logScale = log(scale);

	for (i=0; i<n; i++)
	{
		z = (x[i] - location) / scale;
		logDensity += -0.5*z*z - logScale;
	}

	return logDensity;
}

// Compute sum of log-density for t distribution
// Omits normalizing constants
double dt_log(double * x, int n, double df, double location, double scale)
{
	int i;
	double logDensity = 0, z, logScale;

	logScale = log(scale);

	for (i=0; i<n; i++)
	{
		z = (x[i] - location) / scale;
		logDensity += -(df+1)/2 * log(1 + z*z/df) - logScale;
	}

	return logDensity;
}
