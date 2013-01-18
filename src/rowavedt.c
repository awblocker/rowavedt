
#include "rowavedt.h"

const char * kHelpMessage = "\nUsage:\trowavedt [options] "
  "BASISFILE BASISROWS BASISCOLS\n\tDATAFILE DATAROWS PRIORFILE\n\n"
  "Options:\n"
  "-c\tSet column number for values. Note: This is base 0.\n"
  "\tDefaults to 1.\n"
  "-d\tSet df for t distribution of residuals.\n"
  "\tDefaults to 5.\n"
  "-i\tOptional ID for results. Used as first entry of output.\n"
  "\tDefaults to DATAFILE.\n"
  "-m\tNumeric code for missing values.\n"
  "\tDefaults to 99.999\n"
  "-n\tMinimum number of observations required; else exit\n"
  "\tDefaults to 10\n"
  "-s\tSet dimension of smooth (low-resolution) partial basis.\n"
  "\tMust be power of 2\n"
  "\tDefaults to 8.\n"
  "-t\tSet column number for times. Note: This is base 0.\n"
  "\tDefaults to 0.\n\n"
  "Outputs a single space-delimited line to stdout consisting of\n"
  "7 + BASISCOLS entries:\n"
  " 1  - ID (defaults to DATAFILE)\n"
  " 2  - Number of non-missing observations read from DATAFILE\n"
  " 3  - Number of components in full basis (BASISCOLS)\n"
  " 4  - Number of components in smooth (low-resolution) partial basis\n"
  " 5  - Degrees of freedom set for t distribution of residuals\n"
  " 6  - 2 * difference of maximized log-likelihoods between full and\n"
  "      low-frequency models\n"
  " 7  - Estimated scale for residual variance of full model\n"
  " 8: - Maximum a posteriori estimates of coefficients for full model\n"
  "       (BASISCOLS coefficients in total).\n"
  "\n";

int main(int argc, char * argv[]) {
  // Define constants
  const int nArgs = 6;

  // Process arguments
  int c, timeCol=0, valueCol=1;
  extern char *optarg;
  char * dataFile, * basisFile, * priorFile, * idString=NULL;
  short readID = 0;
  int kSmooth = 8;
  int maxIter = 1e3;
  double nu=5;
  double tol=1e-9;
  double missingCode = 99.999;
  int minObs = 10;
  int basisRows, basisCols, dataRows;
  int i;

  // Parse options
  while ( (c=getopt(argc, argv, "c:d:i:m:n:s:t:h")) != -1 ) {
    switch(c) {
      case 'h':
        puts(kHelpMessage);
        return 0;
      case 'c':
        valueCol = atoi(optarg);
        valueCol = (valueCol < 0) ? 0 : valueCol;
        break;
      case 'd':
        nu = atof(optarg);
        nu = (nu < 1) ? 1 : nu;
        break;
      case 'i':
        idString = optarg;
        readID = 1;
        break;
      case 'm':
        missingCode = atof(optarg);
        break;
      case 'n':
        minObs = atoi(optarg);
        break;
      case 's':
        kSmooth = atoi(optarg);
        kSmooth = (trunc(log2(kSmooth))-log2(kSmooth) > 1e-16) ?
          (int) gsl_pow_int(2, trunc(log2(kSmooth))) : kSmooth;
        break;
      case 't':
        timeCol = atoi(optarg);
        timeCol = (timeCol < 0) ? 1 : timeCol;
        break;
      case '?':
        if ( isprint(optopt) )
          fprintf(stderr, "Unknown argument '%c'\n", optopt);
        else
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        exit(1);
        break;
      default:
        abort();
    }
  }

  // Parse positional arguments
  if (argc-optind < nArgs) {
    fprintf(stderr, "Not enough arguments\n");
    exit(1);
  }

  // Get file names and dimensions
  basisFile = argv[optind];
  basisRows = atoi( argv[optind+1] );
  basisCols = atoi( argv[optind+2] );

  dataFile = argv[optind+3];
  dataRows = atoi( argv[optind+4] );

  priorFile = argv[optind + 5];

  if (readID == 0)
  {
    idString = dataFile;
  }

  // Read basis to allocated matrix
  double* basisMat;
  basisMat = malloc(basisCols * basisRows * sizeof(double));
  checkPtr(basisMat, "out of memory");

  readToDoubleMatrix(basisFile, basisRows, basisCols, basisMat);

  // Read times to vector
  double* timeVec;
  timeVec = malloc(dataRows * sizeof(double));
  if (timeVec==NULL)
  {
    fprintf(stderr, "Error -- out of memory\n");
    exit(1);
  }
  int timesRead = readToDoubleVectorDynamic(dataFile, dataRows,
      timeCol, &timeVec);

  // Check for valid first time
  if (isnan(timeVec[0]))
  {
    fprintf(stderr, "Error -- NaN at first time; aborting\n");
    exit(1);
  }

  // Read y values to vector
  double* yVec;
  yVec = malloc(dataRows * sizeof(double));
  if (yVec==NULL)
  {
    fprintf(stderr, "Error -- out of memory\n");
    exit(1);
  }
  int obsRead = readToDoubleVectorDynamic(dataFile, dataRows,
      valueCol, &yVec);

  // Check for valid first obs
  if (isnan(yVec[0]))
  {
    fprintf(stderr, "Error -- NaN at first observation; aborting\n");
    exit(1);
  }

  // Check for agreement between timesRead and obsRead
  if (timesRead != obsRead)
  {
    fprintf(stderr, "Error -- differing number of entries ");
    fprintf(stderr, "in time and data columns\n");
    exit(1);
  }

  // Read prior to vector
  double * priorVec;
  priorVec = malloc( (basisCols-1) * sizeof(double));
  if (priorVec==NULL)
  {
    fprintf(stderr, "Error -- out of memory\n");
    exit(1);
  }
  readToDoubleVector(priorFile, basisCols-1, 0, priorVec);

  /*
   * Handle coded missing values and restructure times
   */
  int nObs = 0;
  int * validInd;
  validInd = calloc(sizeof(int), obsRead);

  // Find valid indices
  for (i=0; i<obsRead; i++)
  {
    if (yVec[i] != missingCode)
    {
      validInd[nObs]=i;
      nObs++;
    }
  }

  // Check for minimum number of observations
  if (nObs < minObs)
  {
    fprintf(stderr, "Error -- read %d obs, minimum to process is %d\n",
        nObs, minObs);
    exit(1);
  }

  // Move valid data to head of timeVec and yVec
  for (i=0; i<nObs; i++)
  {
    yVec[i] = yVec[validInd[i]];
    timeVec[i] = timeVec[validInd[i]];
  }


  // Setup variables for estimation
  double logPosterior_l, logLikelihood_l, tau_l;
  double logPosterior_m, logLikelihood_m, tau_m;
  double * coef_l, * coef_m;

  coef_l = malloc(kSmooth * sizeof(double));
  if (coef_l==NULL)
  {
    fprintf(stderr, "Error -- out of memory\n");
    exit(1);
  }

  coef_m = malloc(basisCols * sizeof(double));
  if (coef_m==NULL)
  {
    fprintf(stderr, "Error -- out of memory\n");
    exit(1);
  }

  // Run wavelet model with t residuals for full basis
  int iter_m;
  iter_m = lmT(basisMat, basisRows, basisCols,
      yVec, nObs,
      timeVec,
      priorVec,
      nu, basisCols,
      maxIter, tol,
      &logPosterior_m,
      &logLikelihood_m,
      coef_m,
      &tau_m);

  // Run wavelet model with t residuals for restricted (smooth) basis
  int iter_l;
  iter_l = lmT(basisMat, basisRows, basisCols,
      yVec, nObs,
      timeVec,
      priorVec,
      nu, kSmooth,
      maxIter, tol,
      &logPosterior_l,
      &logLikelihood_l,
      coef_l,
      &tau_l);

  // Calculate test statistics (LLR & LPR)
  double llr, lpr;
  llr = logLikelihood_m - logLikelihood_l;
  llr *= 2;
  lpr = logPosterior_m - logPosterior_l;

  /*
   * Print output
   */

  // Basic information
  fprintf(stdout, "%s %d %d %d %g ", idString, nObs, basisCols, kSmooth, nu);

  // Test statistics
  fprintf(stdout, "%g %g ", llr, lpr);

  // Other statistics
  fprintf(stdout, "%g ", sqrt(tau_m));

  // Coefficients for k = 1..m
  for (i=0; i<basisCols; i++) fprintf(stdout, "%g ", coef_m[i]);
  fprintf(stdout, "\n");

  /*
   * Free allocated memory
   */

  // Free basisMat
  free(basisMat);
  basisMat = NULL;

  // Free timeVec & yVec
  free(timeVec);
  timeVec=NULL;

  free(yVec);
  yVec=NULL;

  // Free coef vecs and prior vec
  free(coef_l);
  coef_l=NULL;

  free(coef_m);
  coef_m=NULL;

  free(priorVec);
  priorVec=NULL;

  return 0;
}
