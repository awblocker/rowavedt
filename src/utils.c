
//#include <stdlib.h>
//#include <math.h>
//#include <stdio.h>
//#include <string.h>

#include "rowavedt.h"

// Function to catch NULL pointers and exit if needed
void checkPtr(const void * ptr, const char * msg)
{
  if (ptr==NULL)
  {
    fprintf(stderr, "Error -- %s\n", msg);
    exit(1);
  }
}

int arrayMinMax(double * X, int n, double * min, double * max)
{
  // Initialization
  *min = X[0];
  *max = X[0];

  if (n > 1)
  {
    int i;
    for (i=1; i<n; i++)
    {
      if (X[i] > *max)
      {
        *max = X[i];
      }
      else if (X[i] < *min)
      {
        *min = X[i];
      }
    }
  }

  return 0;
}

int allocateMatrix(double*** matrix, int nRows, int nCols)
{
  int i;
  (*matrix) = malloc(nRows * sizeof(double*));
  if ((*matrix) == NULL)
    return 1;
  for (i=0; i<nRows; i++)
  {
    (*matrix)[i] = malloc(nCols * sizeof(double));
    if ((*matrix)[i] == NULL)
      return 1;
  }
  return 0;
}

void freeMatrix(double*** matrix, int nRows, int nCols)
{
  int i;
  for (i=0; i<nCols; i++)
  {
    free((*matrix)[i]);
  }
  free(*matrix);
  *matrix = NULL;
}

double arraymean(double x[], int start, int stop)
{
  int n,i;
  double xbar=0;

  n = stop-start;

  for (i=start; i<stop; i++)
  {
    xbar += x[i]/n;
  }

  return xbar;
}

double arrayvar(double x[], int start, int stop) {
  // Computes variance of an array using Welford's algorithm
  int n, k;
  double m_new=0, m_old=0, s=0;

  n = stop-start;

  for (k=start; k<stop; k++) {
    m_new = m_old + (x[k]-m_old)/(k+1);
    s = s + (x[k]-m_old)*(x[k]-m_new);
    m_old = m_new;
  }

  return s/n;
}

double bisect(double x, double data[], int start, int end)
  /* Iterative function to execute bisection search. Returns double.
   * On failure to find key, returns midpoint of indices.*/
{
  int mid;
  while (start <= end) {
    mid = (start + end)/2;

    if (x < data[mid])
      end = mid - 1;
    else if (x > data[mid])
      start = mid + 1;
    else
      return mid;
  }
  return (double) (start+end)/2;
}

double quantile(double x, double data[], int n)
  /* Function to calculate quantile of value from sorted value
   * using bisection algorithm*/
{
  if (x < data[0])
    return 0;
  else if (x >= data[n-1])
    return 1;

  return (1+bisect(x, data, 0, n))/n;
}

double quantile_int(double x, double data[], int n)
  /* Function to calculate quantile of value from sorted value
   * using bisection algorithm. Returns index instead of pvalue.*/
{
  if (x < data[0])
    return 0;
  else if (x >= data[n-1])
    return n;

  return (1+bisect(x, data, 0, n));
}

int compare_dbl(const void * a, const void * b)
  /* Comparator function for qsort */
{
  return (*(double*)a - *(double*)b);
}

// Column-major ordering for BLAS/LAPACK compatibility
void readToDoubleMatrix(const char * fname, int nRows, int nCols,
    double * X)
{
  // Initialize pointers
  FILE * infile;
  const char sep[] = "\t ,";

  // Initialize buffers and counters
  char row[(nCols+1)*100];
  char * buf;
  int i=0, j=0;

  // Open file
  infile = fopen(fname, "r");
  checkPtr(infile, fname);

  // Iterate over lines of file
  while(fgets(row, (nCols+1)*100, infile)!=NULL && i<nRows)
  {
    // Skip comments
    if (row[0] == '#') continue;

    // Read columns from data
    j = 0;
    buf = strtok(row, sep);

    while (buf != NULL && j<nCols)
    {
      X[i + j*nRows] = atof(buf);
      j++;
      if (buf != NULL)
      {
        buf = strtok(NULL, sep);
      }
    }
    i++;
  }

  // Close file
  fclose(infile);
}

// Column-major ordering for BLAS/LAPACK compatibility
void readToDoubleMatrixPtr (const char * fname, int nRows, int nCols,
    double ** X)
{
  // Initialize pointers
  FILE * infile;
  const char sep[] = "\t ,";

  // Initialize buffers and counters
  char row[(nCols+1)*100];
  char * buf;
  int i=0, j=0;

  // Open file
  infile = fopen(fname, "r");
  checkPtr(infile, fname);

  // Iterate over lines of file
  while(fgets(row, (nCols+1)*100, infile)!=NULL && i<nRows)
  {
    // Skip comments
    if (row[0] == '#') continue;

    // Read columns from data
    j = 0;
    buf = strtok(row, sep);

    while (buf != NULL && j<nCols)
    {
      buf = strtok(NULL, sep);
      if (buf != NULL)
      {
        X[i][j] = atof(buf);
        j++;
      }
    }
    i++;
  }

  // Close file
  fclose(infile);
}

int readToDoubleVector (const char * fname, int nRows, int col, double* X)
{
  // Initialize pointers
  FILE * infile;
  const char sep[] = "\t ,";

  // Initialize buffers and counters
  char row[(col+1)*100];
  char * buf;
  int i=0, j=0;

  // Open file
  infile = fopen(fname, "r");
  checkPtr(infile, fname);

  // Iterate over lines of file
  while(fgets(row, (col+1)*100, infile)!=NULL && i<nRows)
  {
    // Skip comments
    if (row[0] == '#') continue;

    // Read columns from data
    j = 0;
    buf = strtok(row, sep);

    while (buf != NULL)
    {
      if (buf != NULL && j==col)
      {
        X[i]= atof(buf);
        break;
      }
      j++;
      buf = strtok(NULL, sep);
    }
    i++;
  }

  // Close file
  fclose(infile);

  // Return number of lines read
  return i;
}

// Dynamically resize vector X to hold data; return final size as int
int readToDoubleVectorDynamic(const char * fname, int startRows, int col,
    double** X)
{
  // Initialize pointers
  FILE * infile;
  const char sep[] = "\t ,";

  // Initialize buffers and counters
  char row[(col+1)*100];
  char * buf;
  int i=0, j=0;
  int size = startRows;
  double * tmp;

  // Open file
  infile = fopen(fname, "r");
  checkPtr(infile, fname);

  // Iterate over lines of file
  while(fgets(row, (col+1)*100, infile)!=NULL)
  {
    // Skip comments
    if (row[0] == '#') continue;

    // Read columns from data
    j = 0;
    buf = strtok(row, sep);

    while (buf != NULL)
    {
      if (buf != NULL && j==col)
      {
        // Check for overflow of current array and reallocate if needed
        if (i >= size)
        {
          size = size * 2;
          tmp = realloc((*X), size * sizeof(double));
          if (tmp != NULL)
          {
            (*X) = tmp;
          }
          else
          {
            fprintf(stderr, "Error -- out of memory\n");
            exit(1);
          }
        }
        (*X)[i]= atof(buf);
        break;
      }
      j++;
      buf = strtok(NULL, sep);
    }
    i++;
  }

  // Close file
  fclose(infile);

  // Return number of lines read
  return i;
}
