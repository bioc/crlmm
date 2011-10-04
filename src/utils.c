#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

int genotypeConfidence(const double *prob){
  int K=1000;
  if (*prob == 1.0) return(INT_MAX);
  else return( (int) round(-K*log2(1- *prob)));
}

int intInSet(const int *x, const int *set, const int *n){
  int i;
  for (i=0; i < *n; i++)
    if (*x == set[i]) return(1);
  return(0);
}

int genotypeCall(const double *pAA, const double *pAB, const double *pBB){
  if( *pAA >= *pAB && *pAA >= *pBB ) return(1);
  else if (*pAB > *pAA && *pAB >= *pBB) return(2);
  else return(3);
}

double sdCorrection(const int *n){
  return(sqrt(1.0+1.0/fmax(1.0, (double) *n)));
}

int sort_double(const double *a1, const double *a2){
  if (*a1 < *a2)
    return (-1);
  if (*a1 > *a2)
    return (1);
  return 0;
}

void trimmed_mean(double *datavec, int *classvec, int class, double trim, int cols, int rows, double *m1, double *m2, double *m3, int i_ext){
  double sum=0, sum2=0;
  int i, j=0, n_ignore, n=0;

  for (i = 0; i < cols; i++)
    if (classvec[i] == class)
      n++;

  double *buffer=Calloc(n, double);
  for (i = 0; i < cols; i++)
    if (classvec[i] == class){
      buffer[j]=datavec[i];
      j++;
    }
  qsort(buffer, n, sizeof(double), (int(*)(const void*, const void*))sort_double);
  n_ignore= (int) floor((double) n * trim);
  j=0;
  for (i = n_ignore; i < (n-n_ignore); i++){
    sum+=buffer[i];
    sum2+=pow(buffer[i], 2);
    j++;
  }
  sum/=j;
  sum2-=(pow(sum, 2)*(double) j);
  sum2/=(j-1);
  sum2=sqrt(sum2);
  m1[i_ext + (class-1) * rows]=sum;
  m2[i_ext + (class-1) * rows]=sum2;
  m3[i_ext + (class-1) * rows]=j;
  Free(buffer);
}

void trimmed_stats(double *data, double *m1, double *m2, double *m3, int *class, int rows, int cols, double *trim){
  int i, j, n1, n2, n3;
  double *datvec=Calloc(cols,double);
  int *classvec=Calloc(cols,int);

  for (i=0; i < rows; i++){
    n1=0;
    n2=0;
    n3=0;

    for (j=0; j < cols; j++){
      if (class[j*rows + i] == 1){
	datvec[j]=data[j*rows + i];
	++n1;
	classvec[j] = 1;
      } else if (class[j*rows + i] == 2){
	datvec[j]=data[j*rows + i];
	++n2;
	classvec[j] = 2;
      } else if (class[j*rows + i] == 3){
	datvec[j]=data[j*rows + i];
	++n3;
	classvec[j] = 3;
      } else {
	// Should be the NA's
	classvec[j] = class[j*rows + i];
      }
    }
    trimmed_mean(datvec, classvec, 1, trim[0], cols, rows, m1, m2, m3, i);
    trimmed_mean(datvec, classvec, 2, trim[0], cols, rows, m1, m2, m3, i);
    trimmed_mean(datvec, classvec, 3, trim[0], cols, rows, m1, m2, m3, i);
  }
  Free(datvec);
  Free(classvec);
}

double  median(double *x, int length){
  int half;
  double med;
  double *buffer = Calloc(length,double);
  memcpy(buffer,x,length*sizeof(double));
  half = (length + 1)/2;
  rPsort(buffer, length, half-1);
  med = buffer[half-1];
  if (length % 2 == 0){
    rPsort(buffer, length, half);
    med = (med + buffer[half])/2.0;
  }
  Free(buffer);
  return med;
}

void mad_median(double *datavec, int *classvec, int class, double trim, int cols, int rows, double *m1, double *m2, double *m3, int i_ext){
  /* trim is ignored for the moment - for compatibility */
  int i, j=0, n_ignore, n=0;

  for (i = 0; i < cols; i++)
    if (classvec[i] == class)
      n++;

  double *buffer=Calloc(n, double);

  for (i = 0; i < cols; i++)
    if (classvec[i] == class){
       buffer[j]=datavec[i];
      j++;
    }
  m1[i_ext + (class-1) * rows] = median(buffer, n);
  for (i = 0; i < n; i++)
    buffer[i] = fabs(buffer[i]-m1[i_ext + (class-1) * rows]);
  m2[i_ext + (class-1) * rows] = median(buffer, n);
  m3[i_ext + (class-1) * rows] = n;
  Free(buffer);
}

SEXP normalizeBAF(SEXP theta, SEXP cTheta){
  /*
    ARGUMENTS:
    theta.: N x C matrix with estimated \theta
    cTheta: N x 3 matrix with canonical \thetas (AA, AB, BB)

    VALUE:
    baf: N x C matrix with normalized \theta
  */

  SEXP baf;
  double *p2baf, *p2theta, *p2ctheta;
  int i, j, idx, rowsT, rowsCT, colsT, colsCT;
  rowsT = INTEGER(getAttrib(theta, R_DimSymbol))[0];
  rowsCT = INTEGER(getAttrib(cTheta, R_DimSymbol))[0];
  if (rowsT != rowsCT)
    error("Number of rows of 'theta' must match number of rows of 'cTheta'\n");
  colsCT = INTEGER(getAttrib(cTheta, R_DimSymbol))[1];
  if (colsCT != 3)
    error("'cTheta' must have 3 columns: AA, AB and BB\n");
  colsT = INTEGER(getAttrib(theta, R_DimSymbol))[1];

  PROTECT(baf = allocMatrix(REALSXP, rowsT, colsT));
  p2baf = NUMERIC_POINTER(baf);
  p2theta = NUMERIC_POINTER(theta);
  p2ctheta = NUMERIC_POINTER(cTheta);
  for (i=0; i < rowsT; i++){
    for (j=0; j < colsT; j++){
      idx = i + j*rowsT;
      if (ISNA(p2theta[idx]) || ISNA(p2ctheta[i]) || ISNA(p2ctheta[i+rowsT]) || ISNA(p2ctheta[i+2*rowsT])){
	p2baf[idx] = NA_REAL;
      }else if (p2theta[idx] < p2ctheta[i]){
	p2baf[idx] = 0;
      }else if (p2theta[idx] >= p2ctheta[i] & p2theta[idx] < p2ctheta[i + rowsT]){
	p2baf[idx] = .5*(p2theta[idx]-p2ctheta[i])/(p2ctheta[i+rowsT]-p2ctheta[i]);
      }else if(p2theta[idx] >= p2ctheta[i+rowsT] & p2theta[idx] < p2ctheta[i + 2*rowsT]){
	p2baf[idx] = .5+.5*(p2theta[idx]-p2ctheta[i+rowsT])/(p2ctheta[i+2*rowsT]-p2ctheta[i+rowsT]);
      }else{
	p2baf[idx] = 1;
      }
    }
  }

  UNPROTECT(1);
  return(baf);
}


/* Pieces below are for testing */

static void mad_stats(double *data, double *m1, double *m2, double *m3, int *class, int rows, int cols, double *trim){
  int i, j, n1, n2, n3;
  double *datvec=Calloc(cols,double);
  int *classvec=Calloc(cols,int);

  for (i=0; i < rows; i++){
    n1=0;
    n2=0;
    n3=0;

    for (j=0; j < cols; j++){
      if (class[j*rows + i] == 1){
	datvec[j]=data[j*rows + i];
	++n1;
	classvec[j] = 1;
      } else if (class[j*rows + i] == 2){
	datvec[j]=data[j*rows + i];
	++n2;
	classvec[j] = 2;
      } else if (class[j*rows + i] == 3){
	datvec[j]=data[j*rows + i];
	++n3;
	classvec[j] = 3;
      } else {
	// Should be the NA's
	classvec[j] = class[j*rows + i];
      }
    }
    mad_median(datvec, classvec, 1, trim[0], cols, rows, m1, m2, m3, i);
    mad_median(datvec, classvec, 2, trim[0], cols, rows, m1, m2, m3, i);
    mad_median(datvec, classvec, 3, trim[0], cols, rows, m1, m2, m3, i);
  }
  Free(datvec);
  Free(classvec);
}


SEXP test_mad_median(SEXP X, SEXP Y, SEXP trim){
  SEXP dim1;
  SEXP estimates1, estimates2, estimates3, output;
  double *Xptr, *Mptr1, *Mptr2, *Mptr3, *Tptr;
  int *Yptr;
  int rows, cols;

  PROTECT(dim1 = getAttrib(X,R_DimSymbol));
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1];

  Xptr = NUMERIC_POINTER(AS_NUMERIC(X));
  Yptr = INTEGER_POINTER(AS_INTEGER(Y));
  Tptr = NUMERIC_POINTER(AS_NUMERIC(trim));

  PROTECT(estimates1 = allocMatrix(REALSXP, rows, 3));
  PROTECT(estimates2 = allocMatrix(REALSXP, rows, 3));
  PROTECT(estimates3 = allocMatrix(REALSXP, rows, 3));

  Mptr1 = NUMERIC_POINTER(estimates1);
  Mptr2 = NUMERIC_POINTER(estimates2);
  Mptr3 = NUMERIC_POINTER(estimates3);

  mad_stats(Xptr, Mptr1, Mptr2, Mptr3, Yptr, rows, cols, Tptr);

  PROTECT(output = allocVector(VECSXP,3));
  SET_VECTOR_ELT(output, 0, estimates1);
  SET_VECTOR_ELT(output, 1, estimates2);
  SET_VECTOR_ELT(output, 2, estimates3);

  UNPROTECT(5);

  return output;
}
