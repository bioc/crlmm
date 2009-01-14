#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

SEXP test (SEXP x){
  double *xptr, buffer;
  int n, i;
  SEXP res;
  n = GET_LENGTH(x);
  Rprintf("There are %d elements.\n");
  PROTECT(res = allocVector(INTSXP, n));
  xptr = NUMERIC_POINTER(AS_NUMERIC(x));
  for (i=0; i<n; i++){
    buffer = xptr[i];
    INTEGER(res)[i] = genotypeConfidence(&buffer);
  }
  UNPROTECT(1);
  return(res);
}
