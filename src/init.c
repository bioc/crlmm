#include <R.h>
#include <R_ext/Rdynload.h>
#include "crlmm.h"
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

static const R_CallMethodDef CallEntries[] = {
    {"gtypeCallerPart1", (DL_FUNC)&gtypeCallerPart1, 17},
    {"gtypeCallerPart2", (DL_FUNC)&gtypeCallerPart2, 19},
    {"normalizeBAF", (DL_FUNC)&normalizeBAF, 2},
    {NULL, NULL, 0}
};

void R_init_crlmm(DllInfo *dll){
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
}


SEXP subColSummarizeMedianPP(SEXP RMatrix, SEXP R_rowIndexList){
  static SEXP(*fun)(SEXP, SEXP) = NULL;
  if (fun == NULL)
    fun =  (SEXP(*)(SEXP, SEXP))R_GetCCallable("preprocessCore","R_subColSummarize_median");
  return fun(RMatrix, R_rowIndexList);
}
