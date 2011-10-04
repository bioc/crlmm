#include <R.h>
#include <R_ext/Rdynload.h>
#include "crlmm.h"

static const R_CallMethodDef CallEntries[] = {
    {"gtypeCallerPart1", (DL_FUNC)&gtypeCallerPart1, 17},
    {"gtypeCallerPart2", (DL_FUNC)&gtypeCallerPart2, 19},
    {"normalizeBAF", (DL_FUNC)&normalizeBAF, 2},
    {NULL, NULL, 0}
};

void R_init_crlmm(DllInfo *dll){
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
}
