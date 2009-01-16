#include <R.h>
#include <R_ext/Rdynload.h>
#include "crlmm.h"

static const R_CallMethodDef CallEntries[] = {
    {"test", (DL_FUNC)&test, 1},
    {"gtypeCallerPart1nm", (DL_FUNC)&gtypeCallerPart1nm, 17},
    {"gtypeCallerPart2nm", (DL_FUNC)&gtypeCallerPart2nm, 19},
    {"gtypeCallerPart1NormalNoN", (DL_FUNC)&gtypeCallerPart1NormalNoN, 17},
    {"gtypeCallerPart1", (DL_FUNC)&gtypeCallerPart1, 17},
    {"gtypeCallerPart2", (DL_FUNC)&gtypeCallerPart2, 19},
    {"gtypeCallerPart1TNoN", (DL_FUNC)&gtypeCallerPart1TNoN, 17},
    {"gtypeCallerPart2TNoN", (DL_FUNC)&gtypeCallerPart2TNoN, 19},
    {NULL, NULL, 0}
};

void R_init_crlmm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);    
}



