#include "SparseM.h"

#include <R_ext/Rdynload.h>

static const R_FortranMethodDef FortranEntries[] = {
    {"aedib",    (DL_FUNC) &F77_SUB(aedib),    16},
    {"aeexpb",   (DL_FUNC) &F77_SUB(aeexpb),   16},
    {"aemub",    (DL_FUNC) &F77_SUB(aemub),    15},
    {"amub",     (DL_FUNC) &F77_SUB(amub),     15},
    {"amubdg",   (DL_FUNC) &F77_SUB(amubdg),   10},
    {"amux",     (DL_FUNC) &F77_SUB(amux),      6},
    {"aplsb",    (DL_FUNC) &F77_SUB(aplsb),    16},
    {"bckslb",   (DL_FUNC) &F77_SUB(bckslb),   16},
    {"bckslf",   (DL_FUNC) &F77_SUB(bckslf),   16},
    {"bckslv",   (DL_FUNC) &F77_SUB(bckslv),   16},
    {"chol",     (DL_FUNC) &F77_SUB(chol),     32},
    {"chol2csr", (DL_FUNC) &F77_SUB(chol2csr), 12},
    {"coocsr",   (DL_FUNC) &F77_SUB(coocsr),    8},
    {"cscssc",   (DL_FUNC) &F77_SUB(cscssc),    9},
    {"csr",      (DL_FUNC) &F77_SUB(csr),       8},
    {"csrcoo",   (DL_FUNC) &F77_SUB(csrcoo),   11},
    {"csrcsc2",  (DL_FUNC) &F77_SUB(csrcsc2),  10},
    {"csrdns",   (DL_FUNC) &F77_SUB(csrdns),    8},
    {"csrssr",   (DL_FUNC) &F77_SUB(csrssr),    9},
    {"filter1",  (DL_FUNC) &F77_SUB(filter1),  11},
    {"nzero",    (DL_FUNC) &F77_SUB(nzero),    11},
    {"ssrcsr",   (DL_FUNC) &F77_SUB(ssrcsr),   13},
    {NULL, NULL, 0}
};

void R_init_SparseM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
