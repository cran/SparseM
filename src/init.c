#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(aedib)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(aeexpb)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(aemub)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(amub)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(amubdg)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(amux)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(aplsb)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(bckslb)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(bckslf)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(bckslv)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(chol)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(chol2csr)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(coocsr)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(cscssc)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(csr)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(csrcoo)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(csrcsc2)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(csrdns)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(csrssr)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(filter1)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(nzero)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ssrcsr)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"aedib",    (DL_FUNC) &F77_NAME(aedib),    16},
    {"aeexpb",   (DL_FUNC) &F77_NAME(aeexpb),   16},
    {"aemub",    (DL_FUNC) &F77_NAME(aemub),    15},
    {"amub",     (DL_FUNC) &F77_NAME(amub),     15},
    {"amubdg",   (DL_FUNC) &F77_NAME(amubdg),   10},
    {"amux",     (DL_FUNC) &F77_NAME(amux),      6},
    {"aplsb",    (DL_FUNC) &F77_NAME(aplsb),    16},
    {"bckslb",   (DL_FUNC) &F77_NAME(bckslb),   16},
    {"bckslf",   (DL_FUNC) &F77_NAME(bckslf),   16},
    {"bckslv",   (DL_FUNC) &F77_NAME(bckslv),   16},
    {"chol",     (DL_FUNC) &F77_NAME(chol),     30},
    {"chol2csr", (DL_FUNC) &F77_NAME(chol2csr), 12},
    {"coocsr",   (DL_FUNC) &F77_NAME(coocsr),    8},
    {"cscssc",   (DL_FUNC) &F77_NAME(cscssc),    9},
    {"csr",      (DL_FUNC) &F77_NAME(csr),       8},
    {"csrcoo",   (DL_FUNC) &F77_NAME(csrcoo),   11},
    {"csrcsc2",  (DL_FUNC) &F77_NAME(csrcsc2),  10},
    {"csrdns",   (DL_FUNC) &F77_NAME(csrdns),    8},
    {"csrssr",   (DL_FUNC) &F77_NAME(csrssr),    9},
    {"filter1",  (DL_FUNC) &F77_NAME(filter1),  11},
    {"nzero",    (DL_FUNC) &F77_NAME(nzero),    11},
    {"ssrcsr",   (DL_FUNC) &F77_NAME(ssrcsr),   13},
    {NULL, NULL, 0}
};

void R_init_SparseM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
