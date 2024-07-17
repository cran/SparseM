#include <R.h>
/* for int_least32_t : */
#include <stdint.h>

/* .Fortran() calls from R code :*/

// bckslv.f :
void F77_NAME(bckslb)(int *m, int *nsubmax, int *nsuper, int *nrhs, int *lindx, int *xlindx, int *nnzlmax, double *lnz,
		      int *xlnz, int *invp, int *perm, int *xsuper, double *newrhs, double *sol, double *b, double *timed);
void F77_NAME(bckslf)(int *m, int *nsubmax, int *nsuper, int *nrhs, int *lindx, int *xlindx, int *nnzlmax, double *lnz,
		      int *xlnz, int *invp, int *perm, int *xsuper, double *newrhs, double *sol, double *b, double *timed);
void F77_NAME(bckslv)(int *m, int *nsubmax, int *nsuper, int *nrhs, int *lindx, int *xlindx, int *nnzlmax, double *lnz,
		      int *xlnz, int *invp, int *perm, int *xsuper, double *newrhs, double *sol, double *b, double *timed);
// chol.f :
void F77_NAME(chol)(int *m, int *nnzdmax, double *d, int *jd, int *id, /* 5 */
		    int *nnzdsm, double *dsub, int *jdsub, int *nsub, int *nsubmax, /* 10 */
		    int *lindx, int *xlindx, int *nsuper, int *nnzlmax, double *lnz, int *xlnz, int *invp, int *perm,/* 18 */
		    int *iwmax, int *iwork, int *colcnt, int *snode, int *xsuper, int *split, int *tmpmax, /* 25 */
		    double *tmpvec, int *cachsz, int *level, int *ierr, /* 29 */
		    double *timed, double *tiny, double *Large); /* 32 */
// chol2csr.f :
void F77_NAME(chol2csr)(int *nrow, int *nnzlindx, int *nsuper, int *lindx, int *xlindx, int *nnzl, double *lnz,
			int *xlnz, int *dim, double *ra, int *ia, int *ja);
// csr.f :
void F77_NAME(csr)(double *a, double *ra, int *ja, int *ia, int *m, int *n, int *nnz, double *eps);
void F77_NAME(nzero)(int *ja, int *ia, int *nrow, int *ncol, int *nnz, int *nz,
		     double *rao, int *jao, int *iao);

// sparskit.f : ------------------------------------------------------------

// 0) Unary subroutines module
void F77_NAME(filter1)(int *n, int *rel, double *drptol, double *a, int *ja, int *ia, double *b, int *jb, int *ib,
		       int *len, int *ierr);
void F77_NAME(amubdg)(int *nrow, int *ncol, int *ncolb, int *ja, int *ia, int *jb, int *ib, int *ndegr, int *nnz, int *iw);
// 1) BASIC LINEAR ALGEBRA for Sparse Matrices.   BLASSM module
void F77_NAME(aedib)(int *nrow, int *ncol, int *job, double *a, int *ja, int *ia, double *b, int *jb, int *ib,
		     double *c, int *jc, int *ic, int *nzmax, int *iw, double *aw, int *ierr);
void F77_NAME(aeexpb)(int *nrow, int *ncol, int *job, double *a, int *ja, int *ia, double *b, int *jb, int *ib,
		      double *c, int *jc, int *ic, int *nzmax, int *iw, double *aw, int *ierr);
void F77_NAME(aemub)(int *nrow, int *ncol, double *a, int *ja, int *ia, double *amask, int *jmask, int *imask,
		     double *c, int *jc, int *ic, int *nzmax, int *ierr);
void F77_NAME(amub)(int *nrow, int *ncol, int *job, double *a, int *ja, int *ia, double *b, int *jb, int *ib,
		    double *c, int *jc, int *ic, int *nzmax, int *iw, int *ierr);
void F77_NAME(amux)(int *n, double *x, double *y, double *a, int *ja, int *ia);
void F77_NAME(aplsb)(int *nrow, int *ncol, int *job, double *a, int *ja,int *ia, double *s, double *b, int *jb,int *ib,
		     double *c, int *jc, int *ic, int *nzmax, int *iw, int *ierr);
// 2) Format Conversion module
void F77_NAME(coocsr)(int *nrow, int *nnz, double *a, int *ir, int *jc, double *ao, int *jao, int *iao);
void F77_NAME(cscssc)(int *ncol, double *a, int *ja,int *ia, int *nzmax, double *ao, int *jao,int *iao, int *ierr);
void F77_NAME(csrcoo)(int *nrow, int *job, int *nzmax, double *a, int *ja, int *ia, int *nnz,
		      double *ao, int *ir, int *jc, int *ierr);
void F77_NAME(csrcsc2)(int *n, int *n2, int *job, int *ipos, double *a, int *ja, int *ia,
		       double *ao, int *jao, int *iao);
void F77_NAME(csrdns)(int *nrow, int *ncol, double *a, int *ja, int *ia, double *dns, int *ndns, int *ierr);
void F77_NAME(csrssr)(int *nrow, double *a, int *ja, int *ia, int *nzmax, double *ao, int *jao, int *iao, int *ierr);
void F77_NAME(ssrcsr)(int *job, int *value2, int *nrow, double *a, int *ja, int *ia, int *nzmax,
		      double *ao, int *jao, int *iao, int *indu, int *iwk, int *ierr);
