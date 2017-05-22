#ifndef _LINCOA_H
#define _LINCOA_H 1

/* Basic types. */
#define INTEGER int     /* used for indexing */
#define REAL    double  /* floating-point type */


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

extern REAL calfun(const INTEGER n, const REAL x[]);

extern void lincoa(const INTEGER n,
                   const INTEGER npt,
                   const INTEGER m,
                   const REAL a[],
                   const INTEGER ia,
                   const REAL b[],
                   REAL x[],
                   const REAL rhobeg,
                   const REAL rhoend,
                   const INTEGER iprint,
                   const INTEGER maxfun,
                   void* ws);

/**
 * Get the size of the workspace for LINCOA.
 *
 * @param n   the number of variables.
 * @param npt the number of points.
 * @param m   the number of constraints.
 *
 * @return The number of bytes needed by LINCOA workspace (assuming the first
 *         bytes are correctly aligned for integers).
 */
extern size_t
lincoa_storage(INTEGER n, INTEGER npt, INTEGER m);

#if 0
/* For the record, prototypes of the FORTRAN subroutines. */
extern int
calfun_(const INTEGER *n, const REAL *x, REAL *f);

extern int
update_(const INTEGER *n,
        const INTEGER *npt,
        const REAL *xpt,
        REAL *bmat,
        REAL *zmat,
        INTEGER *idz,
        const INTEGER *ndim,
        const REAL *sp,
        const REAL *step,
        const INTEGER *kopt,
        INTEGER *knew,
        REAL *vlag,
        REAL *w);

extern int
prelim_(const INTEGER *n,
        const INTEGER *npt,
        const INTEGER *m,
        const REAL *amat,
        REAL *b,
        REAL *x,
        const REAL *rhobeg,
        const INTEGER *iprint,
        REAL *xbase,
        REAL *xpt,
        REAL *fval,
        REAL *xsav,
        REAL *xopt,
        REAL *gopt,
        INTEGER *kopt,
        REAL *hq,
        REAL *pq,
        REAL *bmat,
        REAL *zmat,
        INTEGER *idz,
        const INTEGER *ndim,
        REAL *sp,
        REAL *rescon,
        REAL *step,
        REAL *pqw,
        REAL *w);

extern int
getact_(const INTEGER *n,
        const INTEGER *m,
        const REAL *amat,
        const REAL *b,
        INTEGER *nact,
        INTEGER *iact,
        REAL *qfac,
        REAL *rfac,
        const REAL *snorm,
        REAL *resnew,
        REAL *resact,
        const REAL *g,
        REAL *dw,
        REAL *vlam,
        REAL *w);

extern int
trstep_(const INTEGER *n,
        const INTEGER *npt,
        const INTEGER *m,
        const REAL *amat,
        const REAL *b,
        const REAL *xpt,
        const REAL *hq,
        const REAL *pq,
        INTEGER *nact,
        INTEGER *iact,
        const REAL *rescon,
        REAL *qfac,
        REAL *rfac,
        REAL *snorm,
        REAL *step,
        REAL *g,
        REAL *resnew,
        REAL *resact,
        REAL *d,
        REAL *dw,
        REAL *w);

extern int
qmstep_(const INTEGER *n,
        const INTEGER *npt,
        const INTEGER *m,
        const REAL *amat,
        const REAL *b,
        const REAL *xpt,
        const REAL *xopt,
        const INTEGER *nact,
        const INTEGER *iact,
        const REAL *rescon,
        const REAL *qfac,
        const INTEGER *kopt,
        const INTEGER *knew,
        const REAL *del,
        REAL *step,
        REAL *gl,
        const REAL *pqw,
        REAL *rstat,
        REAL *w,
        INTEGER *ifeas);

extern int
lincob_(const INTEGER *n,
        const INTEGER *npt,
        const INTEGER *m,
        const REAL *amat,
        REAL *b,
        REAL *x,
        const REAL *rhobeg,
        const REAL *rhoend,
        const INTEGER *iprint,
        const INTEGER *maxfun,
        REAL *xbase,
        REAL *xpt,
        REAL *fval,
        REAL *xsav,
        REAL *xopt,
        REAL *gopt,
        REAL *hq,
        REAL *pq,
        REAL *bmat,
        REAL *zmat,
        const INTEGER *ndim,
        REAL *step,
        REAL *sp,
        REAL *xnew,
        INTEGER *iact,
        REAL *rescon,
        REAL *qfac,
        REAL *rfac,
        REAL *pqw,
        REAL *w);

extern int
lincoa_(const INTEGER *n,
        const INTEGER *npt,
        const INTEGER *m,
        const REAL *a,
        const INTEGER *ia,
        const REAL *b,
        REAL *x,
        const REAL *rhobeg,
        const REAL *rhoend,
        const INTEGER *iprint,
        const INTEGER *maxfun,
        REAL *w);
#endif

  #ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _LINCOA_H */
