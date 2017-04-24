#ifndef _TOLMIN_H
#define _TOLMIN_H 1

/* Basic types. */
#define INTEGER int     /* used for indexing */
#define REAL    double  /* floating-point type */


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  typedef REAL tolmin_objective(void* ctx, const REAL x[], REAL g[]);

  extern void getmin(tolmin_objective fg, void* ctx, const INTEGER n,
                     const INTEGER m, const INTEGER meq, const REAL a[],
                     const INTEGER ia, const REAL b[], const REAL xl[],
                     const REAL xu[], REAL x[], const REAL acc,
                     INTEGER iact[], INTEGER* nact, REAL par[],
                     const INTEGER iprint, INTEGER* info, REAL w[]);

#if 0

  /* For the record, prototypes of the FORTRAN subroutines. */

  extern void fgcalc_(const INTEGER *n, const REAL *x, REAL *f, REAL *g);

  extern void addcon_(const INTEGER *n, const INTEGER *m, const REAL *a,
                      const INTEGER *ia, INTEGER *iact, INTEGER *nact,
                      REAL *z, REAL *u, const REAL *relacc,
                      const INTEGER *indxbd, REAL *ztc, REAL *cgrad);

  extern void newcon_(const INTEGER *n, const INTEGER *m, const REAL *a,
                      const INTEGER *ia, INTEGER *iact, INTEGER *nact,
                      REAL *z, REAL *u, REAL *d, const REAL *relacc,
                      const INTEGER *mdeg, REAL *zzdiag, REAL *gmnew,
                      REAL *cgrad);

  extern void delcon_(const INTEGER *n, const INTEGER *m, const REAL *a,
                      const INTEGER *ia, INTEGER *iact, INTEGER *nact,
                      REAL *z, REAL *u, const REAL *relacc,
                      const INTEGER *idrop);

  extern void sdirn_(const INTEGER *n, const INTEGER *nact, const REAL *z,
                     REAL *d, REAL *ztg, const REAL *gm, const REAL *relacc,
                     REAL *ddotgm);

  extern void sdegen_(const INTEGER *n, const INTEGER *m, const REAL *a,
                      const INTEGER *ia, INTEGER *iact, INTEGER *nact,
                      REAL *par, REAL *z, REAL *u, REAL *d, REAL *ztg, REAL *gm,
                      const REAL *relacc, REAL *ddotgm, const INTEGER *meql,
                      const INTEGER *mdeg, REAL *gmnew, REAL *parnew,
                      REAL *cgrad);

  extern void getd_(const INTEGER *n, const INTEGER *m, const REAL *a,
                    const INTEGER *ia, INTEGER *iact, INTEGER *nact, REAL *par,
                    const REAL *g, REAL *z, REAL *u, REAL *d, REAL *ztg,
                    const REAL *relacc, REAL *ddotg, const INTEGER *meql,
                    const INTEGER *mdeg, REAL *gm, REAL *gmnew, REAL *parnew,
                    REAL *cgrad);

  extern void stepbd_(const INTEGER *n, const INTEGER *m, const REAL *a,
                      const INTEGER *ia, INTEGER *iact, REAL *bres,
                      const REAL *d, REAL *stepcb, REAL *ddotg,
                      const INTEGER *mdeg, INTEGER *msat, const INTEGER *mtot,
                      INTEGER *indxbd);

  extern void conres_(const INTEGER *n, const INTEGER *m, const REAL *a,
                      const INTEGER *ia, const REAL *b, const REAL *xl,
                      const REAL *xu, const REAL *x, INTEGER *iact,
                      INTEGER *nact, REAL *par, REAL *g, REAL *z, REAL *u,
                      const REAL *xbig, REAL *bres, REAL *d, REAL *ztg,
                      const REAL *relacc, const REAL *tol, REAL *stepcb,
                      REAL *sumres, const INTEGER *meql, INTEGER *msat,
                      const INTEGER *mtot, INTEGER *indxbd, REAL *gm,
                      REAL *gmnew, REAL *parnew, REAL *cgrad);

  extern void ktvec_(const INTEGER *n, const INTEGER *m, const REAL *a,
                     const INTEGER *ia, const INTEGER *iact,
                     const INTEGER *nact, REAL *par, const REAL *g, REAL *reskt,
                     const REAL *z, const REAL *u, const REAL *bres,
                     REAL *relaxf, const INTEGER *meql, REAL *ssqkt,
                     REAL *parw, REAL *resktw);

  extern void lsrch_(const INTEGER *n, REAL *x, REAL *g, const REAL *d,
                     REAL *xs, REAL *gs, const REAL *relacc, const REAL *stepcb,
                     const REAL *ddotg, REAL *f, REAL *step, INTEGER *nfvals,
                     const INTEGER *nfmax, REAL *gopt);

  extern void zbfgs_(const INTEGER *n, const REAL *x, const INTEGER *nact,
                     const REAL *g, REAL *z, REAL *ztg, REAL *xs, REAL *gs,
                     REAL *zznorm);

  extern void minfun_(const INTEGER *n, const INTEGER *m, const REAL *a,
                      const INTEGER *ia, const REAL *b, const REAL *xl,
                      const REAL *xu, REAL *x, const REAL *acc, INTEGER *iact,
                      INTEGER *nact, REAL *par, const INTEGER *iprint,
                      INTEGER *info, REAL *g, REAL *z, REAL *u, REAL *xbig,
                      const REAL *relacc, REAL *zznorm, const REAL *tol,
                      const INTEGER *meql, const INTEGER *mtot, INTEGER *iterc,
                      INTEGER *nfvals, const INTEGER *nfmax, REAL *reskt,
                      REAL *bres, REAL *d, REAL *ztg, REAL *gm, REAL *xs,
                      REAL *gs);

  extern void satact_(const INTEGER *n, const INTEGER *m, const REAL *a,
                      const INTEGER *ia, const REAL *b, const REAL *xl,
                      const REAL *xu, REAL *x, INTEGER *iact, INTEGER *nact,
                      INTEGER *info, REAL *z, REAL *u, REAL *xbig,
                      const REAL *relacc, const REAL *tol, const INTEGER *meql);

  extern void adjtol_(const INTEGER *n, const INTEGER *m, const REAL *a,
                      const INTEGER *ia, const REAL *b, const REAL *xl,
                      const REAL *xu, const REAL *x, const INTEGER *iact,
                      const INTEGER *nact, REAL *xbig, const REAL *relacc,
                      REAL *tol, const INTEGER *meql);

  extern void getfes_(const INTEGER *n, const INTEGER *m, const REAL *a,
                      const INTEGER *ia, const REAL *b, const REAL *xl,
                      const REAL *xu, REAL *x, INTEGER *iact, INTEGER *nact,
                      REAL *par, INTEGER *info, REAL *g, REAL *z, REAL *u,
                      REAL *xbig, const REAL *relacc, REAL *tol,
                      const INTEGER *meql, INTEGER *msat, const INTEGER *mtot,
                      REAL *bres, REAL *d, REAL *ztg, REAL *gm, REAL *gmnew,
                      REAL *parnew, REAL *cgrad);

  extern void initzu_(const INTEGER *n, const INTEGER *m, const REAL *xl,
                      const REAL *xu, REAL *x, INTEGER *iact, INTEGER *meql,
                      INTEGER *info, REAL *z, REAL *u, REAL *xbig,
                      REAL *relacc);

  extern void eqcons_(const INTEGER *n, const INTEGER *m, const INTEGER *meq,
                      const REAL *a, const INTEGER *ia, const REAL *b,
                      const REAL *xu, INTEGER *iact, INTEGER *meql,
                      INTEGER *info, REAL *z, REAL *u, const REAL *relacc,
                      REAL *am, REAL *cgrad);

  extern void minflc_(const INTEGER *n, const INTEGER *m, const INTEGER *meq,
                      const REAL *a, const INTEGER *ia, const REAL *b,
                      const REAL *xl, const REAL *xu, REAL *x, const REAL *acc,
                      INTEGER *iact, INTEGER *nact, REAL *par,
                      const INTEGER *iprint, INTEGER *info, REAL *g, REAL *z,
                      REAL *u, REAL *xbig, REAL *reskt, REAL *bres, REAL *d,
                      REAL *ztg, REAL *gm, REAL *xs, REAL *gs);


  extern void getmin_(const INTEGER *n, const INTEGER *m, const INTEGER *meq,
                      const REAL *a, const INTEGER *ia, const REAL *b,
                      const REAL *xl, const REAL *xu, REAL *x, const REAL *acc,
                      INTEGER *iact, INTEGER *nact, REAL *par,
                      const INTEGER *iprint, INTEGER *info, REAL *w);

#endif

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _TOLMIN_H */
