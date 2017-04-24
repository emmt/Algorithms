#ifndef _TOLMIN_H
#define _TOLMIN_H 1

/* Basic types. */
#define INTEGER long    /* used for indexing */
#define REAL    double  /* floating-point type */


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  typedef REAL tolmin_objective(void* ctx, const REAL x[], REAL g[]);

  extern int tolmin(tolmin_objective fg, void* ctx, const INTEGER n,
                    const INTEGER m, const INTEGER meq, const REAL a[],
                    const INTEGER ia, const REAL b[], const REAL xl[],
                    const REAL xu[], REAL x[], const REAL acc,
                    INTEGER iact[], INTEGER* nact, REAL par[],
                    const INTEGER iprint, INTEGER nfmax, REAL w[]);

  /* FORTRAN WRAPPER */

  extern void fgcalc_(const INTEGER *n, const REAL *x, REAL *f, REAL *g);


  /* This is the entry point to a package of subroutines that calculate the
   *    the least value of a differentiable function of several variables
   *    subject to linear constraints on the values of the variables, using
   *    the method that is described in the paper "A tolerant algorithm for
   *    linearly constrained optimization calculations", Math. Programming B,
   *    Vol. 45, pp. 547-566 (1989).
   *
   * N is the number of variables and must be set by the user.
   * M is the number of linear constraints (excluding simple bounds) and
   *    must be set by the user.
   * MEQ is the number of constraints that are equalities and must be set
   *    by the user.
   * A(.,.) is a 2-dimensional array whose columns are the gradients of the
   *    M constraint functions.  Its entries must be set by the user and
   *    its dimensions must be at least N by M.
   * IA is the actual first dimension of the array A that is supplied by the
   *    user, so its value may not be less than N.
   * B(.) is a vector of constraint right hand sides that must also be set
   *    by the user.  Specifically the constraints on the variables X(I)
   *    I=1(1)N are
   *         A(1,K)*X(1)+...+A(N,K)*X(N) .EQ. B(K)  K=1,...,MEQ
   *         A(1,K)*X(1)+...+A(N,K)*X(N) .LE. B(K)  K=MEQ+1,...,M  .
   *    Note that the data that define the equality constraints come before
   *    the data of the inequalities.
   * XL(.) and XU(.) are vectors whose components must be set to lower and
   *    upper bounds on the variables.  Choose very large negative and
   *    positive entries if a component should be unconstrained, or set
   *    XL(I)=XU(I) to freeze the I-th variable.  Specifically these simple
   *    bounds are
   *         XL(I) .LE. X(I) and X(I) .LE. XU(I)  I=1,...,N  .
   * X(.) is the vector of variables of the optimization calculation.  Its
   *    initial elements must be set by the user to an estimate of the
   *    required solution.  The subroutines can usually cope with poor
   *    estimates, and there is no need for X(.) to be feasible initially.
   *    These variables are adjusted automatically and the values that give
   *    the least feasible calculated value of the objective function are
   *    available in X(.) on the return from GETMIN.
   * ACC is a tolerance on the first order conditions at the calculated
   *    solution of the optimization problem.  These first order conditions
   *    state that, if X(.) is a solution, then there is a set of active
   *    constraints with indices IACT(K) K=1(1)NACT, say, such that X(.) is
   *    on the boundaries of these constraints, and the gradient of the
   *    objective function can be expressed in the form
   *         GRAD(F)=PAR(1)*GRAD(C(IACT(1)))+...
   *                       ...+PAR(NACT)*GRAD(C(IACT(NACT)))  .
   *    Here PAR(K) K=1(1)NACT are Lagrange multipliers that are nonpositive
   *    for inequality constraints, and GRAD(C(IACT(K))) is the gradient of
   *    the IACT(K)-th constraint function, so it is A(.,IACT(K)) if IACT(K)
   *    .LE. M, and it is minus or plus the J-th coordinate vector if the
   *    constraint is the lower or upper bound on X(J) respectively.  The
   *    normal return from the calculation occurs when X(.) is feasible and
   *    the sum of squares of components of the vector RESKT(.) is at most
   *    ACC**2, where RESKT(.) is the N-component vector of residuals of
   *    the first order condition that is displayed above.  Sometimes the
   *    package cannot satisfy this condition, because noise in the function
   *    values can prevent a change to the variables, no line search being
   *    allowed to increase the objective function.
   * IACT(.) is a working space array of integer variables that must be
   *    provided by the user.  Its length must be at least (M+2*N).  Its
   *    leading entries on the return from the subroutine are the indices
   *    IACT(K) K=1(1)NACT that are mentioned in the previous paragraph:
   *    in other words they are the indices of the final active constraints.
   *    Here the indices M+1,...,M+N and M+N+1,...,M+2*N denote the lower
   *    and upper bounds respectively.
   * NACT is set automatically to the integer variable of this ilk that has
   *    been introduced already.
   * PAR(.) is a one-dimensional array that will hold the Lagrange
   *    multipliers PAR(K) K=1(1)NACT on the return from GETMIN, these
   *    parameters being defined in the above paragraph on ACC.  The length
   *    of PAR should be at least N.
   * IPRINT must be set by the user to specify the frequency of printing
   *    during the execution of the optimization package.  There is no
   *    printed output if IPRINT=0.  Otherwise, after ensuring feasibility,
   *    information is given every IABS(IPRINT) iterations and whenever a
   *    parameter called TOL is reduced.  The printing provides the values
   *    of X(.), F(.) and G(.)=GRAD(F) if IPRINT is positive, while if
   *    IPRINT is negative this information is augmented by the current
   *    values of IACT(K) K=1(1)NACT, PAR(K) K=1(1)NACT and RESKT(I)
   *    I=1(1)N.  The reason for returning to the calling program is also
   *    displayed when IPRINT is nonzero.
   * INFO is an integer variable that should be set to zero initially,
   *    unless the user wishes to impose an upper bound on the number of
   *    evaluations of the objective function and its gradient, in which
   *    case INFO should be set to the value of this bound.  On the exit
   *    from GETMIN it will have one of the following integer values to
   *    indicate the reason for leaving the optimization package:
   *         INFO=1   X(.) is feasible and the condition that depends on
   *    ACC is satisfied.
   *         INFO=2   X(.) is feasible and rounding errors are preventing
   *    further progress.
   *         INFO=3   X(.) is feasible but the objective function fails to
   *    decrease although a decrease is predicted by the current gradient
   *    vector.  If this return occurs and RESKT(.) has large components
   *    then the user's calculation of the gradient of the objective
   *    function may be incorrect.  One should also question the coding of
   *    the gradient when the final rate of convergence is slow.
   *         INFO=4   In this case the calculation cannot begin because IA
   *    is less than N or because the lower bound on a variable is greater
   *    than the upper bound.
   *         INFO=5   This value indicates that the equality constraints
   *    are inconsistent.   These constraints include any components of
   *    X(.) that are frozen by setting XL(I)=XU(I).
   *         INFO=6   In this case there is an error return because the
   *    equality constraints and the bounds on the variables are found to
   *    be inconsistent.
   *         INFO=7   This value indicates that there is no vector of
   *    variables that satisfies all of the constraints.  Specifically,
   *    when this return or an INFO=6 return occurs, the current active
   *    constraints (whose indices are IACT(K) K=1(1)NACT) prevent any
   *    change in X(.) that reduces the sum of constraint violations,
   *    where only bounds are included in this sum if INFO=6.
   *         INFO=8   In this case the limit on the number of calls of
   *    subroutine FGCALC (see below) has been reached, and there would
   *    have been further calculation otherwise.
   * W(.) is a working space array of real variables that must be provided
   *    by the user.  Its length must be at least (M+11*N+N**2).  On exit
   *    from the package one can find the final components of GRAD(F) and
   *    RESKT(.) in W(1),...,W(N) and W(N+1),...,W(2*N) respectively.
   * Note 1.   The variables N, M, MEQ, IA, ACC and IPRINT and the elements
   *    of the arrays A(,.,), B(.), XL(.) and XU(.) are not altered by the
   *    optimization procedure.  Their values, the value of INFO and the
   *    initial components of X(.) must be set on entry to GETMIN.
   * Note 2.   Of course the package needs the objective function and its
   *    gradient.  Therefore the user must provide a subroutine called
   *    FGCALC, its first two lines being
   *         SUBROUTINE FGCALC (N,X,F,G)
   *         DIMENSION X(*),G(*)   .
   *    It is called automatically with N set as above and with X(.) set
   *    to a feasible vector of variables.  It should calculate the value
   *    of the objective function and its gradient for this X(.) and should
   *    set them in F and G(I) I=1(1)N respectively, without disturbing N
   *    or any of the components of X(.).
   */
  extern void getmin_(const INTEGER* n, const INTEGER* m, const INTEGER* meq,
                      const REAL a[], const INTEGER* ia, const REAL b[],
                      const REAL xl[], const REAL xu[], REAL x[],
                      const REAL* acc, INTEGER iact[], INTEGER* nact,
                      REAL par[], const INTEGER* iprint, INTEGER* info,
                      REAL w[]);

#if 0

  /* For the record, prototypes of the FORTRAN subroutines. */

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

#endif

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _TOLMIN_H */
