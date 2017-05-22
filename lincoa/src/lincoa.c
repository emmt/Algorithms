/*
 * lincoa.c --
 *
 * Implementation of Mike Powell's LINCOA algorithm for minimizing a function
 * of many variables subject to linear constraints.  The algorithm is described
 * in:
 *
 *     M.J.D. Powell, "On fast trust region methods for quadratic models with
 *     linear constraints", Report of the Department of Applied Mathematics and
 *     Theoretical Physics, Cambridge University, DAMTP 2014/NA02 (2014).
 *
 * The present code is based on the original FORTRAN version written by Mike
 * Powell who released his code under the GNU Lesser General Public License.
 * His original code is available at CCPForge
 * <https://ccpforge.cse.rl.ac.uk/gf/project/powell/>
 *
 * ----------------------------------------------------------------------------
 *
 * Copyright (c) 2014, Mike J. D. Powell (FORTRAN version released under the
 * GNU Lesser General Public License).
 *
 * Copyright (c) 2017, Éric Thiébaut (C version).
 *
 * ----------------------------------------------------------------------------
 *
 * History of changes:
 *  - Convert original FORTRAN code to C.
 *  - Make all possible arguments passed by address `const`.
 *  - Constant scalar arguments passed by address now passed by value.
 *  - Remove all array adjustments and use macros to mimic FORTRAN indexing.
 *  - Replace FORTRAN intrinsics by fast macros.
 *  - Replace all FORTRAN i/o by calls to the standard C library (the outputs
 *    with IPRINT=-1 are identical between the original FORTRAN code and the C
 *    version for the provided example).
 *  - Change workspace layout to accomodate to different sizes for integers and
 *    reals.
 *  - Pass objective function as argument.
 *
 * Things to do:
 *  - Use reverse communication.
 *  - Define status constants and related messages.
 *  - Write FORTRAN wrapper with the same syntax as the original code.
 *  - Cleanup goto statements.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "lincoa.h"

/* Basic types. */
#define LOGICAL int     /* type for boolean value */
#define INTEGER int     /* used for indexing */
#define REAL    double  /* floating-point type */

/* Where to print the results. */
#define OUTPUT stdout


/* Set some constants. */
static const REAL ZERO  = 0.0;
static const REAL HALF  = 0.5;
static const REAL ONE   = 1.0;
static const REAL PI    = M_PI;
static const REAL TENTH = 0.1;
static const REAL TINY  = 1e-60;


/* Fast macros to emulate FORTRAN intrinsics. */
#define sign(x,y)         _sign(typeof(x), x, y)
#define abs(x)            _abs(typeof(x), x)
#define pow2(x)           _pow2(typeof(x), x)
#define min(a,b)          _min(typeof((a) + (b)), a, b)
#define max(a,b)          _max(typeof((a) + (b)), a, b)
#define clamp(x,lo,hi)    _clamp(typeof((x) + (lo) + (hi)), x, lo, hi)
#define howmany(a,b)      _howmany(typeof((a) + (b)), a, b)
#define roundup(a,b)      _roundup(typeof((a) + (b)), a, b)
#define _abs(T,x)         ({ T _x = (x);  _x >= (T)0 ? _x : -_x; })
#define _pow2(T,x)        ({ T _x = (x);  _x*_x; })
#define _min(T,a,b)       ({ T _a = (a), _b = (b); _a < _b ? _a : _b; })
#define _max(T,a,b)       ({ T _a = (a), _b = (b); _a > _b ? _a : _b; })
#define _clamp(T,x,lo,hi) _min(T, _max(T, x, lo), hi)
#define _sign(T,x,y)      ({ T _x = (x); (_x < 0) != ((y) < 0) ? -_x : _x; })
#define _howmany(T,a,b)   ({ T _a = (a), _b = (b); ((_a + _b - 1)/_b); })
#define _roundup(T,a,b)   ({ T _a = (a), _b = (b); ((_a + _b - 1)/_b)*_b; })


/* A macro to print vectors (given their expression). */
#define PRTVECT(out, idx, beg, end, expr, fmt, ncols)                   \
  do {                                                                  \
    fputs("  ", out);                                                   \
    for (idx = (beg); idx <= (end); ++idx) {                            \
      if ((ncols) > 0 && idx > (beg) && (idx - (beg)) % (ncols) == 0) { \
        fputs("\n  ", out);                                             \
      }                                                                 \
      fprintf(out, fmt, (double)(expr));                                \
    }                                                                   \
    fputs("\n", out);                                                   \
  } while (0)


/* Define macros to mimic FORTRAN indexing. */
#define A(i1,i2)    a[(i1)-1 + ia*((i2)-1)]
#define AMAT(i1,i2) amat[(i1)-1 + n*((i2)-1)]
#define B(i1)       b[(i1)-1]
#define BMAT(i1,i2) bmat[(i1)-1 + ndim*((i2)-1)]
#define D(i1)       d[(i1)-1]
#define DW(i1)      dw[(i1)-1]
#define FVAL(i1)    fval[(i1)-1]
#define G(i1)       g[(i1)-1]
#define GL(i1)      gl[(i1)-1]
#define GOPT(i1)    gopt[(i1)-1]
#define HQ(i1)      hq[(i1)-1]
#define IACT(i1)    iact[(i1)-1]
#define PQ(i1)      pq[(i1)-1]
#define PQW(i1)     pqw[(i1)-1]
#define QFAC(i1,i2) qfac[(i1)-1 + n*((i2)-1)]
#define RESACT(i1)  resact[(i1)-1]
#define RESCON(i1)  rescon[(i1)-1]
#define RESNEW(i1)  resnew[(i1)-1]
#define RFAC(i1)    rfac[(i1)-1]
#define RSTAT(i1)   rstat[(i1)-1]
#define SP(i1)      sp[(i1)-1]
#define STEP(i1)    step[(i1)-1]
#define VLAG(i1)    vlag[(i1)-1]
#define VLAM(i1)    vlam[(i1)-1]
#define W(i1)       w[(i1)-1]
#define X(i1)       x[(i1)-1]
#define XP(i1)      xp[(i1)-1]
#define XBASE(i1)   xbase[(i1)-1]
#define XNEW(i1)    xnew[(i1)-1]
#define XOPT(i1)    xopt[(i1)-1]
#define XPT(i1,i2)  xpt[(i1)-1 + npt*((i2)-1)]
#define XSAV(i1)    xsav[(i1)-1]
#define YP(i1)      yp[(i1)-1]
#define ZP(i1)      zp[(i1)-1]
#define ZMAT(i1,i2) zmat[(i1)-1 + npt*((i2)-1)]
#define WB(i1)      wb[(i1)-1]


/*
 * The arguments N, NPT, XPT, BMAT, ZMAT, IDZ, NDIM ,SP and STEP are
 *   identical to the corresponding arguments in SUBROUTINE LINCOB.
 * KOPT is such that XPT(KOPT,.) is the current trust region centre.
 * KNEW on exit is usually positive, and then it is the index of an
 *   interpolation point to be moved to the position XPT(KOPT,.)+STEP(.).
 *   It is set on entry either to its final value or to 0. In the latter
 *   case, the final value of KNEW is chosen to maximize the denominator
 *   of the matrix updating formula times a weighting factor.
 * VLAG and W are used for working space, the first NPT+N elements of
 *   both of these vectors being required.
 *
 * The arrays BMAT and ZMAT with IDZ are updated, the new matrices being
 *   the ones that are suitable after the shift of the KNEW-th point to
 *   the new position XPT(KOPT,.)+STEP(.). A return with KNEW set to zero
 *   occurs if the calculation fails due to a zero denominator in the
 *   updating formula, which should never happen.
 */
static void
update(const INTEGER n,
       const INTEGER npt,
       const REAL xpt[],
       REAL bmat[],
       REAL zmat[],
       INTEGER *idz,             /* inp/out */
       const INTEGER ndim,
       const REAL sp[],
       const REAL step[],
       const INTEGER kopt,
       INTEGER *knew,            /* inp/out */
       REAL vlag[],
       REAL w[])
{
  REAL dx, tau, sum, ssq, beta, temp, bsum, hdiag, scala, scalb,
    alpha, denom, tempa, tempb, tausq, denabs, denmax, distsq, sqrtdn;
  INTEGER i, j, k, ja, jb, jl, jp, nptm;
  int iflag;

  /* Set some constants. */
  nptm = npt - n - 1;

  /* Calculate VLAG and BETA for the current choice of STEP. The first NPT
   *   elements of VLAG are set to the values of the Lagrange functions at
   *   XPT(KOPT,.)+STEP(.). The first NPT components of W_check are held
   *   in W, where W_check is defined in a paper on the updating method. */
  for (k = 1; k <= npt; ++k) {
    W(k) = SP(npt+k)*(HALF*SP(npt+k) + SP(k));
    sum = ZERO;
    for (j = 1; j <= n; ++j) {
      sum += BMAT(k,j)*STEP(j);
    }
    VLAG(k) = sum;
  }
  beta = ZERO;
  for (k = 1; k <= nptm; ++k) {
    sum = ZERO;
    for (i = 1; i <= npt; ++i) {
      sum += ZMAT(i,k)*W(i);
    }
    if (k < *idz) {
      beta += sum*sum;
      sum = -sum;
    } else {
      beta -= sum*sum;
    }
    for (i = 1; i <= npt; ++i) {
      VLAG(i) += sum*ZMAT(i,k);
    }
  }
  bsum = ZERO;
  dx = ZERO;
  ssq = ZERO;
  for (j = 1; j <= n; ++j) {
    sum = ZERO;
    for (i = 1; i <= npt; ++i) {
      sum += W(i)*BMAT(i,j);
    }
    bsum += sum*STEP(j);
    jp = npt + j;
    for (k = 1; k <= n; ++k) {
      sum += BMAT(jp,k)*STEP(k);
    }
    VLAG(jp) = sum;
    bsum += sum*STEP(j);
    dx += STEP(j)*XPT(kopt,j);
    ssq += pow2(STEP(j));
  }
  beta = dx*dx + ssq*(SP(kopt) + dx + dx + HALF*ssq) + beta - bsum;
  VLAG(kopt) += ONE;

  /* If KNEW is zero initially, then pick the index of the interpolation
   *   point to be deleted, by maximizing the absolute value of the
   *   denominator of the updating formula times a weighting factor. */
  if (*knew == 0) {
    denmax = ZERO;
    for (k = 1; k <= npt; ++k) {
      hdiag = ZERO;
      for (j = 1; j <= nptm; ++j) {
        temp = ONE;
        if (j < *idz) {
          temp = -ONE;
        }
        hdiag += temp*pow2(ZMAT(k,j));
      }
      denabs = abs(beta*hdiag + pow2(VLAG(k)));
      distsq = ZERO;
      for (j = 1; j <= n; ++j) {
        distsq += pow2(XPT(k,j) - XPT(kopt,j));
      }
      temp = denabs*distsq*distsq;
      if (temp > denmax) {
        denmax = temp;
        *knew = k;
      }
    }
  }

  /* Apply the rotations that put zeros in the KNEW-th row of ZMAT. */
  jl = 1;
  if (nptm >= 2) {
    for (j = 2; j <= nptm; ++j) {
      if (j == *idz) {
        jl = *idz;
      } else if (ZMAT(*knew,j) != ZERO) {
        temp = sqrt(pow2(ZMAT(*knew,jl)) + pow2(ZMAT(*knew,j)));
        tempa = ZMAT(*knew,jl)/temp;
        tempb = ZMAT(*knew,j)/temp;
        for (i = 1; i <= npt; ++i) {
          temp = tempa*ZMAT(i,jl) + tempb*ZMAT(i,j);
          ZMAT(i,j) = tempa*ZMAT(i,j) - tempb*ZMAT(i,jl);
          ZMAT(i,jl) = temp;
        }
        ZMAT(*knew,j) = ZERO;
      }
    }
  }

  /* Put the first NPT components of the KNEW-th column of the Z Z^T matrix
   *   into W, and calculate the parameters of the updating formula. */
  tempa = ZMAT(*knew,1);
  if (*idz >= 2) {
    tempa = -tempa;
  }
  if (jl > 1) {
    tempb = ZMAT(*knew,jl);
  }
  for (i = 1; i <= npt; ++i) {
    W(i) = tempa*ZMAT(i,1);
    if (jl > 1) {
      W(i) += tempb*ZMAT(i,jl);
    }
  }
  alpha = W(*knew);
  tau = VLAG(*knew);
  tausq = tau*tau;
  denom = alpha*beta + tausq;
  VLAG(*knew) -= ONE;
  if (denom == ZERO) {
    *knew = 0;
    return;
  }
  sqrtdn = sqrt((abs(denom)));

  /* Complete the updating of ZMAT when there is only one nonzero element
   *   in the KNEW-th row of the new matrix ZMAT. IFLAG is set to one when
   *   the value of IDZ is going to be reduced. */
  iflag = 0;
  if (jl == 1) {
    tempa = tau/sqrtdn;
    tempb = ZMAT(*knew,1)/sqrtdn;
    for (i = 1; i <= npt; ++i) {
      ZMAT(i,1) = tempa*ZMAT(i,1) - tempb*VLAG(i);
    }
    if (denom < ZERO) {
      if (*idz == 1) {
        *idz = 2;
      } else {
        iflag = 1;
      }
    }
  } else {
    /* Complete the updating of ZMAT in the alternative case. */
    ja = 1;
    if (beta >= ZERO) {
      ja = jl;
    }
    jb = jl + 1 - ja;
    temp = ZMAT(*knew,jb)/denom;
    tempa = temp*beta;
    tempb = temp*tau;
    temp = ZMAT(*knew,ja);
    scala = ONE/sqrt(abs(beta)*temp*temp + tausq);
    scalb = scala*sqrtdn;
    for (i = 1; i <= npt; ++i) {
      ZMAT(i,ja) = scala*(tau*ZMAT(i,ja) - temp*VLAG(i));
      ZMAT(i,jb) = scalb*(ZMAT(i,jb) - tempa*W(i) - tempb*VLAG(i));
    }
    if (denom <= ZERO) {
      if (beta < ZERO) {
        ++(*idz);
      } else {
        iflag = 1;
      }
    }
  }

  /* Reduce IDZ when the diagonal part of the ZMAT times Diag(DZ) times
   *   ZMAT^T factorization gains another positive element. Then exchange
   *   the first and IDZ-th columns of ZMAT. */
  if (iflag == 1) {
    --(*idz);
    for (i = 1; i <= npt; ++i) {
      temp = ZMAT(i,1);
      ZMAT(i,1) = ZMAT(i,*idz);
      ZMAT(i,*idz) = temp;
    }
  }

  /* Finally, update the matrix BMAT. */
  for (j = 1; j <= n; ++j) {
    jp = npt + j;
    W(jp) = BMAT(*knew,j);
    tempa = (alpha*VLAG(jp) - tau*W(jp))/denom;
    tempb = (-beta*W(jp) - tau*VLAG(jp))/denom;
    for (i = 1; i <= jp; ++i) {
      BMAT(i,j) = BMAT(i,j) + tempa*VLAG(i) + tempb*W(i);
      if (i > npt) {
        BMAT(jp,i-npt) = BMAT(i,j);
      }
    }
  }
}


/*
 * The arguments N, NPT, M, AMAT, B, X, RHOBEG, IPRINT, XBASE, XPT, FVAL,
 *   XSAV, XOPT, GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SP and RESCON are the
 *   same as the corresponding arguments in SUBROUTINE LINCOB.
 * KOPT is set to the integer such that XPT(KOPT,.) is the initial trust
 *   region centre.
 * IDZ is going to be set to one, so that every element of Diag(DZ) is
 *   one in the product ZMAT times Diag(DZ) times ZMAT^T, which is the
 *   factorization of the leading NPT by NPT submatrix of H.
 * STEP, PQW and W are used for working space, the arrays STEP and PQW
 *   being taken from LINCOB. The length of W must be at least N+NPT.
 *
 * SUBROUTINE PRELIM provides the elements of XBASE, XPT, BMAT and ZMAT
 *   for the first iteration, an important feature being that, if any of
 *   of the columns of XPT is an infeasible point, then the largest of
 *   the constraint violations there is at least 0.2*RHOBEG. It also sets
 *   the initial elements of FVAL, XOPT, GOPT, HQ, PQ, SP and RESCON.
 */
static void
prelim(const INTEGER n,
       const INTEGER npt,
       const INTEGER m,
       lincoa_objfun* objfun,
       void* data,
       const REAL amat[],
       REAL b[],
       REAL x[],
       const REAL rhobeg,
       const INTEGER iprint,
       REAL xbase[],
       REAL xpt[],
       REAL fval[],
       REAL xsav[],
       REAL xopt[],
       REAL gopt[],
       INTEGER *kopt,
       REAL hq[],
       REAL pq[],
       REAL bmat[],
       REAL zmat[],
       INTEGER *idz,
       const INTEGER ndim,
       REAL sp[],
       REAL rescon[],
       REAL step[],
       REAL pqw[],
       REAL w[])
{
  REAL f;
  REAL feas, bigv;
  REAL recip, reciq, resid;
  REAL rhosq;
  REAL temp;
  REAL test;
  INTEGER i, j, k, nf, jp;
  INTEGER ipt, jpt;
  INTEGER itemp;
  INTEGER jsav;
  INTEGER kbase;
  INTEGER nptm, nhq;

  /* Set some constants. */
  nptm = npt - n - 1;
  nhq = (n*n + n)/2;
  rhosq = rhobeg*rhobeg;
  recip = ONE/rhosq;
  reciq = sqrt(HALF)/rhosq;
  test = rhobeg*0.2;
  *idz = 1;
  kbase = 1;

  /* FIXME: Fix uninitialized variables. */
  jsav = 0;

  /* Set the initial elements of XPT, BMAT, SP and ZMAT to zero. */
  for (j = 1; j <= n; ++j) {
    XBASE(j) = X(j);
    for (k = 1; k <= npt; ++k) {
      XPT(k,j) = ZERO;
    }
    for (i = 1; i <= ndim; ++i) {
      BMAT(i,j) = ZERO;
    }
  }
  for (k = 1; k <= npt; ++k) {
    SP(k) = ZERO;
    for (j = 1; j <= nptm; ++j) {
      ZMAT(k,j) = ZERO;
    }
  }

  /* Set the nonzero coordinates of XPT(K,.), K=1,2,...,min[2*N+1,NPT],
   *   but they may be altered later to make a constraint violation
   *   sufficiently large. The initial nonzero elements of BMAT and of
   *   the first min[N,NPT-N-1] columns of ZMAT are set also. */
  for (j = 1; j <= n; ++j) {
    XPT(j+1,j) = rhobeg;
    if (j < npt - n) {
      jp = n + j + 1;
      XPT(jp,j) = -(rhobeg);
      BMAT(j+1,j) = HALF/rhobeg;
      BMAT(jp,j) = -HALF/rhobeg;
      ZMAT(1,j) = -reciq - reciq;
      ZMAT(j+1,j) = reciq;
      ZMAT(jp,j) = reciq;
    } else {
      BMAT(1,j) = -ONE/rhobeg;
      BMAT(j+1,j) = ONE/rhobeg;
      BMAT(npt+j,j) = -HALF*rhosq;
    }
  }

  /* Set the remaining initial nonzero elements of XPT and ZMAT when the
   *   number of interpolation points exceeds 2*N+1. */
  if (npt > 2*n + 1) {
    for (k = n + 1; k <= nptm; ++k) {
      itemp = (k - 1)/n;
      ipt = k - itemp*n;
      jpt = ipt + itemp;
      if (jpt > n) {
        jpt -= n;
      }
      XPT(n+k+1,ipt) = rhobeg;
      XPT(n+k+1,jpt) = rhobeg;
      ZMAT(1,k) = recip;
      ZMAT(ipt+1,k) = -recip;
      ZMAT(jpt+1,k) = -recip;
      ZMAT(n+k+1,k) = recip;
    }
  }

  /* Update the constraint right hand sides to allow for the shift XBASE. */
  if (m > 0) {
    for (j = 1; j <= m; ++j) {
      temp = ZERO;
      for (i = 1; i <= n; ++i) {
        temp += AMAT(i,j)*XBASE(i);
      }
      B(j) -= temp;
    }
  }

  /* Go through the initial points, shifting every infeasible point if
   *   necessary so that its constraint violation is at least 0.2*RHOBEG. */
  for (nf = 1; nf <= npt; ++nf) {
    feas = ONE;
    bigv = ZERO;
    j = 0;
  L80:
    ++j;
    if (j <= m && nf >= 2) {
      resid = -B(j);
      for (i = 1; i <= n; ++i) {
        resid += XPT(nf,i)*AMAT(i,j);
      }
      if (resid <= bigv) {
        goto L80;
      }
      bigv = resid;
      jsav = j;
      if (resid <= test) {
        feas = -ONE;
        goto L80;
      }
      feas = ZERO;
    }
    if (feas < ZERO) {
      for (i = 1; i <= n; ++i) {
        STEP(i) = XPT(nf,i) + (test - bigv)*AMAT(i,jsav);
      }
      for (k = 1; k <= npt; ++k) {
        SP(npt+k) = ZERO;
        for (j = 1; j <= n; ++j) {
          SP(npt+k) += XPT(k,j)*STEP(j);
        }
      }
      update(n, npt, xpt, bmat, zmat, idz, ndim, sp, step,
             kbase, &nf, pqw, w);
      for (i = 1; i <= n; ++i) {
        XPT(nf,i) = STEP(i);
      }
    }

    /* Calculate the objective function at the current interpolation point,
     *   and set KOPT to the index of the first trust region centre. */
    for (j = 1; j <= n; ++j) {
      X(j) = XBASE(j) + XPT(nf,j);
    }
    f = objfun(n, x, (feas > ZERO), data);
    if (iprint == 3) {
      fprintf(OUTPUT,
              "    Function number %5ld"
              "    F =%18.10E"
              "    The corresponding X is:\n",
              (long)nf, (double)f);
      PRTVECT(OUTPUT, i, 1, n, X(i), "%15.6E", 5);
    }
    if (nf == 1) {
      *kopt = 1;
    } else if (f < FVAL(*kopt) && feas > ZERO) {
      *kopt = nf;
    }
    FVAL(nf) = f;
  }

  /* Set PQ for the first quadratic model. */
  for (j = 1; j <= nptm; ++j) {
    W(j) = ZERO;
    for (k = 1; k <= npt; ++k) {
      W(j) += ZMAT(k,j)*FVAL(k);
    }
  }
  for (k = 1; k <= npt; ++k) {
    PQ(k) = ZERO;
    for (j = 1; j <= nptm; ++j) {
      PQ(k) += ZMAT(k,j)*W(j);
    }
  }

  /* Set XOPT, SP, GOPT and HQ for the first quadratic model. */
  for (j = 1; j <= n; ++j) {
    XOPT(j) = XPT(*kopt,j);
    XSAV(j) = XBASE(j) + XOPT(j);
    GOPT(j) = ZERO;
  }
  for (k = 1; k <= npt; ++k) {
    SP(k) = ZERO;
    for (j = 1; j <= n; ++j) {
      SP(k) += XPT(k,j)*XOPT(j);
    }
    temp = PQ(k)*SP(k);
    for (j = 1; j <= n; ++j) {
      GOPT(j) = GOPT(j) + FVAL(k)*BMAT(k,j) + temp*XPT(k,j);
    }
  }
  for (i = 1; i <= nhq; ++i) {
    HQ(i) = ZERO;
  }

  /* Set the initial elements of RESCON. */
  for (j = 1; j <= m; ++j) {
    temp = B(j);
    for (i = 1; i <= n; ++i) {
      temp -= XOPT(i)*AMAT(i,j);
    }
    temp = max(temp,ZERO);
    if (temp >= rhobeg) {
      temp = -temp;
    }
    RESCON(j) = temp;
  }
}


/*
 * N, M, AMAT, B, NACT, IACT, QFAC and RFAC are the same as the terms
 *   with these names in SUBROUTINE LINCOB. The current values must be
 *   set on entry. NACT, IACT, QFAC and RFAC are kept up to date when
 *   GETACT changes the current active set.
 * SNORM, RESNEW, RESACT, G and DW are the same as the terms with these
 *   names in SUBROUTINE TRSTEP. The elements of RESNEW and RESACT are
 *   also kept up to date.
 * VLAM and W are used for working space, the vector VLAM being reserved
 *   for the Lagrange multipliers of the calculation. Their lengths must
 *   be at least N.
 * The main purpose of GETACT is to pick the current active set. It is
 *   defined by the property that the projection of -G into the space
 *   orthogonal to the active constraint normals is as large as possible,
 *   subject to this projected steepest descent direction moving no closer
 *   to the boundary of every constraint whose current residual is at most
 *   0.2*SNORM. On return, the settings in NACT, IACT, QFAC and RFAC are
 *   all appropriate to this choice of active set.
 * Occasionally this projected direction is zero, and then the final value
 *   of W(1) is set to zero. Otherwise, the direction itself is returned
 *   in DW, and W(1) is set to the square of the length of the direction.
 */
static void
getact(const INTEGER n,
       const INTEGER m,
       const REAL amat[],
       const REAL b[],
       INTEGER *nact,
       INTEGER iact[],
       REAL qfac[],
       REAL rfac[],
       const REAL *snorm,
       REAL resnew[],
       REAL resact[],
       const REAL g[],
       REAL dw[],
       REAL vlam[],
       REAL w[])
{
  REAL dd;
  REAL dnorm, sprod, vmult, violmx;
  REAL rdiag, ddsav;
  REAL sum, cval, tdel, ctol, temp, sval, cosv, test, sinv;
  INTEGER i, j, k, l;
  INTEGER ic, jc, jw, jcp;
  INTEGER idiag, jdiag, iflag;
  INTEGER nactp, numb;

  /* Set some constants and a temporary VLAM. */
  tdel = (*snorm)*0.2;
  ddsav = ZERO;
  for (i = 1; i <= n; ++i) {
    ddsav += pow2(G(i));
    VLAM(i) = ZERO;
  }
  ddsav += ddsav;

  /* Set the initial QFAC to the identity matrix in the case NACT=0. */
  if (*nact == 0) {
    for (i = 1; i <= n; ++i) {
      for (j = 1; j <= n; ++j) {
        QFAC(i,j) = ZERO;
      }
      QFAC(i,i) = ONE;
    }
    goto L100;
  }

  /* Remove any constraints from the initial active set whose residuals
   *   exceed TDEL. */
  iflag = 1;
  ic = *nact;
 L40:
  if (RESACT(ic) > tdel) {
    goto L800;
  }
 L50:
  --ic;
  if (ic > 0) {
    goto L40;
  }

  /* Remove any constraints from the initial active set whose Lagrange
   *   multipliers are nonnegative, and set the surviving multipliers. */
  iflag = 2;
 L60:
  if (*nact == 0) {
    goto L100;
  }
  ic = *nact;
 L70:
  temp = ZERO;
  for (i = 1; i <= n; ++i) {
    temp += QFAC(i,ic)*G(i);
  }
  idiag = (ic*ic + ic)/2;
  if (ic < *nact) {
    jw = idiag + ic;
    numb = *nact;
    for (j = ic + 1; j <= numb; ++j) {
      temp -= RFAC(jw)*VLAM(j);
      jw += j;
    }
  }
  if (temp >= ZERO) {
    goto L800;
  }
  VLAM(ic) = temp/RFAC(idiag);
  --ic;
  if (ic > 0) {
    goto L70;
  }

  /* Set the new search direction D. Terminate if the 2-norm of D is zero
   *   or does not decrease, or if NACT=N holds. The situation NACT=N
   *   occurs for sufficiently large SNORM if the origin is in the convex
   *   hull of the constraint gradients. */
 L100:
  if (*nact == n) {
    goto L290;
  }
  for (j = (*nact) + 1; j <= n; ++j) {
    W(j) = ZERO;
    for (i = 1; i <= n; ++i) {
      W(j) += QFAC(i,j)*G(i);
    }
  }
  dd = ZERO;
  for (i = 1; i <= n; ++i) {
    DW(i) = ZERO;
    for (j = (*nact) + 1; j <= n; ++j) {
      DW(i) -= W(j)*QFAC(i,j);
    }
    dd += pow2(DW(i));
  }
  if (dd >= ddsav) {
    goto L290;
  }
  if (dd == ZERO) {
    goto L300;
  }
  ddsav = dd;
  dnorm = sqrt(dd);

  /* Pick the next integer L or terminate, a positive value of L being
   *   the index of the most violated constraint. The purpose of CTOL
   *   below is to estimate whether a positive value of VIOLMX may be
   *   due to computer rounding errors. */
  l = 0;
  if (m > 0) {
    test = dnorm/(*snorm);
    violmx = ZERO;
    for (j = 1; j <= m; ++j) {
      if (RESNEW(j) > ZERO && RESNEW(j) <= tdel) {
        sum = ZERO;
        for (i = 1; i <= n; ++i) {
          sum += AMAT(i,j)*DW(i);
        }
        if (sum > test*RESNEW(j)) {
          if (sum > violmx) {
            l = j;
            violmx = sum;
          }
        }
      }
    }
    ctol = ZERO;
    temp = dnorm*0.01;
    if (violmx > ZERO && violmx < temp) {
      if (*nact > 0) {
        numb = *nact;
        for (k = 1; k <= numb; ++k) {
          j = IACT(k);
          sum = ZERO;
          for (i = 1; i <= n; ++i) {
            sum += DW(i)*AMAT(i,j);
          }
          ctol = max(ctol, abs(sum));
        }
      }
    }
  }
  W(1) = ONE;
  if (l == 0) {
    goto L300;
  }
  if (violmx <= ctol*10.0) {
    goto L300;
  }

  /* Apply Givens rotations to the last (N-NACT) columns of QFAC so that
   *   the first (NACT+1) columns of QFAC are the ones required for the
   *   addition of the L-th constraint, and add the appropriate column
   *   to RFAC. */
  nactp = (*nact) + 1;
  idiag = (nactp*nactp - nactp)/2;
  rdiag = ZERO;
  for (j = n; j >= 1; --j) {
    sprod = ZERO;
    for (i = 1; i <= n; ++i) {
      sprod += QFAC(i,j)*AMAT(i, l);
    }
    if (j <= *nact) {
      RFAC(idiag+j) = sprod;
    } else {
      if (abs(rdiag) <= abs(sprod)*1e-20) {
        rdiag = sprod;
      } else {
        temp = sqrt(sprod*sprod + rdiag*rdiag);
        cosv = sprod/temp;
        sinv = rdiag/temp;
        rdiag = temp;
        for (i = 1; i <= n; ++i) {
          temp = cosv*QFAC(i,j) + sinv*QFAC(i,j+1);
          QFAC(i,j+1) = -sinv*QFAC(i,j) + cosv*QFAC(i,j+1);
          QFAC(i,j) = temp;
        }
      }
    }
  }
  if (rdiag < ZERO) {
    for (i = 1; i <= n; ++i) {
      QFAC(i,nactp) = -QFAC(i,nactp);
    }
  }
  RFAC(idiag+nactp) = abs(rdiag);
  *nact = nactp;
  IACT(*nact) = l;
  RESACT(*nact) = RESNEW(l);
  VLAM(*nact) = ZERO;
  RESNEW(l) = ZERO;

  /* Set the components of the vector VMU in W. */
 L220:
  W(*nact) = ONE/pow2(RFAC(((*nact)*(*nact)+(*nact))/2));
  if (*nact > 1) {
    for (i = (*nact) - 1; i >= 1; --i) {
      idiag = (i*i + i)/2;
      jw = idiag + i;
      sum = ZERO;
      numb = *nact;
      for (j = i + 1; j <= numb; ++j) {
        sum -= RFAC(jw)*W(j);
        jw += j;
      }
      W(i) = sum/RFAC(idiag);
    }
  }

  /* Calculate the multiple of VMU to subtract from VLAM, and update VLAM. */
  vmult = violmx;
  ic = 0;
  j = 1;
 L250:
  if (j < *nact) {
    if (VLAM(j) >= vmult*W(j)) {
      ic = j;
      vmult = VLAM(j)/W(j);
    }
    ++j;
    goto L250;
  }
  numb = *nact;
  for (j = 1; j <= numb; ++j) {
    VLAM(j) -= vmult*W(j);
  }
  if (ic > 0) {
    VLAM(ic) = ZERO;
  }
  violmx = max(violmx - vmult, ZERO);
  if (ic == 0) {
    violmx = ZERO;
  }

  /* Reduce the active set if necessary, so that all components of the
   *   new VLAM are negative, with resetting of the residuals of the
   *   constraints that become inactive. */
  iflag = 3;
  ic = *nact;
 L270:
  if (VLAM(ic) < ZERO) {
    goto L280;
  }
  RESNEW(IACT(ic)) = max(RESACT(ic), TINY);
  goto L800;
 L280:
  --ic;
  if (ic > 0) {
    goto L270;
  }

  /* Calculate the next VMU if VIOLMX is positive. Return if NACT=N holds,
   *   as then the active constraints imply D=0. Otherwise, go to label
   *   100, to calculate the new D and to test for termination. */
  if (violmx > ZERO) {
    goto L220;
  }
  if (*nact < n) {
    goto L100;
  }
 L290:
  dd = ZERO;
 L300:
  W(1) = dd;
  return;

  /* These instructions rearrange the active constraints so that the new
   *   value of IACT(NACT) is the old value of IACT(IC). A sequence of
   *   Givens rotations is applied to the current QFAC and RFAC. Then NACT
   *   is reduced by one. */
 L800:
  RESNEW(IACT(ic)) = max(RESACT(ic), TINY);
  jc = ic;
 L810:
  if (jc < *nact) {
    jcp = jc + 1;
    idiag = jc*jcp/2;
    jw = idiag + jcp;
    temp = sqrt(pow2(RFAC(jw-1)) + pow2(RFAC(jw)));
    cval = RFAC(jw)/temp;
    sval = RFAC(jw-1)/temp;
    RFAC(jw-1) = sval*RFAC(idiag);
    RFAC(jw) = cval*RFAC(idiag);
    RFAC(idiag) = temp;
    if (jcp < *nact) {
      numb = *nact;
      for (j = jcp + 1; j <= numb; ++j) {
        temp = sval*RFAC(jw+jc) + cval*RFAC(jw+jcp);
        RFAC(jw+jcp) = cval*RFAC(jw+jc) - sval*RFAC(jw+jcp);
        RFAC(jw+jc) = temp;
        jw += j;
      }
    }
    jdiag = idiag - jc;
    for (i = 1; i <= n; ++i) {
      if (i < jc) {
        temp = RFAC(idiag+i);
        RFAC(idiag+i) = RFAC(jdiag+i);
        RFAC(jdiag+i) = temp;
      }
      temp = sval*QFAC(i,jc) + cval*QFAC(i,jcp);
      QFAC(i,jcp) = cval*QFAC(i,jc) - sval*QFAC(i,jcp);
      QFAC(i,jc) = temp;
    }
    IACT(jc) = IACT(jcp);
    RESACT(jc) = RESACT(jcp);
    VLAM(jc) = VLAM(jcp);
    jc = jcp;
    goto L810;
  }
  --(*nact);
  switch (iflag) {
  case 1:  goto L50;
  case 2:  goto L60;
  case 3:  goto L280;
  }
}


/*
 * N, NPT, M, AMAT, B, XPT, HQ, PQ, NACT, IACT, RESCON, QFAC and RFAC
 *   are the same as the terms with these names in LINCOB. If RESCON(J)
 *   is negative, then |RESCON(J)| must be no less than the trust region
 *   radius, so that the J-th constraint can be ignored.
 * SNORM is set to the trust region radius DELTA initially. On the
 *   return, however, it is the length of the calculated STEP, which is
 *   set to zero if the constraints do not allow a long enough step.
 * STEP is the total calculated step so far from the trust region centre,
 *   its final value being given by the sequence of CG iterations, which
 *   terminate if the trust region boundary is reached.
 * G must be set on entry to the gradient of the quadratic model at the
 *   trust region centre. It is used as working space, however, and is
 *   always the gradient of the model at the current STEP, except that
 *   on return the value of G(1) is set to ONE instead of to ZERO if
 *   and only if GETACT is called more than once.
 * RESNEW, RESACT, D, DW and W are used for working space. A negative
 *   value of RESNEW(J) indicates that the J-th constraint does not
 *   restrict the CG steps of the current trust region calculation, a
 *   zero value of RESNEW(J) indicates that the J-th constraint is active,
 *   and otherwise RESNEW(J) is set to the greater of TINY and the actual
 *   residual of the J-th constraint for the current STEP. RESACT holds
 *   the residuals of the active constraints, which may be positive.
 *   D is the search direction of each line search. DW is either another
 *   search direction or the change in gradient along D. The length of W
 *   must be at least MAX[M,2*N].
 */
static void
trstep(const INTEGER n,
       const INTEGER npt,
       const INTEGER m,
       const REAL amat[],
       const REAL b[],
       const REAL xpt[],
       const REAL hq[],
       const REAL pq[],
       INTEGER *nact,
       INTEGER iact[],
       const REAL rescon[],
       REAL qfac[],
       REAL rfac[],
       REAL *snorm,
       REAL step[],
       REAL g[],
       REAL resnew[],
       REAL resact[],
       REAL d[],
       REAL dw[],
       REAL w[])
{
  REAL ad, dd, dg;
  REAL alphm, alpht, ctest;
  REAL ds;
  REAL reduct, resmax;
  REAL ss, dgd, adw, wgd, rhs, sum, beta;
  REAL temp, snsq, gamma, alpbd, alpha, scale;
  INTEGER i, j, k;
  INTEGER icount;
  INTEGER ih;
  INTEGER ir;
  INTEGER jsav;
  INTEGER ncall, n1, n2;

  /* Set some numbers for the conjugate gradient iterations. */
  ctest = 0.01;
  snsq = (*snorm)*(*snorm);

  /* Set the initial elements of RESNEW, RESACT and STEP. */
  if (m > 0) {
    for (j = 1; j <= m; ++j) {
      RESNEW(j) = RESCON(j);
      if (RESCON(j) >= *snorm) {
        RESNEW(j) = -ONE;
      } else if (RESCON(j) >= ZERO) {
        RESNEW(j) = max(RESNEW(j), TINY);
      }
    }
    if (*nact > 0) {
      n1 = *nact;
      for (k = 1; k <= n1; ++k) {
        RESACT(k) = RESCON(IACT(k));
        RESNEW(IACT(k)) = ZERO;
      }
    }
  }
  for (i = 1; i <= n; ++i) {
    STEP(i) = ZERO;
  }
  ss = ZERO;
  reduct = ZERO;
  ncall = 0;

  /* GETACT picks the active set for the current STEP. It also sets DW to
   *   the vector closest to -G that is orthogonal to the normals of the
   *   active constraints. DW is scaled to have length 0.2*SNORM, as then
   *   a move of DW from STEP is allowed by the linear constraints. */
 L40:
  ++ncall;
  getact(n, m, amat, b, nact, iact, qfac, rfac, snorm, resnew,
         resact, g, dw, w, w + n);
  if (W(n+1) == ZERO) {
    goto L320;
  }
  scale = (*snorm)*0.2/sqrt(W(n+1));
  for (i = 1; i <= n; ++i) {
    DW(i) = scale*DW(i);
  }

  /* If the modulus of the residual of an active constraint is substantial,
   *   then set D to the shortest move from STEP to the boundaries of the
   *   active constraints. */
  resmax = ZERO;
  if (*nact > 0) {
    n1 = *nact;
    for (k = 1; k <= n1; ++k) {
      resmax = max(resmax, RESACT(k));
    }
  }
  gamma = ZERO;
  if (resmax > (*snorm)*1e-4) {
    ir = 0;
    n1 = *nact;
    for (k = 1; k <= n1; ++k) {
      temp = RESACT(k);
      if (k >= 2) {
        n2 = k - 1;
        for (i = 1; i <= n2; ++i) {
          ++ir;
          temp -= RFAC(ir)*W(i);
        }
      }
      ++ir;
      W(k) = temp/RFAC(ir);
    }
    for (i = 1; i <= n; ++i) {
      D(i) = ZERO;
      n2 = *nact;
      for (k = 1; k <= n2; ++k) {
        D(i) += W(k)*QFAC(i,k);
      }
    }

    /* The vector D that has just been calculated is also the shortest move
     *   from STEP+DW to the boundaries of the active constraints. Set GAMMA
     *   to the greatest steplength of this move that satisfies the trust
     *   region bound. */
    rhs = snsq;
    ds = ZERO;
    dd = ZERO;
    for (i = 1; i <= n; ++i) {
      sum = STEP(i) + DW(i);
      rhs -= sum*sum;
      ds += D(i)*sum;
      dd += pow2(D(i));
    }
    if (rhs > ZERO) {
      temp = sqrt(ds*ds + dd*rhs);
      if (ds <= ZERO) {
        gamma = (temp - ds)/dd;
      } else {
        gamma = rhs/(temp + ds);
      }
    }

    /* Reduce the steplength GAMMA if necessary so that the move along D
     *   also satisfies the linear constraints. */
    j = 0;
  L110:
    if (gamma > ZERO) {
      ++j;
      if (RESNEW(j) > ZERO) {
        ad = ZERO;
        adw = ZERO;
        for (i = 1; i <= n; ++i) {
          ad += AMAT(i,j)*D(i);
          adw += AMAT(i,j)*DW(i);
        }
        if (ad > ZERO) {
          temp = max((RESNEW(j) - adw)/ad, ZERO);
          gamma = min(gamma,temp);
        }
      }
      if (j < m) {
        goto L110;
      }
    }
    gamma = min(gamma,ONE);
  }

  /* Set the next direction for seeking a reduction in the model function
   *   subject to the trust region bound and the linear constraints. */
  if (gamma <= ZERO) {
    for (i = 1; i <= n; ++i) {
      D(i) = DW(i);
    }
    icount = *nact;
  } else {
    for (i = 1; i <= n; ++i) {
      D(i) = DW(i) + gamma*D(i);
    }
    icount = (*nact) - 1;
  }
  alpbd = ONE;

  /* Set ALPHA to the steplength from STEP along D to the trust region
   *   boundary. Return if the first derivative term of this step is
   *   sufficiently small or if no further progress is possible. */
 L150:
  ++icount;
  rhs = snsq - ss;
  if (rhs <= ZERO) {
    goto L320;
  }
  dg = ZERO;
  ds = ZERO;
  dd = ZERO;
  for (i = 1; i <= n; ++i) {
    dg += D(i)*G(i);
    ds += D(i)*STEP(i);
    dd += pow2(D(i));
  }
  if (dg >= ZERO) {
    goto L320;
  }
  temp = sqrt(rhs*dd + ds*ds);
  if (ds <= ZERO) {
    alpha = (temp - ds)/dd;
  } else {
    alpha = rhs/(temp + ds);
  }
  if (-alpha*dg <= ctest*reduct) {
    goto L320;
  }

  /* Set DW to the change in gradient along D. */
  ih = 0;
  for (j = 1; j <= n; ++j) {
    DW(j) = ZERO;
    for (i = 1; i <= j; ++i) {
      ++ih;
      if (i < j) {
        DW(j) += HQ(ih)*D(i);
      }
      DW(i) += HQ(ih)*D(j);
    }
  }
  for (k = 1; k <= npt; ++k) {
    temp = ZERO;
    for (j = 1; j <= n; ++j) {
      temp += XPT(k,j)*D(j);
    }
    temp = PQ(k)*temp;
    for (i = 1; i <= n; ++i) {
      DW(i) += temp*XPT(k,i);
    }
  }

  /* Set DGD to the curvature of the model along D. Then reduce ALPHA if
   *   necessary to the value that minimizes the model. */
  dgd = ZERO;
  for (i = 1; i <= n; ++i) {
    dgd += D(i)*DW(i);
  }
  alpht = alpha;
  if (dg + alpha*dgd > ZERO) {
    alpha = -dg/dgd;
  }

  /* Make a further reduction in ALPHA if necessary to preserve feasibility,
   *   and put some scalar products of D with constraint gradients in W. */
  alphm = alpha;
  jsav = 0;
  if (m > 0) {
    for (j = 1; j <= m; ++j) {
      ad = ZERO;
      if (RESNEW(j) > ZERO) {
        for (i = 1; i <= n; ++i) {
          ad += AMAT(i,j)*D(i);
        }
        if (alpha*ad > RESNEW(j)) {
          alpha = RESNEW(j)/ad;
          jsav = j;
        }
      }
      W(j) = ad;
    }
  }
  alpha = max(alpha,alpbd);
  alpha = min(alpha,alphm);
  if (icount == *nact) {
    alpha = min(alpha,ONE);
  }

  /* Update STEP, G, RESNEW, RESACT and REDUCT. */
  ss = ZERO;
  for (i = 1; i <= n; ++i) {
    STEP(i) += alpha*D(i);
    ss += pow2(STEP(i));
    G(i) += alpha*DW(i);
  }
  if (m > 0) {
    for (j = 1; j <= m; ++j) {
      if (RESNEW(j) > ZERO) {
        RESNEW(j) = max(RESNEW(j) - alpha*W(j), TINY);
      }
    }
  }
  if (icount == *nact && *nact > 0) {
    n2 = *nact;
    for (k = 1; k <= n2; ++k) {
      RESACT(k) = (ONE - gamma)*RESACT(k);
    }
  }
  reduct -= alpha*(dg + HALF*alpha*dgd);

  /* Test for termination. Branch to label 40 if there is a new active
   *   constraint and if the distance from STEP to the trust region
   *   boundary is at least 0.2*SNORM. */
  if (alpha == alpht) {
    goto L320;
  }
  temp = -alphm*(dg + HALF*alphm*dgd);
  if (temp <= ctest*reduct) {
    goto L320;
  }
  if (jsav > 0) {
    if (ss <= snsq*0.64) {
      goto L40;
    }
    goto L320;
  }
  if (icount == n) {
    goto L320;
  }

  /* Calculate the next search direction, which is conjugate to the
   *   previous one except in the case ICOUNT=NACT. */
  if (*nact > 0) {
    for (j = (*nact) + 1; j <= n; ++j) {
      W(j) = ZERO;
      for (i = 1; i <= n; ++i) {
        W(j) += G(i)*QFAC(i,j);
      }
    }
    for (i = 1; i <= n; ++i) {
      temp = ZERO;
      for (j = (*nact) + 1; j <= n; ++j) {
        temp += QFAC(i,j)*W(j);
      }
      W(n+i) = temp;
    }
  } else {
    for (i = 1; i <= n; ++i) {
      W(n+i) = G(i);
    }
  }
  if (icount == *nact) {
    beta = ZERO;
  } else {
    wgd = ZERO;
    for (i = 1; i <= n; ++i) {
      wgd += W(n+i)*DW(i);
    }
    beta = wgd/dgd;
  }
  for (i = 1; i <= n; ++i) {
    D(i) = -W(n+i) + beta*D(i);
  }
  alpbd = ZERO;
  goto L150;

  /* Return from the subroutine. */
 L320:
  *snorm = ZERO;
  if (reduct > ZERO) {
    *snorm = sqrt(ss);
  }
  G(1) = ZERO;
  if (ncall > 1) {
    G(1) = ONE;
  }
}


/*
 * N, NPT, M, AMAT, B, XPT, XOPT, NACT, IACT, RESCON, QFAC, KOPT are the
 *   same as the terms with these names in SUBROUTINE LINCOB.
 * KNEW is the index of the interpolation point that is going to be moved.
 * DEL is the current restriction on the length of STEP, which is never
 *   greater than the current trust region radius DELTA.
 * STEP will be set to the required step from XOPT to the new point.
 * GL must be set on entry to the gradient of LFUNC at XBASE, where LFUNC
 *   is the KNEW-th Lagrange function. It is used also for some other
 *   gradients of LFUNC.
 * PQW provides the second derivative parameters of LFUNC.
 * RSTAT and W are used for working space. Their lengths must be at least
 *   M and N, respectively. RSTAT(J) is set to -1.0, 0.0, or 1.0 if the
 *   J-th constraint is irrelevant, active, or both inactive and relevant,
 *   respectively.
 * IFEAS will be set to 0 or 1 if XOPT+STEP is infeasible or feasible.
 *
 * STEP is chosen to provide a relatively large value of the modulus of
 *   LFUNC(XOPT+STEP), subject to ||STEP|| .LE. DEL. A projected STEP is
 *   calculated too, within the trust region, that does not alter the
 *   residuals of the active constraints. The projected step is preferred
 *   if its value of | LFUNC(XOPT+STEP) | is at least one fifth of the
 *   original one, but the greatest violation of a linear constraint must
 *   be at least 0.2*DEL, in order to keep the interpolation points apart.
 *   The remedy when the maximum constraint violation is too small is to
 *   restore the original step, which is perturbed if necessary so that
 *   its maximum constraint violation becomes 0.2*DEL.
 */
static void
qmstep(const INTEGER n,
       const INTEGER npt,
       const INTEGER m,
       const REAL amat[],
       const REAL b[],
       const REAL xpt[],
       const REAL xopt[],
       const INTEGER nact,
       const INTEGER iact[],
       const REAL rescon[],
       const REAL qfac[],
       const INTEGER kopt,
       const INTEGER knew,
       const REAL del,
       REAL step[],
       REAL gl[],
       const REAL pqw[],
       REAL rstat[],
       REAL w[],
       INTEGER *ifeas)
{
  REAL gg, sp, ss, ww, ghg, sum, stp, vbig, bigv, vlag, ctol;
  REAL temp;
  REAL test, vnew;
  REAL vgrad, resmax, stpsav;
  INTEGER i, j, k;
  INTEGER jsav;
  INTEGER ksav;
  INTEGER n1, n2;

  /* FIXME: Fix uninitialized variables. */
  jsav = 0;
  ksav = 0;
  stpsav = 0;

  /* Set some constants. */
  test = del*0.2;

  /* Replace GL by the gradient of LFUNC at the trust region centre, and
   *   set the elements of RSTAT. */
  for (k = 1; k <= npt; ++k) {
    temp = ZERO;
    for (j = 1; j <= n; ++j) {
      temp += XPT(k,j)*XOPT(j);
    }
    temp = PQW(k)*temp;
    for (i = 1; i <= n; ++i) {
      GL(i) += temp*XPT(k,i);
    }
  }
  if (m > 0) {
    for (j = 1; j <= m; ++j) {
      RSTAT(j) = ONE;
      if (abs(RESCON(j)) >= del) {
        RSTAT(j) = -ONE;
      }
    }
    n2 = nact;
    for (k = 1; k <= n2; ++k) {
      RSTAT(IACT(k)) = ZERO;
    }
  }

  /* Find the greatest modulus of LFUNC on a line through XOPT and
   *   another interpolation point within the trust region. */
  vbig = ZERO;
  for (k = 1; k <= npt; ++k) {
    if (k == kopt) {
      goto L60;
    }
    ss = ZERO;
    sp = ZERO;
    for (i = 1; i <= n; ++i) {
      temp = XPT(k,i) - XOPT(i);
      ss += temp*temp;
      sp += GL(i)*temp;
    }
    stp = -del/sqrt(ss);
    if (k == knew) {
      if (sp*(sp - ONE) < ZERO) {
        stp = -stp;
      }
      vlag = abs(stp*sp) + stp*stp*abs(sp - ONE);
    } else {
      vlag = abs(stp*(ONE - stp)*sp);
    }
    if (vlag > vbig) {
      ksav = k;
      stpsav = stp;
      vbig = vlag;
    }
  L60:
    ;
  }

  /* Set STEP to the move that gives the greatest modulus calculated above.
   *   This move may be replaced by a steepest ascent step from XOPT. */
  gg = ZERO;
  for (i = 1; i <= n; ++i) {
    gg += pow2(GL(i));
    STEP(i) = stpsav*(XPT(ksav,i) - XOPT(i));
  }
  vgrad = del*sqrt(gg);
  if (vgrad <= TENTH*vbig) {
    goto L220;
  }

  /* Make the replacement if it provides a larger value of VBIG. */
  ghg = ZERO;
  for (k = 1; k <= npt; ++k) {
    temp = ZERO;
    for (j = 1; j <= n; ++j) {
      temp += XPT(k,j)*GL(j);
    }
    ghg += PQW(k)*temp*temp;
  }
  vnew = vgrad + abs(HALF*del*del*ghg/gg);
  if (vnew > vbig) {
    vbig = vnew;
    stp = del/sqrt(gg);
    if (ghg < ZERO) {
      stp = -stp;
    }
    for (i = 1; i <= n; ++i) {
      STEP(i) = stp*GL(i);
    }
  }
  if (nact == 0 || nact == n) {
    goto L220;
  }

  /* Overwrite GL by its projection. Then set VNEW to the greatest
   *   value of |LFUNC| on the projected gradient from XOPT subject to
   *   the trust region bound. If VNEW is sufficiently large, then STEP
   *   may be changed to a move along the projected gradient. */
  for (k = nact + 1; k <= n; ++k) {
    W(k) = ZERO;
    for (i = 1; i <= n; ++i) {
      W(k) += GL(i)*QFAC(i,k);
    }
  }
  gg = ZERO;
  for (i = 1; i <= n; ++i) {
    GL(i) = ZERO;
    for (k = nact + 1; k <= n; ++k) {
      GL(i) += QFAC(i,k)*W(k);
    }
    gg += pow2(GL(i));
  }
  vgrad = del*sqrt(gg);
  if (vgrad <= TENTH*vbig) {
    goto L220;
  }
  ghg = ZERO;
  for (k = 1; k <= npt; ++k) {
    temp = ZERO;
    for (j = 1; j <= n; ++j) {
      temp += XPT(k,j)*GL(j);
    }
    ghg += PQW(k)*temp*temp;
  }
  vnew = vgrad + abs(HALF*del*del*ghg/gg);

  /* Set W to the possible move along the projected gradient. */
  stp = del/sqrt(gg);
  if (ghg < ZERO) {
    stp = -stp;
  }
  ww = ZERO;
  for (i = 1; i <= n; ++i) {
    W(i) = stp*GL(i);
    ww += pow2(W(i));
  }

  /* Set STEP to W if W gives a sufficiently large value of the modulus
   *   of the Lagrange function, and if W either preserves feasibility
   *   or gives a constraint violation of at least 0.2*DEL. The purpose
   *   of CTOL below is to provide a check on feasibility that includes
   *   a tolerance for contributions from computer rounding errors. */
  if (vnew/vbig >= 0.2) {
    *ifeas = 1;
    bigv = ZERO;
    j = 0;
  L170:
    ++j;
    if (j <= m) {
      if (RSTAT(j) == ONE) {
        temp = -RESCON(j);
        for (i = 1; i <= n; ++i) {
          temp += W(i)*AMAT(i,j);
        }
        bigv = max(bigv,temp);
      }
      if (bigv < test) {
        goto L170;
      }
      *ifeas = 0;
    }
    ctol = ZERO;
    temp = sqrt(ww)*0.01;
    if (bigv > ZERO && bigv < temp) {
      n1 = nact;
      for (k = 1; k <= n1; ++k) {
        j = IACT(k);
        sum = ZERO;
        for (i = 1; i <= n; ++i) {
          sum += W(i)*AMAT(i,j);
        }
        ctol = max(ctol, abs(sum));
      }
    }
    if (bigv <= ctol*10.0 || bigv >= test) {
      for (i = 1; i <= n; ++i) {
        STEP(i) = W(i);
      }
      goto L260;
    }
  }

  /* Calculate the greatest constraint violation at XOPT+STEP with STEP at
   *   its original value. Modify STEP if this violation is unacceptable. */
 L220:
  *ifeas = 1;
  bigv = ZERO;
  resmax = ZERO;
  j = 0;
 L230:
  ++j;
  if (j <= m) {
    if (RSTAT(j) < ZERO) {
      goto L230;
    }
    temp = -RESCON(j);
    for (i = 1; i <= n; ++i) {
      temp += STEP(i)*AMAT(i,j);
    }
    resmax = max(resmax,temp);
    if (temp < test) {
      if (temp <= bigv) {
        goto L230;
      }
      bigv = temp;
      jsav = j;
      *ifeas = -1;
      goto L230;
    }
    *ifeas = 0;
  }
  if (*ifeas == -1) {
    for (i = 1; i <= n; ++i) {
      STEP(i) += (test - bigv)*AMAT(i,jsav);
    }
    *ifeas = 0;
  }

  /* Return the calculated STEP and the value of IFEAS. */
 L260:
  return;
}


/*
 * The arguments N, NPT, M, X, RHOBEG, RHOEND, IPRINT and MAXFUN are
 *   identical to the corresponding arguments in SUBROUTINE LINCOA.
 * AMAT is a matrix whose columns are the constraint gradients, scaled
 *   so that they have unit length.
 * B contains on entry the right hand sides of the constraints, scaled
 *   as above, but later B is modified for variables relative to XBASE.
 * XBASE holds a shift of origin that should reduce the contributions
 *   from rounding errors to values of the model and Lagrange functions.
 * XPT contains the interpolation point coordinates relative to XBASE.
 * FVAL holds the values of F at the interpolation points.
 * XSAV holds the best feasible vector of variables so far, without any
 *   shift of origin.
 * XOPT is set to XSAV-XBASE, which is the displacement from XBASE of
 *   the feasible vector of variables that provides the least calculated
 *   F so far, this vector being the current trust region centre.
 * GOPT holds the gradient of the quadratic model at XSAV = XBASE+XOPT.
 * HQ holds the explicit second derivatives of the quadratic model.
 * PQ contains the parameters of the implicit second derivatives of the
 *   quadratic model.
 * BMAT holds the last N columns of the big inverse matrix H.
 * ZMAT holds the factorization of the leading NPT by NPT submatrix
 *   of H, this factorization being ZMAT times Diag(DZ) times ZMAT^T,
 *   where the elements of DZ are plus or minus one, as specified by IDZ.
 * NDIM is the first dimension of BMAT and has the value NPT+N.
 * STEP is employed for trial steps from XOPT. It is also used for working
 *   space when XBASE is shifted and in PRELIM.
 * SP is reserved for the scalar products XOPT^T XPT(K,.), K=1,2,...,NPT,
 *   followed by STEP^T XPT(K,.), K=1,2,...,NPT.
 * XNEW is the displacement from XBASE of the vector of variables for
 *   the current calculation of F, except that SUBROUTINE TRSTEP uses it
 *   for working space.
 * IACT is an integer array for the indices of the active constraints.
 * RESCON holds useful information about the constraint residuals. Every
 *   nonnegative RESCON(J) is the residual of the J-th constraint at the
 *   current trust region centre. Otherwise, if RESCON(J) is negative, the
 *   J-th constraint holds as a strict inequality at the trust region
 *   centre, its residual being at least |RESCON(J)|; further, the value
 *   of |RESCON(J)| is at least the current trust region radius DELTA.
 * QFAC is the orthogonal part of the QR factorization of the matrix of
 *   active constraint gradients, these gradients being ordered in
 *   accordance with IACT. When NACT is less than N, columns are added
 *   to QFAC to complete an N by N orthogonal matrix, which is important
 *   for keeping calculated steps sufficiently close to the boundaries
 *   of the active constraints.
 * RFAC is the upper triangular part of this QR factorization, beginning
 *   with the first diagonal element, followed by the two elements in the
 *   upper triangular part of the second column and so on.
 * PQW is used for working space, mainly for storing second derivative
 *   coefficients of quadratic functions. Its length is NPT+N.
 * The array W is also used for working space. The required number of
 *   elements, namely MAX[M+3*N,2*M+N,2*NPT], is set in LINCOA.
 */
static void
lincob(const INTEGER n,
       const INTEGER npt,
       const INTEGER m,
       lincoa_objfun* objfun,
       void* data,
       const REAL amat[],
       REAL b[],
       REAL x[],
       const REAL rhobeg,
       const REAL rhoend,
       const INTEGER iprint,
       const INTEGER maxfun,
       REAL xbase[],
       REAL xpt[],
       REAL fval[],
       REAL xsav[],
       REAL xopt[],
       REAL gopt[],
       REAL hq[],
       REAL pq[],
       REAL bmat[],
       REAL zmat[],
       const INTEGER ndim,
       REAL step[],
       REAL sp[],
       REAL xnew[],
       INTEGER iact[],
       REAL rescon[],
       REAL qfac[],
       REAL rfac[],
       REAL pqw[],
       REAL w[])
{
  REAL del, delsav, delta, xdiff, distsq, f, fsave, qoptsq, xoptsq, ratio,
    vquad, vqalt, rho, sum, ssq, diff, snorm, dffalt, temp, fopt, sumz;
  INTEGER i, j, k, ih, nf, nh, ip, np, idz, ifeas, itest, kopt, nptm, ksave,
    nact, knew, nvala, nvalb;


  /* FIXME: Fix uninitialized variables. */
  dffalt = 0;
  ratio = 0;
  snorm = 0;

  /* Set some constants. */
  np = n + 1;
  nh = n*np/2;
  nptm = npt - np;

  /* Set the elements of XBASE, XPT, FVAL, XSAV, XOPT, GOPT, HQ, PQ, BMAT,
   *   ZMAT and SP for the first iteration. An important feature is that,
   *   if the interpolation point XPT(K,.) is not feasible, where K is any
   *   integer from [1,NPT], then a change is made to XPT(K,.) if necessary
   *   so that the constraint violation is at least 0.2*RHOBEG. Also KOPT
   *   is set so that XPT(KOPT,.) is the initial trust region centre. */
  prelim(n, npt, m, objfun, data, amat, b, x, rhobeg, iprint, xbase, xpt,
         fval, xsav, xopt, gopt, &kopt, hq, pq, bmat, zmat, &idz, ndim,
         sp, rescon, step, pqw, w);

  /* Begin the iterative procedure. */
  nf = npt;
  f = fopt = FVAL(kopt);
  rho = rhobeg;
  delta = rho;
  ifeas = 0;
  nact = 0;
  itest = 3;
 L10:
  knew = 0;
  nvala = 0;
  nvalb = 0;

  /* Shift XBASE if XOPT may be too far from XBASE. First make the changes
   *   to BMAT that do not depend on ZMAT. */
 L20:
  fsave = fopt;
  xoptsq = ZERO;
  for (i = 1; i <= n; ++i) {
    xoptsq += pow2(XOPT(i));
  }
  if (xoptsq >= delta*1e4*delta) {
    qoptsq = xoptsq*0.25;
    for (k = 1; k <= npt; ++k) {
      sum = ZERO;
      for (i = 1; i <= n; ++i) {
        sum += XPT(k,i)*XOPT(i);
      }
      sum -= HALF*xoptsq;
      W(npt+k) = sum;
      SP(k) = ZERO;
      for (i = 1; i <= n; ++i) {
        XPT(k,i) = XPT(k,i) - HALF*XOPT(i);
        STEP(i) = BMAT(k,i);
        W(i) = sum*XPT(k,i) + qoptsq*XOPT(i);
        ip = npt + i;
        for (j = 1; j <= i; ++j) {
          BMAT(ip,j) = BMAT(ip,j) + STEP(i)*W(j) + W(i)*STEP(j);
        }
      }
    }

    /* Then the revisions of BMAT that depend on ZMAT are calculated. */
    for (k = 1; k <= nptm; ++k) {
      sumz = ZERO;
      for (i = 1; i <= npt; ++i) {
        sumz += ZMAT(i,k);
        W(i) = W(npt+i)*ZMAT(i,k);
      }
      for (j = 1; j <= n; ++j) {
        sum = qoptsq*sumz*XOPT(j);
        for (i = 1; i <= npt; ++i) {
          sum += W(i)*XPT(i,j);
        }
        STEP(j) = sum;
        if (k < idz) {
          sum = -sum;
        }
        for (i = 1; i <= npt; ++i) {
          BMAT(i,j) = BMAT(i,j) + sum*ZMAT(i,k);
        }
      }
      for (i = 1; i <= n; ++i) {
        ip = i + npt;
        temp = STEP(i);
        if (k < idz) {
          temp = -temp;
        }
        for (j = 1; j <= i; ++j) {
          BMAT(ip,j) = BMAT(ip,j) + temp*STEP(j);
        }
      }
    }

    /* Update the right hand sides of the constraints. */
    if (m > 0) {
      for (j = 1; j <= m; ++j) {
        temp = ZERO;
        for (i = 1; i <= n; ++i) {
          temp += AMAT(i,j)*XOPT(i);
        }
        B(j) -= temp;
      }
    }

    /* The following instructions complete the shift of XBASE, including the
     *   changes to the parameters of the quadratic model. */
    ih = 0;
    for (j = 1; j <= n; ++j) {
      W(j) = ZERO;
      for (k = 1; k <= npt; ++k) {
        W(j) += PQ(k)*XPT(k,j);
        XPT(k,j) = XPT(k,j) - HALF*XOPT(j);
      }
      for (i = 1; i <= j; ++i) {
        ++ih;
        HQ(ih) = HQ(ih) + W(i)*XOPT(j) + XOPT(i)*W(j);
        BMAT(npt+i,j) = BMAT(npt+j,i);
      }
    }
    for (j = 1; j <= n; ++j) {
      XBASE(j) += XOPT(j);
      XOPT(j) = ZERO;
      XPT(kopt,j) = ZERO;
    }
  }

  /* In the case KNEW=0, generate the next trust region step by calling
   *   TRSTEP, where SNORM is the current trust region radius initially.
   *   The final value of SNORM is the length of the calculated step,
   *   except that SNORM is zero on return if the projected gradient is
   *   unsuitable for starting the conjugate gradient iterations. */
  delsav = delta;
  ksave = knew;
  if (knew == 0) {
    snorm = delta;
    for (i = 1; i <= n; ++i) {
      XNEW(i) = GOPT(i);
    }
    trstep(n, npt, m, amat, b, xpt, hq, pq, &nact, iact, rescon, qfac,
           rfac, &snorm, step, xnew, w, &W(m+1), pqw, &PQW(np), &W(m+np));

    /* A trust region step is applied whenever its length, namely SNORM, is at
     *   least HALF*DELTA. It is also applied if its length is at least 0.1999
     *   times DELTA and if a line search of TRSTEP has caused a change to the
     *   active set. Otherwise there is a branch below to label 530 or 560. */
    temp = HALF*delta;
    if (XNEW(1) >= HALF) {
      temp = delta*0.1999;
    }
    if (snorm <= temp) {
      delta = HALF*delta;
      if (delta <= rho*1.4) {
        delta = rho;
      }
      ++nvala;
      ++nvalb;
      temp = snorm/rho;
      if (delsav > rho) {
        temp = ONE;
      }
      if (temp >= HALF) {
        nvala = 0;
      }
      if (temp >= TENTH) {
        nvalb = 0;
      }
      if (delsav > rho) {
        goto L530;
      }
      if (nvala < 5 && nvalb < 3) {
        goto L530;
      }
      if (snorm > ZERO) {
        ksave = -1;
      }
      goto L560;
    }
    nvala = 0;
    nvalb = 0;

    /* Alternatively, KNEW is positive. Then the model step is calculated
     *   within a trust region of radius DEL, after setting the gradient at
     *   XBASE and the second derivative parameters of the KNEW-th Lagrange
     *   function in W(1) to W(N) and in PQW(1) to PQW(NPT), respectively. */
  } else {
    del = max(TENTH*delta, rho);
    for (i = 1; i <= n; ++i) {
      W(i) = BMAT(knew,i);
    }
    for (k = 1; k <= npt; ++k) {
      PQW(k) = ZERO;
    }
    for (j = 1; j <= nptm; ++j) {
      temp = ZMAT(knew,j);
      if (j < idz) {
        temp = -temp;
      }
      for (k = 1; k <= npt; ++k) {
        PQW(k) += temp*ZMAT(k,j);
      }
    }
    qmstep(n, npt, m, amat, b, xpt, xopt, nact, iact, rescon, qfac,
           kopt, knew, del, step, w, pqw, &W(np), &W(np+m), &ifeas);
  }

  /* Set VQUAD to the change to the quadratic model when the move STEP is
   *   made from XOPT. If STEP is a trust region step, then VQUAD should be
   *   negative. If it is nonnegative due to rounding errors in this case,
   *   there is a branch to label 530 to try to improve the model. */
  vquad = ZERO;
  ih = 0;
  for (j = 1; j <= n; ++j) {
    vquad += STEP(j)*GOPT(j);
    for (i = 1; i <= j; ++i) {
      ++ih;
      temp = STEP(i)*STEP(j);
      if (i == j) {
        temp = HALF*temp;
      }
      vquad += temp*HQ(ih);
    }
  }
  for (k = 1; k <= npt; ++k) {
    temp = ZERO;
    for (j = 1; j <= n; ++j) {
      temp += XPT(k,j)*STEP(j);
      SP(npt+k) = temp;
    }
    vquad += HALF*PQ(k)*temp*temp;
  }
  if (ksave == 0 && vquad >= ZERO) {
    goto L530;
  }

  /* Calculate the next value of the objective function. The difference
   *   between the actual new value of F and the value predicted by the
   *   model is recorded in DIFF. */
 L220:
  ++nf;
  if (nf > maxfun) {
    --nf;
    if (iprint > 0) {
      fputs("Return from LINCOA because objective function has "
            "been called MAXFUN times.\n", stderr);
    }
    goto L600;
  }
  xdiff = ZERO;
  for (i = 1; i <= n; ++i) {
    XNEW(i) = XOPT(i) + STEP(i);
    X(i) = XBASE(i) + XNEW(i);
    xdiff += pow2(X(i) - XSAV(i));
  }
  xdiff = sqrt(xdiff);
  if (ksave == -1) {
    xdiff = rho;
  }
  if (xdiff <= TENTH*rho || xdiff >= delta + delta) {
    ifeas = 0;
    if (iprint > 0) {
      fputs("Return from LINCOA because rounding errors "
            "prevent reasonable changes to X.\n", stderr);
    }
    goto L600;
  }
  if (ksave <= 0) {
    ifeas = 1;
  }
  f = objfun(n, x, (ifeas > 0), data);
  if (iprint == 3) {
      fprintf(OUTPUT,
              "    Function number %5ld"
              "    F =%18.10E"
              "    The corresponding X is:\n",
              (long)nf, (double)f);
      PRTVECT(OUTPUT, i, 1, n, X(i), "%15.6E", 5);
  }
  if (ksave == -1) {
    goto L600;
  }
  diff = f - fopt - vquad;

  /* If X is feasible, then set DFFALT to the difference between the new
   *   value of F and the value predicted by the alternative model. */
  if (ifeas == 1 && itest < 3) {
    for (k = 1; k <= npt; ++k) {
      PQW(k) = ZERO;
      W(k) = FVAL(k) - FVAL(kopt);
    }
    for (j = 1; j <= nptm; ++j) {
      sum = ZERO;
      for (i = 1; i <= npt; ++i) {
        sum += W(i)*ZMAT(i,j);
      }
      if (j < idz) {
        sum = -sum;
      }
      for (k = 1; k <= npt; ++k) {
        PQW(k) += sum*ZMAT(k,j);
      }
    }
    vqalt = ZERO;
    for (k = 1; k <= npt; ++k) {
      sum = ZERO;
      for (j = 1; j <= n; ++j) {
        sum += BMAT(k,j)*STEP(j);
      }
      vqalt += sum*W(k);
      vqalt += PQW(k)*SP(npt+k)*(HALF*SP(npt+k) + SP(k));
    }
    dffalt = f - fopt - vqalt;
  }
  if (itest == 3) {
    dffalt = diff;
    itest = 0;
  }

  /* Pick the next value of DELTA after a trust region step. */
  if (ksave == 0) {
    ratio = (f - fopt)/vquad;
    if (ratio <= TENTH) {
      delta = HALF*delta;
    } else if (ratio <= 0.7) {
      delta = max(HALF*delta, snorm);
    } else {
      temp = sqrt(2.0)*delta;
      delta = max(HALF*delta, snorm + snorm);
      delta = min(delta, temp);
    }
    if (delta <= rho*1.4) {
      delta = rho;
    }
  }

  /* Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
   *   can be moved. If STEP is a trust region step, then KNEW is zero at
   *   present, but a positive value is picked by subroutine UPDATE. */
  update(n, npt, xpt, bmat, zmat, &idz, ndim, sp, step, kopt,
         &knew, pqw, w);
  if (knew == 0) {
    if (iprint > 0) {
      fputs("Return from LINCOA because the denominator of "
            "the updating formula is zero.\n", stderr);
    }
    goto L600;
  }

  /* If ITEST is increased to 3, then the next quadratic model is the
   *   one whose second derivative matrix is least subject to the new
   *   interpolation conditions. Otherwise the new model is constructed
   *   by the symmetric Broyden method in the usual way. */
  if (ifeas == 1) {
    ++itest;
    if (abs(dffalt) >= TENTH*abs(diff)) {
      itest = 0;
    }
  }

  /* Update the second derivatives of the model by the symmetric Broyden
   *   method, using PQW for the second derivative parameters of the new
   *   KNEW-th Lagrange function. The contribution from the old parameter
   *   PQ(KNEW) is included in the second derivative matrix HQ. W is used
   *   later for the gradient of the new KNEW-th Lagrange function. */
  if (itest < 3) {
    for (k = 1; k <= npt; ++k) {
      PQW(k) = ZERO;
    }
    for (j = 1; j <= nptm; ++j) {
      temp = ZMAT(knew,j);
      if (temp != ZERO) {
        if (j < idz) {
          temp = -temp;
        }
        for (k = 1; k <= npt; ++k) {
          PQW(k) += temp*ZMAT(k,j);
        }
      }
    }
    ih = 0;
    for (i = 1; i <= n; ++i) {
      W(i) = BMAT(knew,i);
      temp = PQ(knew)*XPT(knew,i);
      for (j = 1; j <= i; ++j) {
        ++ih;
        HQ(ih) += temp*XPT(knew,j);
      }
    }
    PQ(knew) = ZERO;
    for (k = 1; k <= npt; ++k) {
      PQ(k) += diff*PQW(k);
    }
  }

  /* Include the new interpolation point with the corresponding updates of
   *   SP. Also make the changes of the symmetric Broyden method to GOPT at
   *   the old XOPT if ITEST is less than 3. */
  FVAL(knew) = f;
  SP(knew) = SP(kopt) + SP(npt+kopt);
  ssq = ZERO;
  for (i = 1; i <= n; ++i) {
    XPT(knew,i) = XNEW(i);
    ssq += pow2(STEP(i));
  }
  SP(npt+knew) = SP(npt+kopt) + ssq;
  if (itest < 3) {
    for (k = 1; k <= npt; ++k) {
      temp = PQW(k)*SP(k);
      for (i = 1; i <= n; ++i) {
        W(i) += temp*XPT(k,i);
      }
    }
    for (i = 1; i <= n; ++i) {
      GOPT(i) += diff*W(i);
    }
  }

  /* Update FOPT, XSAV, XOPT, KOPT, RESCON and SP if the new F is the
   *   least calculated value so far with a feasible vector of variables. */
  if (f < fopt && ifeas == 1) {
    fopt = f;
    for (j = 1; j <= n; ++j) {
      XSAV(j) = X(j);
      XOPT(j) = XNEW(j);
    }
    kopt = knew;
    snorm = sqrt(ssq);
    for (j = 1; j <= m; ++j) {
      if (RESCON(j) >= delta + snorm) {
        RESCON(j) = snorm - RESCON(j);
      } else {
        RESCON(j) += snorm;
        if (RESCON(j) + delta > ZERO) {
          temp = B(j);
          for (i = 1; i <= n; ++i) {
            temp -= XOPT(i)*AMAT(i,j);
          }
          temp = max(temp,ZERO);
          if (temp >= delta) {
            temp = -temp;
          }
          RESCON(j) = temp;
        }
      }
    }
    for (k = 1; k <= npt; ++k) {
      SP(k) += SP(npt+k);
    }

    /* Also revise GOPT when symmetric Broyden updating is applied. */
    if (itest < 3) {
      ih = 0;
      for (j = 1; j <= n; ++j) {
        for (i = 1; i <= j; ++i) {
          ++ih;
          if (i < j) {
            GOPT(j) += HQ(ih)*STEP(i);
          }
          GOPT(i) += HQ(ih)*STEP(j);
        }
      }
      for (k = 1; k <= npt; ++k) {
        temp = PQ(k)*SP(npt+k);
        for (i = 1; i <= n; ++i) {
          GOPT(i) += temp*XPT(k,i);
        }
      }
    }
  }

  /* Replace the current model by the least Frobenius norm interpolant if
   *   this interpolant gives substantial reductions in the predictions
   *   of values of F at feasible points. */
  if (itest == 3) {
    for (k = 1; k <= npt; ++k) {
      PQ(k) = ZERO;
      W(k) = FVAL(k) - FVAL(kopt);
    }
    for (j = 1; j <= nptm; ++j) {
      sum = ZERO;
      for (i = 1; i <= npt; ++i) {
        sum += W(i)*ZMAT(i,j);
      }
      if (j < idz) {
        sum = -sum;
      }
      for (k = 1; k <= npt; ++k) {
        PQ(k) += sum*ZMAT(k,j);
      }
    }
    for (j = 1; j <= n; ++j) {
      temp = ZERO;
      for (i = 1; i <= npt; ++i) {
       temp += W(i)*BMAT(i,j);
      }
      GOPT(j) = temp;
    }
    for (k = 1; k <= npt; ++k) {
      temp = PQ(k)*SP(k);
      for (i = 1; i <= n; ++i) {
        GOPT(i) += temp*XPT(k,i);
      }
    }
    for (ih = 1; ih <= nh; ++ih) {
      HQ(ih) = ZERO;
    }
  }

  /* If a trust region step has provided a sufficient decrease in F, then
   *   branch for another trust region calculation. Every iteration that
   *   takes a model step is followed by an attempt to take a trust region
   *   step. */
  knew = 0;
  if (ksave > 0) {
    goto L20;
  }
  if (ratio >= TENTH) {
    goto L20;
  }

  /* Alternatively, find out if the interpolation points are close enough
   *   to the best point so far. */
 L530:
  distsq = max(delta*delta, rho*4.0*rho);
  for (k = 1; k <= npt; ++k) {
    sum = ZERO;
    for (j = 1; j <= n; ++j) {
      sum += pow2(XPT(k,j) - XOPT(j));
    }
    if (sum > distsq) {
      knew = k;
      distsq = sum;
    }
  }

  /* If KNEW is positive, then branch back for the next iteration, which
   *   will generate a "model step". Otherwise, if the current iteration
   *   has reduced F, or if DELTA was above its lower bound when the last
   *   trust region step was calculated, then try a "trust region" step
   *   instead. */
  if (knew > 0) {
    goto L20;
  }
  knew = 0;
  if (fopt < fsave) {
    goto L20;
  }
  if (delsav > rho) {
    goto L20;
  }

  /* The calculations with the current value of RHO are complete.
   *   Pick the next value of RHO. */
 L560:
  if (rho > rhoend) {
    delta = HALF*rho;
    if (rho > rhoend*250.0) {
      rho = TENTH*rho;
    } else if (rho <= rhoend*16.0) {
      rho = rhoend;
    } else {
      rho = sqrt(rho*rhoend);
    }
    delta = max(delta,rho);
    if (iprint >= 2) {
      if (iprint >= 3) {
        fputs("\n", OUTPUT);
      }
      fprintf(OUTPUT,
              "    New RHO =%11.4E "
              "    Number of function values =%6ld\n",
              (double)rho, (long)nf);
      fprintf(OUTPUT,
              "    Least value of F =%23.15E     "
              "    The corresponding X is:\n",
              (double)fopt);
      PRTVECT(OUTPUT, i, 1, n, XBASE(i) + XOPT(i), "%15.6E", 5);
    }
    goto L10;
  }

  /* Return from the calculation, after branching to label 220 for another
   *   Newton-Raphson step if it has not been tried before. */
  if (ksave == -1) {
    goto L220;
  }
 L600:
  if (fopt <= f || ifeas == 0) {
    for (i = 1; i <= n; ++i) {
      X(i) = XSAV(i);
    }
    f = fopt;
  }
  if (iprint >= 1) {
    fprintf(OUTPUT, "\n"
            "    At the return from LINCOA "
            "    Number of function values =%6ld\n", (long)nf);
    fprintf(OUTPUT,
            "    Least value of F =%23.15E     "
            "    The corresponding X is:\n", (double)fopt);
    PRTVECT(OUTPUT, i, 1, n, X(i), "%15.6E", 5);
  }
  W(1) = f;
  W(2) = (REAL) nf + HALF;
}


/* Yields the number of reals needed by the workspace. */
static INTEGER
realworkspace(INTEGER n, INTEGER npt, INTEGER m)
{
  return (max(max(m+3*n, 2*m+n), 2*npt) + m*(2 + n) +
          npt*(4 + n + npt) + n*(8 + 3*n));
}

size_t
lincoa_storage(INTEGER n, INTEGER npt, INTEGER m)
{
  size_t nr = realworkspace(n, npt, m);
  size_t ni = n;
  return roundup(nr*sizeof(REAL), sizeof(INTEGER)) + ni*sizeof(INTEGER);
}


/*
 * This subroutine seeks the least value of a function of many variables,
 *   subject to general linear inequality constraints, by a trust region
 *   method that forms quadratic models by interpolation. Usually there
 *   is much freedom in each new model after satisfying the interpolation
 *   conditions, which is taken up by minimizing the Frobenius norm of
 *   the change to the second derivative matrix of the model. One new
 *   function value is calculated on each iteration, usually at a point
 *   where the current model predicts a reduction in the least value so
 *   far of the objective function subject to the linear constraints.
 *   Alternatively, a new vector of variables may be chosen to replace
 *   an interpolation point that may be too far away for reliability, and
 *   then the new point does not have to satisfy the linear constraints.
 *   The arguments of the subroutine are as follows.
 *
 * N must be set to the number of variables and must be at least two.
 * NPT must be set to the number of interpolation conditions, which is
 *   required to be in the interval [N+2,(N+1)(N+2)/2]. Typical choices
 *   of the author are NPT=N+6 and NPT=2*N+1. Larger values tend to be
 *   highly inefficent when the number of variables is substantial, due
 *   to the amount of work and extra difficulty of adjusting more points.
 * M must be set to the number of linear inequality constraints.
 * A is a matrix whose columns are the constraint gradients, which are
 *   required to be nonzero.
 * IA is the first dimension of the array A, which must be at least N.
 * B is the vector of right hand sides of the constraints, the J-th
 *   constraint being that the scalar product of A(.,J) with X(.) is at
 *   most B(J). The initial vector X(.) is made feasible by increasing
 *   the value of B(J) if necessary.
 * X is the vector of variables. Initial values of X(1),X(2),...,X(N)
 *   must be supplied. If they do not satisfy the constraints, then B
 *   is increased as mentioned above. X contains on return the variables
 *   that have given the least calculated F subject to the constraints.
 * RHOBEG and RHOEND must be set to the initial and final values of a
 *   trust region radius, so both must be positive with RHOEND<=RHOBEG.
 *   Typically, RHOBEG should be about one tenth of the greatest expected
 *   change to a variable, and RHOEND should indicate the accuracy that
 *   is required in the final values of the variables.
 * The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
 *   amount of printing. Specifically, there is no output if IPRINT=0 and
 *   there is output only at the return if IPRINT=1. Otherwise, the best
 *   feasible vector of variables so far and the corresponding value of
 *   the objective function are printed whenever RHO is reduced, where
 *   RHO is the current lower bound on the trust region radius. Further,
 *   each new value of F with its variables are output if IPRINT=3.
 * MAXFUN must be set to an upper bound on the number of calls of CALFUN,
 *   its value being at least NPT+1.
 * W is an array used for working space. Its length must be at least
 *   M*(2+N) + NPT*(4+N+NPT) + N*(9+3*N) + MAX[ M+3*N, 2*M+N, 2*NPT ].
 *   On return, W(1) is set to the final value of F, and W(2) is set to
 *   the total number of function evaluations plus 0.5.
 *
 * SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set
 *   F to the value of the objective function for the variables X(1),
 *   X(2),...,X(N). The value of the argument F is positive when CALFUN
 *   is called if and only if the current X satisfies the constraints
 *   to working accuracy.
 */
void
lincoa(const INTEGER n,
       const INTEGER npt,
       const INTEGER m,
       lincoa_objfun* objfun,
       void* data,
       const REAL a[],
       const INTEGER ia,
       const REAL b[],
       REAL x[],
       const REAL rhobeg,
       const REAL rhoend,
       const INTEGER iprint,
       const INTEGER maxfun,
       void* ws)
{
  REAL smallx, sum, temp;
  INTEGER i, j, np, iflag, nptm;

  /* Check that N, NPT, M and MAXFUN are acceptable. */
  smallx = rhoend*1e-6;
  np = n + 1;
  nptm = npt - np;
  if (n < 2) {
    fputs("Return from LINCOA because N is less than 2.\n", stderr);
    return;
  }
  if (npt < n + 2 || npt > (n + 2)*np/2) {
    fputs("Return from LINCOA because NPT is not in the required interval.\n",
          stderr);
    return;
  }
  if (m < 0) {
    fputs("Return from LINCOA because M is less than 0.\n", stderr);
    return;
  }
  if (maxfun <= npt) {
    fputs("Return from LINCOA because MAXFUN is less than NPT+1.\n",
          stderr);
    return;
  }


  {
    /* Partition the working space array, so that different parts of it
     *   can be treated separately by the subroutine that performs the
     *   main calculation. */
    INTEGER ndim = npt + n;
    size_t nr = realworkspace(n, npt, m);
    REAL* w = (REAL*)ws;
    REAL* amat = &w[max(max(m + 3*n, 2*m + n), 2*npt)];
    REAL* wb = &amat[m*n];
    REAL* xbase = &wb[m];
    REAL* xpt = &xbase[n];
    REAL* fval = &xpt[npt*n];
    REAL* xsav = &fval[npt];
    REAL* xopt = &xsav[n];
    REAL* gopt = &xopt[n];
    REAL* hq = &gopt[n];
    REAL* pq = &hq[n*np/2];
    REAL* bmat = &pq[npt];
    REAL* zmat = &bmat[ndim*n];
    REAL* step = &zmat[npt*nptm];
    REAL* sp = &step[n];
    REAL* xnew = &sp[npt + npt];
    REAL* rescon = &xnew[n];
    REAL* qfac = &rescon[m];
    REAL* rfac = &qfac[n*n];
    REAL* pqw = &rfac[n*np/2];
    INTEGER* iact = (INTEGER*)((char*)ws + roundup(nr*sizeof(REAL),
                                                   sizeof(INTEGER)));

    /* Normalize the constraints, and copy the resultant constraint matrix
     *   and right hand sides into working space, after increasing the right
     *   hand sides if necessary so that the starting point is feasible. */
    iflag = 0;
    if (m > 0) {
      for (j = 1; j <= m; ++j) {
        sum = ZERO;
        temp = ZERO;
        for (i = 1; i <= n; ++i) {
          sum += A(i,j)*X(i);
          temp += pow2(A(i,j));
        }
        if (temp == ZERO) {
          fputs("Return from LINCOA because the gradient "
                "of a constraint is zero.\n", stderr);
          return;
        }
        temp = sqrt(temp);
        if (sum - B(j) > smallx*temp) {
          iflag = 1;
        }
        WB(j) = max(B(j), sum)/temp;
        for (i = 1; i <= n; ++i) {
          AMAT(i,j) = A(i,j)/temp;
        }
      }
    }
    if (iflag == 1) {
      if (iprint > 0) {
        fputs("LINCOA has made the initial X feasible "
              "by increasing part(s) of B.\n", stderr);
      }
    }

    /* The above settings provide a partition of W for subroutine LINCOB. */
    lincob(n, npt, m, objfun, data, amat, wb, x, rhobeg, rhoend,
           iprint, maxfun, xbase, xpt, fval, xsav, xopt, gopt,
           hq, pq, bmat, zmat, ndim, step, sp, xnew, iact, rescon,
           qfac, rfac, pqw, w);
  }
}


#ifdef TESTING

/* Variable used to store the fcuntion value when `x` is infeasible. */
static REAL fmax;

/* Objective function. */
REAL
objfun(const INTEGER n, const REAL x[], LOGICAL feasible, void* data)
{
  REAL f, v12, v13, v14, v23, v24, v34, del1, del2, del3, del4, temp;

  v12 = X(1)*X(5) - X(4)*X(2);
  v13 = X(1)*X(8) - X(7)*X(2);
  v14 = X(1)*X(11) - X(10)*X(2);
  v23 = X(4)*X(8) - X(7)*X(5);
  v24 = X(4)*X(11) - X(10)*X(5);
  v34 = X(7)*X(11) - X(10)*X(8);
  del1 = v23*X(12) - v24*X(9) + v34*X(6);
  if (del1 <= ZERO) {
    return fmax;
  }
  del2 = -v34*X(3) - v13*X(12) + v14*X(9);
  if (del2 <= ZERO) {
    return fmax;
  }
  del3 = -v14*X(6) + v24*X(3) + v12*X(12);
  if (del3 <= ZERO) {
    return fmax;
  }
  del4 = -v12*X(9) + v13*X(6) - v23*X(3);
  if (del4 <= ZERO) {
    return fmax;
  }
  temp = del1 + del2 + del3 + del4;
  temp = temp*(temp*temp)/(del1*del2*del3*del4);
  temp = temp/6.0;
  f = min(temp,fmax);
  return f;
}


/* Calculate the tetrahedron of least volume that encloses the points
 *   (XP(J),YP(J),ZP(J)), J=1,2,...,NP. Our method requires the origin
 *   to be strictly inside the convex hull of these points. There are
 *   twelve variables that define the four faces of each tetrahedron
 *   that is considered. Each face has the form ALPHA*X+BETA*Y+GAMMA*Z=1,
 *   the variables X(3K-2), X(3K-1) and X(3K) being the values of ALPHA,
 *   BETA and GAMMA for the K-th face, K=1,2,3,4. Let the set T contain
 *   all points in three dimensions that can be reached from the origin
 *   without crossing a face. Because the volume of T may be infinite,
 *   the objective function is the smaller of FMAX and the volume of T,
 *   where FMAX is set to an upper bound on the final volume initially.
 *   There are 4*NP linear constraints on the variables, namely that each
 *   of the given points (XP(J),YP(J),ZP(J)) shall be in T. Let XS = min
 *   XP(J), YS = min YP(J), ZS = min ZP(J) and SS = max XP(J)+YP(J)+ZP(J),
 *   where J runs from 1 to NP. The initial values of the variables are
 *   X(1)=1/XS, X(5)=1/YS, X(9)=1/ZS, X(2)=X(3)=X(4)=X(6)=X(7) =X(8)=0
 *   and X(10)=X(11)=X(12)=1/SS, which satisfy the linear constraints,
 *   and which provide the bound FMAX=(SS-XS-YS-ZS)**3/6. Other details
 *   of the test calculation are given below, including the choice of
 *   the data points (XP(J),YP(J),ZP(J)), J=1,2,...,NP. The smaller final
 *   value of the objective function in the case NPT=35 shows that the
 *   problem has local minima.
 */
int main(int argc, char* argv[])
{
  /* Set some constants. */
  const REAL two = 2.0;
  const INTEGER n = 12;
  const INTEGER ia = n;
  const INTEGER np = 50;
  const INTEGER m = 4*np;

  /* Local variables */
  REAL a[ia*m], b[m], x[n];
  REAL rhobeg, rhoend;
  REAL ss, xp[np], yp[np], zp[np], xs, ys, zs;
  REAL theta, temp;
  REAL sumx, sumy, sumz;
  void* ws;
  INTEGER i, j, k;
  INTEGER jcase;
  INTEGER maxfun, iprint;
  INTEGER iw;
  INTEGER npt;
  size_t nbytes;

  /* Set the data points. */
  sumx = ZERO;
  sumy = ZERO;
  sumz = ZERO;
  for (j = 1; j <= np; ++j) {
    theta = (REAL)(j - 1)*PI/(REAL)(np - 1);
    XP(j) = cos(theta)*cos(two*theta);
    sumx += XP(j);
    YP(j) = sin(theta)*cos(two*theta);
    sumy += YP(j);
    ZP(j) = sin(two*theta);
    sumz += ZP(j);
  }
  sumx /= (REAL)np;
  sumy /= (REAL)np;
  sumz /= (REAL)np;
  for (j = 1; j <= np; ++j) {
    XP(j) -= sumx;
    YP(j) -= sumy;
    ZP(j) -= sumz;
  }

  /* Set the linear constraints. */
  for (k = 1; k <= m; ++k) {
    B(k) = ONE;
    for (i = 1; i <= n; ++i) {
      A(i,k) = ZERO;
    }
  }
  for (j = 1; j <= np; ++j) {
    for (i = 1; i <= 4; ++i) {
      k = 4*j + i - 4;
      iw = i*3;
      A(iw-2,k) = XP(j);
      A(iw-1,k) = YP(j);
      A(iw,k) = ZP(j);
    }
  }

  /* Set the initial vector of variables. The JCASE=1,6 loop gives six
   *   different choices of NPT when LINCOA is called. */
  xs = ZERO;
  ys = ZERO;
  zs = ZERO;
  ss = ZERO;
  for (j = 1; j <= np; ++j) {
    xs = min(xs, XP(j));
    ys = min(ys, YP(j));
    zs = min(zs, ZP(j));
    ss = max(ss, XP(j) + YP(j) + ZP(j));
  }
  temp = ss - xs - ys - zs;
  fmax = temp*(temp*temp)/6.0;
  for (jcase = 1; jcase <= 6; ++jcase) {
    for (i = 2; i <= 8; ++i) {
      X(i) = ZERO;
    }
    X(1) = ONE/xs;
    X(5) = ONE/ys;
    X(9) = ONE/zs;
    X(10) = ONE/ss;
    X(11) = ONE/ss;
    X(12) = ONE/ss;

    /* Choose number of points and allocate workspace. */
    npt = jcase*5 + 10;
    nbytes = lincoa_storage(n, npt, m);
    ws = malloc(nbytes);
    if (ws == NULL) {
      fprintf(stderr, "cannot allocate %ld bytes for workspace\n",
              (long)nbytes);
      return 1;
    }

    /* Call of LINCOA, which provides the printing given at the end of this
     *   note. */
    rhobeg = 1.0;
    rhoend = 1e-6;
    iprint = 1;
    maxfun = 10000;
    fprintf(OUTPUT,
          "\n\n    Output from LINCOA with  NPT =%4ld  and  RHOEND =%12.4E\n",
          (long)npt, (double)rhoend);
    lincoa(n, npt, m, objfun, NULL, a, ia, b, x,
           rhobeg, rhoend, iprint, maxfun, ws);

    /* Free workspace. */
    free(ws);
  }

  /* Return. */
  return 0;
}

#endif /* TESTING */
