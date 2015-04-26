/*
 * newuoa.c -
 *
 * C implemtation of Mike Powell's NEWUOA algorithm for minimizing a function
 * of many variables.  The method is "derivatives free" (only the function
 * values are needed).  The algorithm is described in:
 *
 *   M.J.D. Powell, "The NEWUOA software for unconstrained minimization
 *   without derivatives", in Large-Scale Nonlinear Optimization, editors
 *   G. Di Pillo and M. Roma, Springer (2006), pages 255-297.
 *
 * The present code is based on the original FORTRAN version written by Mike
 * Powell who kindly provides his code on demand (at mjdp@cam.ac.uk) and has
 * been converted to C by É. Thiébaut.
 *
 * Copyright (c) 2004, Mike Powell (FORTRAN version).
 * Copyright (c) 2015, Éric Thiébaut (C version).
 *
 * Read the accompanying `LICENSE` file for details.
 */

#include <stdio.h>
#include <math.h>

#include "newuoa.h"

#define OUTPUT       stdout

/* Macros to deal with single/double precision. */
#undef REAL
#ifdef SINGLE_PRECISION
# define FLT(x)    x ## f  /* floating point literal constant */
# define REAL      float
# define ABS(x)    fabsf(x)
# define SQRT(x)   sqrtf(x)
# define HYPOT(x)  hypotf(x)
# define LOG(x)    logf(x)
# define EXP(x)    expf(x)
# define SIN(x)    sinf(x)
# define COS(x)    cosf(x)
# define TAN(x)    tanf(x)
# define ASIN(x)   asinf(x)
# define ACOS(x)   acosf(x)
# define ATAN(x)   atanf(x)
#else
# define FLT(x)    x       /* floating point literal constant */
# define REAL      double
# define ABS(x)    fabs(x)
# define SQRT(x)   sqrt(x)
# define HYPOT(x)  hypot(x)
# define LOG(x)    log(x)
# define EXP(x)    exp(x)
# define SIN(x)    sin(x)
# define COS(x)    cos(x)
# define TAN(x)    tan(x)
# define ASIN(x)   asin(x)
# define ACOS(x)   acos(x)
# define ATAN(x)   atan(x)
#endif

/* Helper macro for simple FORTRAN-like loops. */
#define LOOP(var,num)    for (var = 1; var <= num; ++var)

/* Macros yielding the min/max of two values. */
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define MIN(a,b) ((a) <= (b) ? (a) : (b))

/* Macro `HOW_MANY(a,b)` yields the minimal number of chunks with `b` elements
   needed to store `a` elements.  Both `a` and `b` must be integers. */
#define HOW_MANY(a,b)  ((((b) - 1) + (a))/(b))

/* Macro `ROUND_UP(a,b)` yields the integer `a` rounded up to a multiple of
   integer `b`. */
#define ROUND_UP(a,b)  (HOW_MANY(a,b)*(b))

/* Macro `ADDRESS(type,base,offset)` yields a `type*` address at `offset` (in
   bytes) from `base` address. */
#define ADDRESS(type, base, offset) ((type*)((char*)(base) + (offset)))

/* Macro `OFFSET(type, field)` yields the offset (in bytes) of member `field`
   in structure/class `type`. */
#ifdef offsetof
#  define OFFSET(type, field)  offsetof(type, field)
#else
#  define OFFSET(type, field)  ((char*)&((type*)0)->field - (char*)0)
#endif

/*---------------------------------------------------------------------------*/
/* DECLARATIONS OF PRIVATE FUNCTIONS */

static void
print_error(const char* reason);

static void
print_x(FILE* output, INTEGER n, const REAL x[], const REAL dx[]);


/*---------------------------------------------------------------------------*/
/* TESTING */

static REAL
objfun_test(const INTEGER n, const REAL* x, void* data);

#ifdef TESTING

int
main(int argc, char* argv[])
{
  newuoa_test();
  return 0;
}

int
calfun_(const INTEGER* n, REAL* x, REAL* f)
{
  *f = objfun_test(*n, x, NULL);
  return 0;
}

#endif

void
newuoa_test(void)
{
  REAL rhobeg, rhoend;
  REAL w[10000], x[10];
  INTEGER i, n, npt, maxfun, iprint;

  iprint = 2;
  maxfun = 5000;
  rhoend = 1e-6;
  for (n = 2; n <= 8; n += 2) {
    npt = 2*n + 1;
    LOOP(i,n) {
      x[i - 1] = (REAL)i/(REAL)(n + 1);
    }
    rhobeg = x[0]*0.2;
    fprintf(stdout, "\n\n    Results with N =%2d and NPT =%3d\n",
            (int)n, (int)npt);
#if 1
    /* Test the C-code. */
    newuoa(n, npt, objfun_test, NULL, x, rhobeg, rhoend, iprint, maxfun, w);
#else
   /* Test the FORTRAN wrapper. */
   newuoa_(&n, &npt, x, &rhobeg, &rhoend, &iprint, &maxfun, w);
#endif
  }
} /* newuoa_test */


/* The Chebyquad test problem (Fletcher, 1965) for N = 2,4,6 and 8,
   with NPT = 2N+1. */
static REAL
objfun_test(const INTEGER n, const REAL* x, void* data)
{
  REAL sum, f;
  INTEGER i, iw, j, np;
  REAL y[100];

#define X(a1)    x[a1 - 1]
#define Y(a1,a2) y[(a2)*10 + a1 - 11]
  LOOP(j,n) {
    Y(1,j) = 1.0;
    Y(2,j) = X(j)*2.0 - 1.0;
  }
  for (i = 2; i <= n; ++i) {
    LOOP(j,n) {
      Y(i+1,j) = Y(2,j)*2.0*Y(i,j) - Y(i-1,j);
    }
  }
  f = 0.0;
  np = n + 1;
  iw = 1;
  LOOP(i,np) {
    sum = 0.0;
    LOOP(j,n) {
      sum += Y(i,j);
    }
    sum /= (REAL)n;
    if (iw > 0) {
      sum += 1.0/(REAL)(i*i - 2*i);
    }
    iw = -iw;
    f += sum*sum;
  }
  return f;
} /* objfun_test */

#undef X
#undef Y

/*---------------------------------------------------------------------------*/
/* NEWUOA DRIVER ROUTINES */

static int
newuob(const INTEGER n, const INTEGER npt,
       newuoa_objfun* objfun, void* data,
       REAL* x, const REAL rhobeg, const REAL rhoend,
       const INTEGER iprint, const INTEGER maxfun,
       REAL* xbase, REAL* xopt, REAL* xnew,
       REAL* xpt, REAL* fval, REAL* gq, REAL* hq,
       REAL* pq, REAL* bmat, REAL* zmat, const INTEGER ndim,
       REAL* d, REAL* vlag, REAL* w);

int
newuoa(const INTEGER n, const INTEGER npt,
       newuoa_objfun* objfun, void* data,
       REAL* x, const REAL rhobeg, const REAL rhoend,
       const INTEGER iprint, const INTEGER maxfun,
       REAL* w)
{
  INTEGER id, np, iw, igq, ihq, ixb, ifv, ipq, ivl, ixn, ixo, ixp, ndim,
    nptm, ibmat, izmat;

  /* Partition the working space array, so that different parts of it can be
     treated separately by the subroutine that performs the main
     calculation. */
  np = n + 1;
  nptm = npt - np;
  if (npt < n + 2 || npt > (n + 2)*np/2) {
    if (iprint > 0) {
      print_error("NPT is not in the required interval");
    }
    return NEWUOA_BAD_NPT;
  }
  ndim = npt + n;
  ixb = 0; /* C-indices start at 0 */
  ixo = ixb + n;
  ixn = ixo + n;
  ixp = ixn + n;
  ifv = ixp + n*npt;
  igq = ifv + npt;
  ihq = igq + n;
  ipq = ihq + n*np/2;
  ibmat = ipq + npt;
  izmat = ibmat + ndim*n;
  id = izmat + npt*nptm;
  ivl = id + n;
  iw = ivl + ndim;

  /* The above settings provide a partition of W for subroutine NEWUOB.  The
     partition requires the first NPT*(NPT+N)+5*N*(N+3)/2 elements of W plus
     the space that is needed by the last array of NEWUOB. */
  return newuob(n, npt, objfun, data, x, rhobeg, rhoend, iprint, maxfun,
                &w[ixb], &w[ixo], &w[ixn], &w[ixp], &w[ifv], &w[igq],
                &w[ihq], &w[ipq], &w[ibmat], &w[izmat], ndim, &w[id],
                &w[ivl], &w[iw]);
} /* newuoa */

REAL
newuoa_calfun_wrapper(const INTEGER n, const REAL* x, void* data)
{
  REAL f;
  calfun_(&n, (REAL*)x, &f);
  return f;
}

int
newuoa_(const INTEGER* n, const INTEGER* npt, REAL* x,
        const REAL* rhobeg, const REAL* rhoend,
        const INTEGER* iprint, const INTEGER* maxfun,
        REAL* w)
{
  newuoa(*n, *npt, newuoa_calfun_wrapper, NULL,
         x, *rhobeg, *rhoend, *iprint, *maxfun, w);
  return 0;
}

/*---------------------------------------------------------------------------*/
/* NEWUOA SUBROUTINES */

static void
bigden(const INTEGER n, const INTEGER npt, REAL* xopt,
       REAL* xpt, REAL* bmat, REAL* zmat, const INTEGER idz,
       const INTEGER ndim, const INTEGER kopt, const INTEGER knew, REAL* d,
       REAL* w, REAL* vlag, REAL* beta, REAL* s,
       REAL* wvec, REAL* prod);

static void
biglag(const INTEGER n, const INTEGER npt, REAL* xopt,
       REAL* xpt, REAL* bmat, REAL* zmat, INTEGER* idz,
       const INTEGER ndim, const INTEGER knew, const REAL delta, REAL* d,
       REAL* alpha, REAL* hcol, REAL* gc, REAL* gd,
       REAL* s, REAL* w);

static void
update(const INTEGER n, const INTEGER npt, REAL* bmat,
       REAL* zmat, INTEGER* idz, const INTEGER ndim, REAL* vlag,
       const REAL beta, const INTEGER knew, REAL* w);

static void
trsapp(const INTEGER n, const INTEGER npt, REAL* xopt,
       REAL* xpt, REAL* gq, REAL* hq, REAL* pq,
       const REAL delta, REAL* step, REAL* d, REAL* g,
       REAL* hd, REAL* hs, REAL* crvmin);

static int
newuob(const INTEGER n, const INTEGER npt,
       newuoa_objfun* objfun, void* data,
       REAL* x, const REAL rhobeg, const REAL rhoend,
       const INTEGER iprint, const INTEGER maxfun,
       REAL* xbase, REAL* xopt, REAL* xnew,
       REAL* xpt, REAL* fval, REAL* gq, REAL* hq,
       REAL* pq, REAL* bmat, REAL* zmat, const INTEGER ndim,
       REAL* d, REAL* vlag, REAL* w)
{
  /* The arguments N, NPT, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical
     to the corresponding arguments in SUBROUTINE NEWUOA.

     XBASE will hold a shift of origin that should reduce the contributions
     from rounding errors to values of the model and Lagrange functions.

     XOPT will be set to the displacement from XBASE of the vector of variables
     that provides the least calculated F so far.

     XNEW will be set to the displacement from XBASE of the vector of variables
     for the current calculation of F.

     XPT will contain the interpolation point coordinates relative to XBASE.

     FVAL will hold the values of F at the interpolation points.

     GQ will hold the gradient of the quadratic model at XBASE.

     HQ will hold the explicit second derivatives of the quadratic model.

     PQ will contain the parameters of the implicit second derivatives of the
     quadratic model.

     BMAT will hold the last N columns of H.

     ZMAT will hold the factorization of the leading NPT by NPT submatrix of H,
     this factorization being ZMAT times Diag(DZ) times ZMAT^T, where the
     elements of DZ are plus or minus one, as specified by IDZ.

     NDIM is the first dimension of BMAT and has the value NPT+N.

     D is reserved for trial steps from XOPT.

     VLAG will contain the values of the Lagrange functions at a new point X.
     They are part of a product that requires VLAG to be of length NDIM.

     The array W will be used for working space. Its length must be at least
     10*NDIM = 10*(NPT+N). */

  /* Constants. */
  const REAL one = 1.0;
  const REAL zero = 0.0;
  const REAL half = 0.5;
  const REAL tenth = 0.1;

  /* Local variables. */
  REAL alpha, beta, bsum, crvmin, delta, detrat, diff, diffa, diffb, diffc,
    distsq, dnorm, dsq, dstep, dx, f, fbeg, fopt, fsave, gisq, gqsq, hdiag,
    ratio, recip, reciq, rho, rhosq, sum, suma, sumb, sumz, temp, tempa,
    tempb, tempq,
    vquad, xipt, xjpt, xoptsq;
  INTEGER i, idz, ih, ip, ipt, itemp, itest, j, jp, jpt, k, knew, kopt,
    ksave, ktemp, nf, nfm, nfmm, nfsav, nftest, nh, np, nptm;
  int status;
  const char* reason;

  /* Parameter adjustments to comply with FORTRAN indexing. */
  x     -= 1;
  xbase -= 1;
  xopt  -= 1;
  xnew  -= 1;
  xpt   -= 1 + npt;
  fval  -= 1;
  gq    -= 1;
  hq    -= 1;
  pq    -= 1;
  bmat  -= 1 + ndim;
  zmat  -= 1 + npt;
  d     -= 1;
  vlag  -= 1;
  w     -= 1;
#define XPT(a1,a2) xpt[(a2)*npt + a1]
#define BMAT(a1,a2) bmat[(a2)*ndim + a1]
#define ZMAT(a1,a2) zmat[(a2)*npt + a1]

  /* FIXME: Set uninitialized variables. */
  idz = 0;
  ipt = 0;
  itest = 0;
  jpt = 0;
  knew = 0;
  kopt = 0;
  nfsav = 0;
  alpha = zero;
  crvmin = zero;
  delta = zero;
  diffa = zero;
  diffb = zero;
  diffc = zero;
  dnorm = zero;
  f = zero;
  fbeg = zero;
  fopt = zero;
  ratio = zero;
  rho = zero;
  xoptsq = zero;

  /* Initialization. */
  np = n + 1;
  nh = n*np/2;
  nptm = npt - np;
  nftest = MAX(maxfun,1);
  status = NEWUOA_SUCCESS;
  reason = NULL;

  /* Set the initial elements of XPT, BMAT, HQ, PQ and ZMAT to zero. */
  LOOP(j,n) {
    xbase[j] = x[j];
    LOOP(k,npt) {
      XPT(k,j) = zero;
    }
    LOOP(i,ndim) {
      BMAT(i,j) = zero;
    }
  }
  LOOP(ih,nh) {
    hq[ih] = zero;
  }
  LOOP(k,npt) {
    pq[k] = zero;
    LOOP(j,nptm) {
      ZMAT(k,j) = zero;
    }
  }

  /* Begin the initialization procedure. NF becomes one more than the number of
     function values so far. The coordinates of the displacement of the next
     initial interpolation point from XBASE are set in XPT(NF,.). */
  rhosq = rhobeg*rhobeg;
  recip = one/rhosq;
  reciq = SQRT(half)/rhosq;
  nf = 0;
 L50:
  nfm = nf;
  nfmm = nf - n;
  ++nf;
  if (nfm <= 2*n) {
    if (nfm >= 1 && nfm <= n) {
      XPT(nf,nfm) = rhobeg;
    } else if (nfm > n) {
      XPT(nf,nfmm) = -rhobeg;
    }
  } else {
    itemp = (nfmm - 1)/n;
    jpt = nfm - itemp*n - n;
    ipt = jpt + itemp;
    if (ipt > n) {
      itemp = jpt;
      jpt = ipt - n;
      ipt = itemp;
    }
    xipt = rhobeg;
    if (fval[ipt + np] < fval[ipt + 1]) {
      xipt = -xipt;
    }
    xjpt = rhobeg;
    if (fval[jpt + np] < fval[jpt + 1]) {
      xjpt = -xjpt;
    }
    XPT(nf,ipt) = xipt;
    XPT(nf,jpt) = xjpt;
  }

  /* Calculate the next value of F, label 70 being reached immediately after
     this calculation. The least function value so far and its index are
     required. */
  LOOP(j,n) {
    x[j] = XPT(nf,j) + xbase[j];
  }
  goto L310;
 L70:
  fval[nf] = f;
  if (nf == 1) {
    fbeg = f;
    fopt = f;
    kopt = 1;
  } else if (f < fopt) {
    fopt = f;
    kopt = nf;
  }

  /* Set the nonzero initial elements of BMAT and the quadratic model in the
     cases when NF is at most 2*N+1. */
  if (nfm <= 2*n) {
    if (nfm >= 1 && nfm <= n) {
      gq[nfm] = (f - fbeg)/rhobeg;
      if (npt < nf + n) {
        BMAT(1,nfm) = -one/rhobeg;
        BMAT(nf,nfm) = one/rhobeg;
        BMAT(npt + nfm, nfm) = -half*rhosq;
      }
    } else if (nfm > n) {
      BMAT(nf - n, nfmm) = half/rhobeg;
      BMAT(nf,nfmm) = -half/rhobeg;
      ZMAT(1,nfmm) = -reciq - reciq;
      ZMAT(nf - n, nfmm) = reciq;
      ZMAT(nf,nfmm) = reciq;
      ih = nfmm*(nfmm + 1)/2;
      temp = (fbeg - f)/rhobeg;
      hq[ih] = (gq[nfmm] - temp)/rhobeg;
      gq[nfmm] = half*(gq[nfmm] + temp);
    }
  } else {
    /* Set the off-diagonal second derivatives of the Lagrange functions and
       the initial quadratic model. */
    ih = ipt*(ipt - 1)/2 + jpt;
    if (xipt < zero) {
      ipt += n;
    }
    if (xjpt < zero) {
      jpt += n;
    }
    ZMAT(1,nfmm) = recip;
    ZMAT(nf,nfmm) = recip;
    ZMAT(ipt + 1, nfmm) = -recip;
    ZMAT(jpt + 1, nfmm) = -recip;
    hq[ih] = (fbeg - fval[ipt + 1] - fval[jpt + 1] + f)/(xipt*xjpt);
  }
  if (nf < npt) {
    goto L50;
  }

  /* Begin the iterative procedure, because the initial model is complete. */
  rho = rhobeg;
  delta = rho;
  idz = 1;
  diffa = zero;
  diffb = zero;
  itest = 0;
  xoptsq = zero;
  LOOP(i,n) {
    xopt[i] = XPT(kopt,i);
    xoptsq += xopt[i]*xopt[i];
  }
 L90:
  nfsav = nf;

  /* Generate the next trust region step and test its length. Set KNEW to -1 if
     the purpose of the next F will be to improve the model. */
 L100:
  knew = 0;
  trsapp(n, npt, &xopt[1], &XPT(1,1), &gq[1], &hq[1], &pq[1],
         delta, &d[1], &w[1], &w[np], &w[np + n],
         &w[np + 2*n], &crvmin);
  dsq = zero;
  LOOP(i,n) {
    dsq += d[i]*d[i];
  }
  dnorm = SQRT(dsq);
  dnorm = MIN(dnorm,delta);
  if (dnorm < half*rho) {
    knew = -1;
    delta = tenth*delta;
    ratio = -1.0;
    if (delta <= rho*1.5) {
      delta = rho;
    }
    if (nf <= nfsav + 2) {
      goto L460;
    }
    temp = crvmin*0.125*rho*rho;
    tempa = MAX(diffa,diffb);
    if (temp <= MAX(tempa,diffc)) {
      goto L460;
    }
    goto L490;
  }

  /* Shift XBASE if XOPT may be too far from XBASE. First make the changes
     to BMAT that do not depend on ZMAT. */
 L120:
  if (dsq <= xoptsq*0.001) {
    tempq = xoptsq*0.25;
    LOOP(k,npt) {
      sum = zero;
      LOOP(i,n) {
        sum += XPT(k,i)*xopt[i];
      }
      temp = pq[k]*sum;
      sum -= half*xoptsq;
      w[npt + k] = sum;
      LOOP(i,n) {
        gq[i] += temp*XPT(k,i);
        XPT(k,i) = XPT(k,i) - half*xopt[i];
        vlag[i] = BMAT(k,i);
        w[i] = sum*XPT(k,i) + tempq*xopt[i];
        ip = npt + i;
        LOOP(j,i) {
          BMAT(ip,j) = BMAT(ip,j) + vlag[i]*w[j] + w[i]*vlag[j];
        }
      }
    }

    /* Then the revisions of BMAT that depend on ZMAT are calculated. */
    LOOP(k,nptm) {
      sumz = zero;
      LOOP(i,npt) {
        sumz += ZMAT(i,k);
        w[i] = w[npt + i]*ZMAT(i,k);
      }
      LOOP(j,n) {
        sum = tempq*sumz*xopt[j];
        LOOP(i,npt) {
          sum += w[i]*XPT(i,j);
        }
        vlag[j] = sum;
        if (k < idz) {
          sum = -sum;
        }
        LOOP(i,npt) {
          BMAT(i,j) = BMAT(i,j) + sum*ZMAT(i,k);
        }
      }
      LOOP(i,n) {
        ip = i + npt;
        temp = vlag[i];
        if (k < idz) {
          temp = -temp;
        }
        LOOP(j,i) {
          BMAT(ip,j) = BMAT(ip,j) + temp*vlag[j];
        }
      }
    }

    /* The following instructions complete the shift of XBASE, including the
       changes to the parameters of the quadratic model. */
    ih = 0;
    LOOP(j,n) {
      w[j] = zero;
      LOOP(k,npt) {
        w[j] += pq[k]*XPT(k,j);
        XPT(k,j) = XPT(k,j) - half*xopt[j];
      }
      LOOP(i,j) {
        ++ih;
        if (i < j) {
          gq[j] += hq[ih]*xopt[i];
        }
        gq[i] += hq[ih]*xopt[j];
        hq[ih] = hq[ih] + w[i]*xopt[j] + xopt[i]*w[j];
        BMAT(npt + i, j) = BMAT(npt + j, i);
      }
    }
    LOOP(j,n) {
      xbase[j] += xopt[j];
      xopt[j] = zero;
    }
    xoptsq = zero;
  }

  /* Pick the model step if KNEW is positive. A different choice of D may be
     made later, if the choice of D by BIGLAG causes substantial cancellation
     in DENOM. */
  if (knew > 0) {
    biglag(n, npt, &xopt[1], &XPT(1,1), &BMAT(1,1),
           &ZMAT(1,1), &idz, ndim, knew, dstep, &d[1],
           &alpha, &vlag[1], &vlag[npt + 1], &w[1], &w[np], &w[np + n]);
  }

  /* Calculate VLAG and BETA for the current choice of D. The first NPT
     components of W_check will be held in W. */
  LOOP(k,npt) {
    suma = zero;
    sumb = zero;
    sum = zero;
    LOOP(j,n) {
      suma += XPT(k,j)*d[j];
      sumb += XPT(k,j)*xopt[j];
      sum += BMAT(k,j)*d[j];
    }
    w[k] = suma*(half*suma + sumb);
    vlag[k] = sum;
  }
  beta = zero;
  LOOP(k,nptm) {
    sum = zero;
    LOOP(i,npt) {
      sum += ZMAT(i,k)*w[i];
    }
    if (k < idz) {
      beta += sum*sum;
      sum = -sum;
    } else {
      beta -= sum*sum;
    }
    LOOP(i,npt) {
      vlag[i] += sum*ZMAT(i,k);
    }
  }
  bsum = zero;
  dx = zero;
  LOOP(j,n) {
    sum = zero;
    LOOP(i,npt) {
      sum += w[i]*BMAT(i,j);
    }
    bsum += sum*d[j];
    jp = npt + j;
    LOOP(k,n) {
      sum += BMAT(jp,k)*d[k];
    }
    vlag[jp] = sum;
    bsum += sum*d[j];
    dx += d[j]*xopt[j];
  }
  beta = dx*dx + dsq*(xoptsq + dx + dx + half*dsq) + beta - bsum;
  vlag[kopt] += one;

  /* If KNEW is positive and if the cancellation in DENOM is unacceptable, then
     BIGDEN calculates an alternative model step, XNEW being used for working
     space. */
  if (knew > 0) {
    temp = one + alpha*beta/(vlag[knew]*vlag[knew]);
    if (ABS(temp) <= 0.8) {
      bigden(n, npt, &xopt[1], &XPT(1,1), &BMAT(1,1),
             &ZMAT(1,1), idz, ndim, kopt, knew, &d[1], &w[1],
             &vlag[1], &beta, &xnew[1], &w[ndim + 1], &w[ndim*6 + 1]);
    }
  }

  /* Calculate the next value of the objective function. */
 L290:
  LOOP(i,n) {
    xnew[i] = xopt[i] + d[i];
    x[i] = xbase[i] + xnew[i];
  }
  ++nf;
 L310:
  if (nf > nftest) {
    --nf;
    reason = "CALFUN has been called MAXFUN times";
    status = NEWUOA_TOO_MANY_EVALUATIONS;
    goto done;
  }
  f = objfun(n, &x[1], data);
  if (iprint == 3) {
    fprintf(OUTPUT, "\n"
            "    Function number%6ld    F =%18.10E"
            "    The corresponding X is:\n",
            (long)nf, (double)f);
    print_x(OUTPUT, n, &x[1], NULL);
  }
  if (nf <= npt) {
    goto L70;
  }
  if (knew == -1) {
    goto done;
  }

  /* Use the quadratic model to predict the change in F due to the step D, and
     set DIFF to the error of this prediction. */
  vquad = zero;
  ih = 0;
  LOOP(j,n) {
    vquad += d[j]*gq[j];
    LOOP(i,j) {
      ++ih;
      temp = d[i]*xnew[j] + d[j]*xopt[i];
      if (i == j) {
        temp = half*temp;
      }
      vquad += temp*hq[ih];
    }
  }
  LOOP(k,npt) {
    vquad += pq[k]*w[k];
  }
  diff = f - fopt - vquad;
  diffc = diffb;
  diffb = diffa;
  diffa = ABS(diff);
  if (dnorm > rho) {
    nfsav = nf;
  }

  /* Update FOPT and XOPT if the new F is the least value of the objective
     function so far. The branch when KNEW is positive occurs if D is not a
     trust region step. */
  fsave = fopt;
  if (f < fopt) {
    fopt = f;
    xoptsq = zero;
    LOOP(i,n) {
      xopt[i] = xnew[i];
      xoptsq += xopt[i]*xopt[i];
    }
  }
  ksave = knew;
  if (knew > 0) {
    goto L410;
  }

  /* Pick the next value of DELTA after a trust region step. */
  if (vquad >= zero) {
    reason = "a trust region step has failed to reduce Q";
    status = NEWUOA_STEP_FAILED;
    goto done;
  }
  ratio = (f - fsave)/vquad;
  if (ratio <= tenth) {
    delta = half*dnorm;
  } else if (ratio <= 0.7) {
    delta *= half;
    delta = MAX(delta,dnorm);
  } else {
    tempa = dnorm + dnorm;
    delta *= half;
    delta = MAX(delta,tempa);
  }
  if (delta <= rho*1.5) {
    delta = rho;
  }

  /* Set KNEW to the index of the next interpolation point to be deleted. */
  tempa = tenth*delta;
  tempa = MAX(tempa,rho);
  rhosq = tempa*tempa;
  ktemp = 0;
  detrat = zero;
  if (f >= fsave) {
    ktemp = kopt;
    detrat = one;
  }
  LOOP(k,npt) {
    hdiag = zero;
    LOOP(j,nptm) {
      temp = one;
      if (j < idz) {
        temp = -one;
      }
      hdiag += temp*(ZMAT(k,j)*ZMAT(k,j));
    }
    temp = ABS(beta*hdiag + vlag[k]*vlag[k]);
    distsq = zero;
    LOOP(j,n) {
      tempa = XPT(k,j) - xopt[j];
      distsq += tempa*tempa;
    }
    if (distsq > rhosq) {
      tempa = distsq/rhosq;
      temp *= (tempa*tempa)*tempa;
    }
    if (temp > detrat && k != ktemp) {
      detrat = temp;
      knew = k;
    }
  }
  if (knew == 0) {
    goto L460;
  }

  /* Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point can be
     moved. Begin the updating of the quadratic model, starting with the
     explicit second derivative term. */
 L410:
  update(n, npt, &BMAT(1,1), &ZMAT(1,1), &idz, ndim, &vlag[1],
         beta, knew, &w[1]);
  fval[knew] = f;
  ih = 0;
  LOOP(i,n) {
    temp = pq[knew]*XPT(knew,i);
    LOOP(j,i) {
      ++ih;
      hq[ih] += temp*XPT(knew,j);
    }
  }
  pq[knew] = zero;

  /* Update the other second derivative parameters, and then the gradient
     vector of the model. Also include the new interpolation point. */
  LOOP(j,nptm) {
    temp = diff*ZMAT(knew,j);
    if (j < idz) {
      temp = -temp;
    }
    LOOP(k,npt) {
      pq[k] += temp*ZMAT(k,j);
    }
  }
  gqsq = zero;
  LOOP(i,n) {
    gq[i] += diff*BMAT(knew,i);
    gqsq += gq[i]*gq[i];
    XPT(knew,i) = xnew[i];
  }

  /* If a trust region step makes a small change to the objective function,
     then calculate the gradient of the least Frobenius norm interpolant at
     XBASE, and store it in W, using VLAG for a vector of right hand sides. */
  if (ksave == 0 && delta == rho) {
    if (ABS(ratio) > 0.01) {
      itest = 0;
    } else {
      LOOP(k,npt) {
        vlag[k] = fval[k] - fval[kopt];
      }
      gisq = zero;
      LOOP(i,n) {
        sum = zero;
        LOOP(k,npt) {
          sum += BMAT(k,i)*vlag[k];
        }
        gisq += sum*sum;
        w[i] = sum;
      }

      /* Test whether to replace the new quadratic model by the least Frobenius
         norm interpolant, making the replacement if the test is satisfied. */
      ++itest;
      if (gqsq < gisq*100.0) {
        itest = 0;
      }
      if (itest >= 3) {
        LOOP(i,n) {
          gq[i] = w[i];
        }
        LOOP(ih,nh) {
          hq[ih] = zero;
        }
        LOOP(j,nptm) {
          w[j] = zero;
          LOOP(k,npt) {
            w[j] += vlag[k]*ZMAT(k,j);
          }
          if (j < idz) {
            w[j] = -w[j];
          }
        }
        LOOP(k,npt) {
          pq[k] = zero;
          LOOP(j,nptm) {
            pq[k] += ZMAT(k,j)*w[j];
          }
        }
        itest = 0;
      }
    }
  }
  if (f < fsave) {
    kopt = knew;
  }

  /* If a trust region step has provided a sufficient decrease in F, then
     branch for another trust region calculation. The case KSAVE>0 occurs when
     the new function value was calculated by a model step. */
  if (f <= fsave + tenth*vquad || ksave > 0) {
    goto L100;
  }

  /* Alternatively, find out if the interpolation points are close enough to
     the best point so far. */
  knew = 0;
 L460:
  distsq = delta*4.0*delta;
  LOOP(k,npt) {
    sum = zero;
    LOOP(j,n) {
      tempa = XPT(k,j) - xopt[j];
      sum += tempa*tempa;
    }
    if (sum > distsq) {
      knew = k;
      distsq = sum;
    }
  }

  /* If KNEW is positive, then set DSTEP, and branch back for the next
     iteration, which will generate a "model step". */
  if (knew > 0) {
    tempa = tenth*SQRT(distsq);
    tempb = half*delta;
    tempa = MIN(tempa,tempb);
    dstep = MAX(tempa,rho);
    dsq = dstep*dstep;
    goto L120;
  }
  if (ratio > zero || MAX(delta,dnorm) > rho) {
    goto L100;
  }

  /* The calculations with the current value of RHO are complete. Pick the
     next values of RHO and DELTA. */
 L490:
  if (rho > rhoend) {
    delta = half*rho;
    ratio = rho/rhoend;
    if (ratio <= 16.0) {
      rho = rhoend;
    } else if (ratio <= 250.0) {
      rho = SQRT(ratio)*rhoend;
    } else {
      rho = tenth*rho;
    }
    delta = MAX(delta,rho);
    if (iprint >= 2) {
      if (iprint >= 3) {
        fprintf(OUTPUT, "\n");
      }
      fprintf(OUTPUT, "\n"
              "    New RHO =%11.4E "
              "    Number of function values =%6ld\n"
              "    Least value of F =%23.15E     "
              "    The corresponding X is:\n",
              (double)rho, (long)nf, (double)fopt);
      print_x(OUTPUT, n, &xbase[1], &xopt[1]);
    }
    goto L90;
  }

  /* Return from the calculation, after another Newton-Raphson step, if
     it is too short to have been tried before. */
  if (knew == -1) {
    goto L290;
  }
 done:
  if (fopt <= f) {
    LOOP(i,n) {
      x[i] = xbase[i] + xopt[i];
    }
    f = fopt;
  }
  if (iprint > 0) {
    if (status == NEWUOA_SUCCESS) {
      fprintf(OUTPUT, "\n"
              "    At the return from NEWUOA "
              "    Number of function values =%6ld\n"
              "    Least value of F =%23.15E     "
              "    The corresponding X is:\n",
              (long)nf, (double)f);
      print_x(OUTPUT, n, &x[1], NULL);
    } else {
      print_error(reason);
    }
  }
  return status;
} /* newuob */

#undef ZMAT
#undef BMAT
#undef XPT

static void
bigden(const INTEGER n, const INTEGER npt, REAL* xopt,
       REAL* xpt, REAL* bmat, REAL* zmat, const INTEGER idz,
       const INTEGER ndim, const INTEGER kopt, const INTEGER knew, REAL* d,
       REAL* w, REAL* vlag, REAL* beta, REAL* s,
       REAL* wvec, REAL* prod)
{
  /* N is the number of variables.

     NPT is the number of interpolation equations.

     XOPT is the best interpolation point so far.

     XPT contains the coordinates of the current interpolation points.

     BMAT provides the last N columns of H.

     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.

     NDIM is the first dimension of BMAT and has the value NPT+N.

     KOPT is the index of the optimal interpolation point.

     KNEW is the index of the interpolation point that is going to be moved.

     D will be set to the step from XOPT to the new point, and on entry it
     should be the D that was calculated by the last call of BIGLAG. The length
     of the initial D provides a trust region bound on the final D.

     W will be set to Wcheck for the final choice of D.

     VLAG will be set to Theta*Wcheck+e_b for the final choice of D.

     BETA will be set to the value that will occur in the updating formula when
     the KNEW-th interpolation point is moved to its new position.

     S, WVEC, PROD and the private arrays DEN, DENEX and PAR will be used for
     working space.

     D is calculated in a way that should provide a denominator with a large
     modulus in the updating formula when the KNEW-th interpolation point is
     shifted to the new position XOPT+D. */

  /* Constants. */
  const REAL zero = 0.0;
  const REAL one = 1.0;
  const REAL two = 2.0;
  const REAL half = 0.5;
  const REAL quart = 0.25;
  const REAL twopi = 2.0*M_PI;

  /* Local variables. */
  REAL alpha, angle, dd, denmax, denold, densav, diff, ds, dstemp, dtest,
    ss, ssden, sstemp, step, sum, sumold, tau, temp, tempa, tempb, tempc,
    xoptd, xopts, xoptsq;
  REAL den[9], denex[9], par[9];
  INTEGER i, ip, isave, iterc, iu, j, jc, k, ksav, nptm, nw;

  /* Parameter adjustments to comply with FORTRAN indexing. */
  xopt -= 1;
  xpt  -= 1 + npt;
  bmat -= 1 + ndim;
  zmat -= 1 + npt;
  d    -= 1;
  w    -= 1;
  vlag -= 1;
  s    -= 1;
  wvec -= 1 + ndim;
  prod -= 1 + ndim;
#define XPT(a1,a2)   xpt[(a2)*npt  + a1]
#define BMAT(a1,a2) bmat[(a2)*ndim + a1]
#define WVEC(a1,a2) wvec[(a2)*ndim + a1]
#define PROD(a1,a2) prod[(a2)*ndim + a1]
#define ZMAT(a1,a2) zmat[(a2)*npt  + a1]

  /* Initialization. */
  nptm = npt - n - 1;

  /* Store the first NPT elements of the KNEW-th column of H in W(N+1)
     to W(N+NPT). */
  LOOP(k,npt) {
    w[n + k] = zero;
  }
  LOOP(j,nptm) {
    temp = ZMAT(knew, j);
    if (j < idz) {
      temp = -temp;
    }
    LOOP(k,npt) {
      w[n + k] += temp*ZMAT(k,j);
    }
  }
  alpha = w[n + knew];

  /* The initial search direction D is taken from the last call of BIGLAG,
     and the initial S is set below, usually to the direction from X_OPT
     to X_KNEW, but a different direction to an interpolation point may
     be chosen, in order to prevent S from being nearly parallel to D. */
  dd = zero;
  ds = zero;
  ss = zero;
  xoptsq = zero;
  LOOP(i,n) {
    dd += d[i]*d[i];
    s[i] = XPT(knew, i) - xopt[i];
    ds += d[i]*s[i];
    ss += s[i]*s[i];
    xoptsq += xopt[i]*xopt[i];
  }
  if (ds*ds > dd*0.99*ss) {
    ksav = knew;
    dtest = ds*ds/ss;
    LOOP(k,npt) {
      if (k != kopt) {
        dstemp = zero;
        sstemp = zero;
        LOOP(i,n) {
          diff = XPT(k,i) - xopt[i];
          dstemp += d[i]*diff;
          sstemp += diff*diff;
        }
        if (dstemp*dstemp/sstemp < dtest) {
          ksav = k;
          dtest = dstemp*dstemp/sstemp;
          ds = dstemp;
          ss = sstemp;
        }
      }
    }
    LOOP(i,n) {
      s[i] = XPT(ksav,i) - xopt[i];
    }
  }
  ssden = dd*ss - ds*ds;
  iterc = 0;
  densav = zero;

  /* Begin the iteration by overwriting S with a vector that has the
     required length and direction. */
 L70:
  ++iterc;
  temp = one/SQRT(ssden);
  xoptd = zero;
  xopts = zero;
  LOOP(i,n) {
    s[i] = temp*(dd*s[i] - ds*d[i]);
    xoptd += xopt[i]*d[i];
    xopts += xopt[i]*s[i];
  }

  /* Set the coefficients of the first two terms of BETA. */
  tempa = half*xoptd*xoptd;
  tempb = half*xopts*xopts;
  den[0] = dd*(xoptsq + half*dd) + tempa + tempb;
  den[1] = two*xoptd*dd;
  den[2] = two*xopts*dd;
  den[3] = tempa - tempb;
  den[4] = xoptd*xopts;
  for (i = 6; i <= 9; ++i) {
    den[i - 1] = zero;
  }

  /* Put the coefficients of Wcheck in WVEC. */
  LOOP(k,npt) {
    tempa = zero;
    tempb = zero;
    tempc = zero;
    LOOP(i,n) {
      tempa += XPT(k,i)*d[i];
      tempb += XPT(k,i)*s[i];
      tempc += XPT(k,i)*xopt[i];
    }
    WVEC(k,1) = quart*(tempa*tempa + tempb*tempb);
    WVEC(k,2) = tempa*tempc;
    WVEC(k,3) = tempb*tempc;
    WVEC(k,4) = quart*(tempa*tempa - tempb*tempb);
    WVEC(k,5) = half*tempa*tempb;
  }
  LOOP(i,n) {
    ip = i + npt;
    WVEC(ip,1) = zero;
    WVEC(ip,2) = d[i];
    WVEC(ip,3) = s[i];
    WVEC(ip,4) = zero;
    WVEC(ip,5) = zero;
  }

  /* Put the coefficents of THETA*Wcheck in PROD. */
  for (jc = 1; jc <= 5; ++jc) {
    nw = npt;
    if (jc == 2 || jc == 3) {
      nw = ndim;
    }
    LOOP(k,npt) {
      PROD(k,jc) = zero;
    }
    LOOP(j,nptm) {
      sum = zero;
      LOOP(k,npt) {
        sum += ZMAT(k,j)*WVEC(k,jc);
      }
      if (j < idz) {
        sum = -sum;
      }
      LOOP(k,npt) {
        PROD(k,jc) = PROD(k,jc) + sum*ZMAT(k,j);
      }
    }
    if (nw == ndim) {
      LOOP(k,npt) {
        sum = zero;
        LOOP(j,n) {
          sum += BMAT(k,j)*WVEC(npt + j, jc);
        }
        PROD(k,jc) = PROD(k,jc) + sum;
      }
    }
    LOOP(j,n) {
      sum = zero;
      LOOP(i,nw) {
        sum += BMAT(i,j)*WVEC(i,jc);
      }
      PROD(npt + j, jc) = sum;
    }
  }

  /* Include in DEN the part of BETA that depends on THETA. */
  LOOP(k,ndim) {
    sum = zero;
    for (i = 1; i <= 5; ++i) {
      par[i - 1] = half*PROD(k,i)*WVEC(k,i);
      sum += par[i - 1];
    }
    den[0] = den[0] - par[0] - sum;
    tempa = PROD(k,1)*WVEC(k,2) + PROD(k,2)*WVEC(k,1);
    tempb = PROD(k,2)*WVEC(k,4) + PROD(k,4)*WVEC(k,2);
    tempc = PROD(k,3)*WVEC(k,5) + PROD(k,5)*WVEC(k,3);
    den[1] = den[1] - tempa - half*(tempb + tempc);
    den[5] -= half*(tempb - tempc);
    tempa = PROD(k,1)*WVEC(k,3) + PROD(k,3)*WVEC(k,1);
    tempb = PROD(k,2)*WVEC(k,5) + PROD(k,5)*WVEC(k,2);
    tempc = PROD(k,3)*WVEC(k,4) + PROD(k,4)*WVEC(k,3);
    den[2] = den[2] - tempa - half*(tempb - tempc);
    den[6] -= half*(tempb + tempc);
    tempa = PROD(k,1)*WVEC(k,4) + PROD(k,4)*WVEC(k,1);
    den[3] = den[3] - tempa - par[1] + par[2];
    tempa = PROD(k,1)*WVEC(k,5) + PROD(k,5)*WVEC(k,1);
    tempb = PROD(k,2)*WVEC(k,3) + PROD(k,3)*WVEC(k,2);
    den[4] = den[4] - tempa - half*tempb;
    den[7] = den[7] - par[3] + par[4];
    tempa = PROD(k,4)*WVEC(k,5) + PROD(k,5)*WVEC(k,4);
    den[8] -= half*tempa;
  }

  /* Extend DEN so that it holds all the coefficients of DENOM. */
  sum = zero;
  for (i = 1; i <= 5; ++i) {
    par[i - 1] = half*(PROD(knew,i)*PROD(knew,i));
    sum += par[i - 1];
  }
  denex[0] = alpha*den[0] + par[0] + sum;
  tempa = two*PROD(knew,1)*PROD(knew,2);
  tempb = PROD(knew,2)*PROD(knew,4);
  tempc = PROD(knew,3)*PROD(knew,5);
  denex[1] = alpha*den[1] + tempa + tempb + tempc;
  denex[5] = alpha*den[5] + tempb - tempc;
  tempa = two*PROD(knew,1)*PROD(knew,3);
  tempb = PROD(knew,2)*PROD(knew,5);
  tempc = PROD(knew,3)*PROD(knew,4);
  denex[2] = alpha*den[2] + tempa + tempb - tempc;
  denex[6] = alpha*den[6] + tempb + tempc;
  tempa = two*PROD(knew,1)*PROD(knew,4);
  denex[3] = alpha*den[3] + tempa + par[1] - par[2];
  tempa = two*PROD(knew,1)*PROD(knew,5);
  denex[4] = alpha*den[4] + tempa + PROD(knew,2)*PROD(knew,3);
  denex[7] = alpha*den[7] + par[3] - par[4];
  denex[8] = alpha*den[8] + PROD(knew,4)*PROD(knew,5);

  /* Seek the value of the angle that maximizes the modulus of DENOM. */
  sum = denex[0] + denex[1] + denex[3] + denex[5] + denex[7];
  denold = sum;
  denmax = sum;
  isave = 0;
  iu = 49;
  temp = twopi/(REAL)(iu + 1);
  par[0] = one;
  LOOP(i,iu) {
    angle = (REAL)i*temp;
    par[1] = COS(angle);
    par[2] = SIN(angle);
    for (j = 4; j <= 8; j += 2) {
      par[j - 1] = par[1]*par[j - 3] - par[2]*par[j - 2];
      par[j] = par[1]*par[j - 2] + par[2]*par[j - 3];
    }
    sumold = sum;
    sum = zero;
    for (j = 1; j <= 9; ++j) {
      sum += denex[j - 1]*par[j - 1];
    }
    if (ABS(sum) > ABS(denmax)) {
      denmax = sum;
      isave = i;
      tempa = sumold;
    } else if (i == isave + 1) {
      tempb = sum;
    }
  }
  if (isave == 0) {
    tempa = sum;
  }
  if (isave == iu) {
    tempb = denold;
  }
  step = zero;
  if (tempa != tempb) {
    tempa -= denmax;
    tempb -= denmax;
    step = half*(tempa - tempb)/(tempa + tempb);
  }
  angle = temp*((REAL)isave + step);

  /* Calculate the new parameters of the denominator, the new VLAG vector
     and the new D. Then test for convergence. */
  par[1] = COS(angle);
  par[2] = SIN(angle);
  for (j = 4; j <= 8; j += 2) {
    par[j - 1] = par[1]*par[j - 3] - par[2]*par[j - 2];
    par[j] = par[1]*par[j - 2] + par[2]*par[j - 3];
  }
  *beta = zero;
  denmax = zero;
  for (j = 1; j <= 9; ++j) {
    *beta += den[j - 1]*par[j - 1];
    denmax += denex[j - 1]*par[j - 1];
  }
  LOOP(k,ndim) {
    vlag[k] = zero;
    for (j = 1; j <= 5; ++j) {
      vlag[k] += PROD(k,j)*par[j - 1];
    }
  }
  tau = vlag[knew];
  dd = zero;
  tempa = zero;
  tempb = zero;
  LOOP(i,n) {
    d[i] = par[1]*d[i] + par[2]*s[i];
    w[i] = xopt[i] + d[i];
    dd += d[i]*d[i];
    tempa += d[i]*w[i];
    tempb += w[i]*w[i];
  }
  if (iterc >= n) {
    goto L340;
  }
  if (iterc > 1) {
    densav = MAX(densav,denold);
  }
  if (ABS(denmax) <= ABS(densav)*1.1) {
    goto L340;
  }
  densav = denmax;

  /* Set S to half the gradient of the denominator with respect to D.
     Then branch for the next iteration. */
  LOOP(i,n) {
    temp = tempa*xopt[i] + tempb*d[i] - vlag[npt + i];
    s[i] = tau*BMAT(knew, i) + alpha*temp;
  }
  LOOP(k,npt) {
    sum = zero;
    LOOP(j,n) {
      sum += XPT(k,j)*w[j];
    }
    temp = (tau*w[n + k] - alpha*vlag[k])*sum;
    LOOP(i,n) {
      s[i] += temp*XPT(k,i);
    }
  }
  ss = zero;
  ds = zero;
  LOOP(i,n) {
    ss += s[i]*s[i];
    ds += d[i]*s[i];
  }
  ssden = dd*ss - ds*ds;
  if (ssden >= dd*1e-8*ss) {
    goto L70;
  }

  /* Set the vector W before the RETURN from the subroutine. */
 L340:
  LOOP(k,ndim) {
    w[k] = zero;
    for (j = 1; j <= 5; ++j) {
      w[k] += WVEC(k,j)*par[j - 1];
    }
  }
  vlag[kopt] += one;
  return;
} /* bigden */

#undef ZMAT
#undef PROD
#undef WVEC
#undef BMAT
#undef XPT

static void
biglag(const INTEGER n, const INTEGER npt, REAL* xopt,
       REAL* xpt, REAL* bmat, REAL* zmat, INTEGER* idz,
       const INTEGER ndim, const INTEGER knew, const REAL delta, REAL* d,
       REAL* alpha, REAL* hcol, REAL* gc, REAL* gd,
       REAL* s, REAL* w)
{
  /* N is the number of variables.

     NPT is the number of interpolation equations.

     XOPT is the best interpolation point so far.

     XPT contains the coordinates of the current interpolation points.

     BMAT provides the last N columns of H.

     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.

     NDIM is the first dimension of BMAT and has the value NPT+N.

     KNEW is the index of the interpolation point that is going to be moved.

     DELTA is the current trust region bound.

     D will be set to the step from XOPT to the new point.

     ALPHA will be set to the KNEW-th diagonal element of the H matrix.

     HCOL, GC, GD, S and W will be used for working space.

     The step D is calculated in a way that attempts to maximize the modulus
     of LFUNC(XOPT+D), subject to the bound ||D|| .LE. DELTA, where LFUNC is
     the KNEW-th Lagrange function. */

  /* Constants. */
  const REAL half = 0.5;
  const REAL one = 1.0;
  const REAL zero = 0.0;
  const REAL twopi = 2.0*M_PI;

  /* Local variables. */
  REAL angle, cf1, cf2, cf3, cf4, cf5, cth, dd, delsq, denom, dhd, gg, scale,
    sp, ss, step, sth, sum, tau, taubeg, taumax, tauold, temp, tempa, tempb;
  INTEGER i, isave, iterc, iu, j, k, nptm;

  /* Parameter adjustments adjustments */
  xopt -= 1;
  xpt  -= 1 + npt;
  bmat -= 1 + ndim;
  zmat -= 1 + npt;
  d    -= 1;
  hcol -= 1;
  gc   -= 1;
  gd   -= 1;
  s    -= 1;
  w    -= 1;
#define XPT(a1,a2) xpt[(a2)*npt + a1]
#define BMAT(a1,a2) bmat[(a2)*ndim + a1]
#define ZMAT(a1,a2) zmat[(a2)*npt + a1]

  /* FIXME: Set uninitialized variables. */
  tempa = zero;
  tempb = zero;

  /* Set some constants. */
  delsq = delta*delta;
  nptm = npt - n - 1;

  /* Set the first NPT components of HCOL to the leading elements of the
     KNEW-th column of H. */
  iterc = 0;
  LOOP(k,npt) {
    hcol[k] = zero;
  }
  LOOP(j,nptm) {
    temp = ZMAT(knew, j);
    if (j < *idz) {
      temp = -temp;
    }
    LOOP(k,npt) {
      hcol[k] += temp*ZMAT(k,j);
    }
  }
  *alpha = hcol[knew];

  /* Set the unscaled initial direction D. Form the gradient of LFUNC at
     XOPT, and multiply D by the second derivative matrix of LFUNC. */
  dd = zero;
  LOOP(i,n) {
    d[i] = XPT(knew, i) - xopt[i];
    gc[i] = BMAT(knew, i);
    gd[i] = zero;
    dd += d[i]*d[i];
  }
  LOOP(k,npt) {
    temp = zero;
    sum = zero;
    LOOP(j,n) {
      temp += XPT(k,j)*xopt[j];
      sum += XPT(k,j)*d[j];
    }
    temp = hcol[k]*temp;
    sum = hcol[k]*sum;
    LOOP(i,n) {
      gc[i] += temp*XPT(k,i);
      gd[i] += sum*XPT(k,i);
    }
  }

  /* Scale D and GD, with a sign change if required. Set S to another
     vector in the initial two dimensional subspace. */
  gg = zero;
  sp = zero;
  dhd = zero;
  LOOP(i,n) {
    gg += gc[i]*gc[i];
    sp += d[i]*gc[i];
    dhd += d[i]*gd[i];
  }
  scale = delta/SQRT(dd);
  if (sp*dhd < zero) {
    scale = -scale;
  }
  temp = zero;
  if (sp*sp > dd*0.99*gg) {
    temp = one;
  }
  tau = scale*(ABS(sp) + half*scale*ABS(dhd));
  if (gg*delsq < tau*0.01*tau) {
    temp = one;
  }
  LOOP(i,n) {
    d[i] = scale*d[i];
    gd[i] = scale*gd[i];
    s[i] = gc[i] + temp*gd[i];
  }

  /* Begin the iteration by overwriting S with a vector that has the
     required length and direction, except that termination occurs if
     the given D and S are nearly parallel. */
 L80:
  ++iterc;
  dd = zero;
  sp = zero;
  ss = zero;
  LOOP(i,n) {
    dd += d[i]*d[i];
    sp += d[i]*s[i];
    ss += s[i]*s[i];
  }
  temp = dd*ss - sp*sp;
  if (temp <= dd*1e-8*ss) {
    goto L160;
  }
  denom = SQRT(temp);
  LOOP(i,n) {
    s[i] = (dd*s[i] - sp*d[i])/denom;
    w[i] = zero;
  }

  /* Calculate the coefficients of the objective function on the circle,
     beginning with the multiplication of S by the second derivative matrix. */
  LOOP(k,npt) {
    sum = zero;
    LOOP(j,n) {
      sum += XPT(k,j)*s[j];
    }
    sum = hcol[k]*sum;
    LOOP(i,n) {
      w[i] += sum*XPT(k,i);
    }
  }
  cf1 = zero;
  cf2 = zero;
  cf3 = zero;
  cf4 = zero;
  cf5 = zero;
  LOOP(i,n) {
    cf1 += s[i]*w[i];
    cf2 += d[i]*gc[i];
    cf3 += s[i]*gc[i];
    cf4 += d[i]*gd[i];
    cf5 += s[i]*gd[i];
  }
  cf1 = half*cf1;
  cf4 = half*cf4 - cf1;

  /* Seek the value of the angle that maximizes the modulus of TAU. */
  taubeg = cf1 + cf2 + cf4;
  taumax = taubeg;
  tauold = taubeg;
  isave = 0;
  iu = 49;
  temp = twopi/(REAL)(iu + 1);
  LOOP(i,iu) {
    angle = (REAL)i*temp;
    cth = COS(angle);
    sth = SIN(angle);
    tau = cf1 + (cf2 + cf4*cth)*cth + (cf3 + cf5*cth)*sth;
    if (ABS(tau) > ABS(taumax)) {
      taumax = tau;
      isave = i;
      tempa = tauold;
    } else if (i == isave + 1) {
      tempb = tau;
    }
    tauold = tau;
  }
  if (isave == 0) {
    tempa = tau;
  }
  if (isave == iu) {
    tempb = taubeg;
  }
  step = zero;
  if (tempa != tempb) {
    tempa -= taumax;
    tempb -= taumax;
    step = half*(tempa - tempb)/(tempa + tempb);
  }
  angle = temp*((REAL)isave + step);

  /* Calculate the new D and GD. Then test for convergence. */
  cth = COS(angle);
  sth = SIN(angle);
  tau = cf1 + (cf2 + cf4*cth)*cth + (cf3 + cf5*cth)*sth;
  LOOP(i,n) {
    d[i] = cth*d[i] + sth*s[i];
    gd[i] = cth*gd[i] + sth*w[i];
    s[i] = gc[i] + gd[i];
  }
  if (ABS(tau) <= ABS(taubeg)*1.1) {
    goto L160;
  }
  if (iterc < n) {
    goto L80;
  }
 L160:
  return;
} /* biglag */

#undef ZMAT
#undef BMAT
#undef XPT

static void
trsapp(const INTEGER n, const INTEGER npt, REAL* xopt,
       REAL* xpt, REAL* gq, REAL* hq, REAL* pq,
       const REAL delta, REAL* step, REAL* d, REAL* g,
       REAL* hd, REAL* hs, REAL* crvmin)
{
  /* N is the number of variables of a quadratic objective function, Q say.
     The arguments NPT, XOPT, XPT, GQ, HQ and PQ have their usual meanings, in
     order to define the current quadratic model Q.

     DELTA is the trust region radius, and has to be positive.

     STEP will be set to the calculated trial step.  The arrays D, G, HD and HS
     will be used for working space.

     CRVMIN will be set to the least curvature of H along the conjugate
     directions that occur, except that it is set to zero if STEP goes all the
     way to the trust region boundary.

     The calculation of STEP begins with the truncated conjugate gradient
     method. If the boundary of the trust region is reached, then further
     changes to STEP may be made, each one being in the 2D space spanned by the
     current STEP and the corresponding gradient of Q. Thus STEP should provide
     a substantial reduction to Q within the trust region. */

  /* Constants. */
  const REAL half = 0.5;
  const REAL zero = 0.0;
  const REAL twopi = 2.0*M_PI;

  /* Local variables. */
  REAL alpha, angle, angtest, bstep, cf, cth, dd, delsq, dg, dhd, dhs, ds,
    gg, ggbeg, ggsav, qadd, qbeg, qmin, qnew, qred, qsav, ratio, reduc,
    sg, sgk, shs, ss, sth, temp, tempa, tempb;
  INTEGER i, ih, isave, iterc, itermax, itersw, iu, j, k;

  /* Parameter adjustments to comply with FORTRAN indexing. */
  xopt -= 1;
  xpt  -= 1 + npt;
  gq   -= 1;
  hq   -= 1;
  pq   -= 1;
  step -= 1;
  d    -= 1;
  g    -= 1;
  hd   -= 1;
  hs   -= 1;
#define XPT(a1,a2) xpt[(a2)*npt + a1]

  /* FIXME: Set uninitialized variables. */
  tempa = zero;
  tempb = zero;
  shs = zero;
  sg = zero;
  qred = zero;
  ggbeg = zero;
  gg = zero;
  dd = zero;
  bstep = zero;

  /* Initialization, which includes setting HD to H times XOPT. */
  delsq = delta*delta;
  iterc = 0;
  itermax = n;
  itersw = itermax;
  LOOP(i,n) {
    d[i] = xopt[i];
  }
  goto L170;

  /* Prepare for the first line search. */
 L20:
  qred = zero;
  dd = zero;
  LOOP(i,n) {
    step[i] = zero;
    hs[i] = zero;
    g[i] = gq[i] + hd[i];
    d[i] = -g[i];
    dd += d[i]*d[i];
  }
  *crvmin = zero;
  if (dd == zero) {
    goto L160;
  }
  ds = zero;
  ss = zero;
  gg = dd;
  ggbeg = gg;

  /* Calculate the step to the trust region boundary and the product HD. */
 L40:
  ++iterc;
  temp = delsq - ss;
  bstep = temp/(ds + SQRT(ds*ds + dd*temp));
  goto L170;
 L50:
  dhd = zero;
  LOOP(j,n) {
    dhd += d[j]*hd[j];
  }

  /* Update CRVMIN and set the step-length ALPHA. */
  alpha = bstep;
  if (dhd > zero) {
    temp = dhd/dd;
    if (iterc == 1) {
      *crvmin = temp;
    } else {
      *crvmin = MIN(*crvmin,temp);
    }
    temp = gg/dhd;
    alpha = MIN(alpha,temp);
  }
  qadd = alpha*(gg - half*alpha*dhd);
  qred += qadd;

  /* Update STEP and HS. */
  ggsav = gg;
  gg = zero;
  LOOP(i,n) {
    step[i] += alpha*d[i];
    hs[i] += alpha*hd[i];
    temp = g[i] + hs[i];
    gg += temp*temp;
  }

  /* Begin another conjugate direction iteration if required. */
  if (alpha < bstep) {
    if (qadd <= qred*0.01) {
      goto L160;
    }
    if (gg <= ggbeg*1e-4) {
      goto L160;
    }
    if (iterc == itermax) {
      goto L160;
    }
    temp = gg/ggsav;
    dd = zero;
    ds = zero;
    ss = zero;
    LOOP(i,n) {
      d[i] = temp*d[i] - g[i] - hs[i];
      dd += d[i]*d[i];
      ds += d[i]*step[i];
      ss += step[i]*step[i];
    }
    if (ds <= zero) {
      goto L160;
    }
    if (ss < delsq) {
      goto L40;
    }
  }
  *crvmin = zero;
  itersw = iterc;

  /* Test whether an alternative iteration is required. */
 L90:
  if (gg <= ggbeg*1e-4) {
    goto L160;
  }
  sg = zero;
  shs = zero;
  LOOP(i,n) {
    sg += step[i]*g[i];
    shs += step[i]*hs[i];
  }
  sgk = sg + shs;
  angtest = sgk/SQRT(gg*delsq);
  if (angtest <= -0.99) {
    goto L160;
  }

  /* Begin the alternative iteration by calculating D and HD and some
     scalar products. */
  ++iterc;
  temp = SQRT(delsq*gg - sgk*sgk);
  tempa = delsq/temp;
  tempb = sgk/temp;
  LOOP(i,n) {
    d[i] = tempa*(g[i] + hs[i]) - tempb*step[i];
  }
  goto L170;
 L120:
  dg = zero;
  dhd = zero;
  dhs = zero;
  LOOP(i,n) {
    dg += d[i]*g[i];
    dhd += hd[i]*d[i];
    dhs += hd[i]*step[i];
  }

  /* Seek the value of the angle that minimizes Q. */
  cf = half*(shs - dhd);
  qbeg = sg + cf;
  qsav = qbeg;
  qmin = qbeg;
  isave = 0;
  iu = 49;
  temp = twopi/(REAL)(iu + 1);
  LOOP(i,iu) {
    angle = (REAL)i*temp;
    cth = COS(angle);
    sth = SIN(angle);
    qnew = (sg + cf*cth)*cth + (dg + dhs*cth)*sth;
    if (qnew < qmin) {
      qmin = qnew;
      isave = i;
      tempa = qsav;
    } else if (i == isave + 1) {
      tempb = qnew;
    }
    qsav = qnew;
  }
  if ((REAL)isave == zero) {
    tempa = qnew;
  }
  if (isave == iu) {
    tempb = qbeg;
  }
  angle = zero;
  if (tempa != tempb) {
    tempa -= qmin;
    tempb -= qmin;
    angle = half*(tempa - tempb)/(tempa + tempb);
  }
  angle = temp*((REAL)isave + angle);

  /* Calculate the new STEP and HS. Then test for convergence. */
  cth = COS(angle);
  sth = SIN(angle);
  reduc = qbeg - (sg + cf*cth)*cth - (dg + dhs*cth)*sth;
  gg = zero;
  LOOP(i,n) {
    step[i] = cth*step[i] + sth*d[i];
    hs[i] = cth*hs[i] + sth*hd[i];
    temp = g[i] + hs[i];
    gg += temp*temp;
  }
  qred += reduc;
  ratio = reduc/qred;
  if (iterc < itermax && ratio > 0.01) {
    goto L90;
  }
 L160:
  return;

  /* The following instructions act as a subroutine for setting the vector
     HD to the vector D multiplied by the second derivative matrix of Q.
     They are called from three different places, which are distinguished
     by the value of ITERC. */
 L170:
  LOOP(i,n) {
    hd[i] = zero;
  }
  LOOP(k,npt) {
    temp = zero;
    LOOP(j,n) {
      temp += XPT(k,j)*d[j];
    }
    temp *= pq[k];
    LOOP(i,n) {
      hd[i] += temp*XPT(k,i);
    }
  }
  ih = 0;
  LOOP(j,n) {
    LOOP(i,j) {
      ++ih;
      if (i < j) {
        hd[j] += hq[ih]*d[i];
      }
      hd[i] += hq[ih]*d[j];
    }
  }
  if (iterc == 0) {
    goto L20;
  }
  if (iterc <= itersw) {
    goto L50;
  }
  goto L120;
} /* trsapp */

#undef XPT

static void
update(const INTEGER n, const INTEGER npt, REAL* bmat,
        REAL* zmat, INTEGER* idz, const INTEGER ndim, REAL* vlag,
        const REAL beta, const INTEGER knew, REAL* w)
{
  /* The arrays BMAT and ZMAT with IDZ are updated, in order to shift the
     interpolation point that has index KNEW. On entry, VLAG contains the
     components of the vector Theta*Wcheck+e_b of the updating formula (6.11),
     and BETA holds the value of the parameter that has this name.  The vector
     W is used for working space. */

  /* Constants. */
  const REAL one = 1.0;
  const REAL zero = 0.0;

  /* Local variables. */
  REAL alpha, denom, scala, scalb, tau, tausq, temp, tempa, tempb;
  INTEGER i, iflag, j, ja, jb, jl, jp, nptm;

  /* Parameter adjustments to comply with FORTRAN indexing. */
  zmat -= 1 + npt;
  bmat -= 1 + ndim;
  vlag -= 1;
  w    -= 1;
#define BMAT(a1,a2) bmat[(a2)*ndim + a1]
#define ZMAT(a1,a2) zmat[(a2)*npt + a1]

  /* Set some constants. */
  nptm = npt - n - 1;

  /* Apply the rotations that put zeros in the KNEW-th row of ZMAT. */
  jl = 1;
  for (j = 2; j <= nptm; ++j) {
    if (j == *idz) {
      jl = *idz;
    } else if (ZMAT(knew, j) != zero) {
      tempa = ZMAT(knew,jl);
      tempb = ZMAT(knew,j);
      temp = SQRT(tempa*tempa + tempb*tempb);
      tempa /= temp;
      tempb /= temp;
      LOOP(i,npt) {
        temp = tempa*ZMAT(i,jl) + tempb*ZMAT(i,j);
        ZMAT(i,j) = tempa*ZMAT(i,j) - tempb*ZMAT(i,jl);
        ZMAT(i,jl) = temp;
      }
      ZMAT(knew, j) = zero;
    }
  }

  /* Put the first NPT components of the KNEW-th column of HLAG into W,
     and calculate the parameters of the updating formula. */
  tempa = ZMAT(knew,1);
  if (*idz >= 2) {
    tempa = -tempa;
  }
  if (jl > 1) {
    tempb = ZMAT(knew, jl);
  }
  LOOP(i,npt) {
    w[i] = tempa*ZMAT(i,1);
    if (jl > 1) {
      w[i] += tempb*ZMAT(i,jl);
    }
  }
  alpha = w[knew];
  tau = vlag[knew];
  tausq = tau*tau;
  denom = alpha*beta + tausq;
  vlag[knew] -= one;

  /* Complete the updating of ZMAT when there is only one nonzero element
     in the KNEW-th row of the new matrix ZMAT, but, if IFLAG is set to one,
     then the first column of ZMAT will be exchanged with another one later. */
  iflag = 0;
  if (jl == 1) {
    temp = SQRT((ABS(denom)));
    tempb = tempa/temp;
    tempa = tau/temp;
    LOOP(i,npt) {
      ZMAT(i,1) = tempa*ZMAT(i,1) - tempb*vlag[i];
    }
    if (*idz == 1 && temp < zero) {
      *idz = 2;
    }
    if (*idz >= 2 && temp >= zero) {
      iflag = 1;
    }
  } else {

    /* Complete the updating of ZMAT in the alternative case. */
    ja = 1;
    if (beta >= zero) {
      ja = jl;
    }
    jb = jl + 1 - ja;
    temp = ZMAT(knew, jb)/denom;
    tempa = temp*beta;
    tempb = temp*tau;
    temp = ZMAT(knew, ja);
    scala = one/SQRT(ABS(beta)*temp*temp + tausq);
    scalb = scala*SQRT((ABS(denom)));
    LOOP(i,npt) {
      ZMAT(i,ja) = scala*(tau*ZMAT(i,ja) - temp*vlag[i]);
      ZMAT(i,jb) = scalb*(ZMAT(i,jb) - tempa*w[i] - tempb*vlag[i]);
    }
    if (denom <= zero) {
      if (beta < zero) {
        ++(*idz);
      }
      if (beta >= zero) {
        iflag = 1;
      }
    }
  }

  /* IDZ is reduced in the following case, and usually the first column
     of ZMAT is exchanged with a later one. */
  if (iflag == 1) {
    --(*idz);
    LOOP(i,npt) {
      temp = ZMAT(i,1);
      ZMAT(i,1) = ZMAT(i, *idz);
      ZMAT(i, *idz) = temp;
    }
  }

  /* Finally, update the matrix BMAT. */
  LOOP(j,n) {
    jp = npt + j;
    w[jp] = BMAT(knew, j);
    tempa = (alpha*vlag[jp] - tau*w[jp])/denom;
    tempb = (-beta*w[jp] - tau*vlag[jp])/denom;
    LOOP(i,jp) {
      BMAT(i,j) = BMAT(i,j) + tempa*vlag[i] + tempb*w[i];
      if (i > npt) {
        BMAT(jp, i - npt) = BMAT(i,j);
      }
    }
  }
  return;
} /* update */

#undef ZMAT
#undef BMAT

static void
print_error(const char* reason)
{
  fprintf(stderr, "\n    Return from NEWUOA because %s.\n", reason);
}

static void
print_x(FILE* output, INTEGER n, const REAL x[], const REAL dx[])
{
  INTEGER i;
  for (i = 0; i < n; ++i) {
    fprintf(output, "%s%15.6E%s",
            ((i%5 == 0) ? "  " : ""),
            (double)(dx == NULL ? x[i] : (x[i] + dx[i])),
            ((i == n - 1 || i%5 == 4) ? "\n" : ""));
  }
}

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * fill-column: 79
 * coding: utf-8
 * End:
 */
