/*
 * bobyqa.c -
 *
 * Implementation of Mike Powell's BOBYQA algorithm for minimizing a function
 * of many variables.  The method is "derivatives free" (only the function
 * values are needed) and accounts for bound constraints on the variables.  The
 * algorithm is described in:
 *
 *   M.J.D. Powell, "The BOBYQA Algorithm for Bound Constrained Optimization
 *   Without Derivatives."  Technical report, Department of Applied Mathematics
 *   and Theoretical Physics, University of Cambridge (2009).
 *
 * The present code is based on the original FORTRAN version written by Mike
 * Powell who kindly provides his code on demand (at mjdp@cam.ac.uk) and has
 * been converted to C by É. Thiébaut.
 *
 * Copyright (c) 2009, Mike Powell (FORTRAN version).
 * Copyright (c) 2015, Éric Thiébaut (C version).
 *
 * Read the accompanying `LICENSE` file for details.
 */

#include <stdio.h>
#include <math.h>

#include "bobyqa.h"

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
print_x(FILE* output, const INTEGER n, const REAL x[], const REAL dx[]);

static int
bobyqb(const INTEGER n, const INTEGER npt,
       bobyqa_objfun* objfun, void* data,
       REAL* x, const REAL* xl, const REAL* xu,
       const REAL rhobeg, const REAL rhoend,
       const INTEGER iprint, const INTEGER maxfun,
       REAL* xbase, REAL* xpt, REAL* fval, REAL* xopt,
       REAL* gopt, REAL* hq, REAL* pq, REAL* bmat, REAL* zmat,
       const INTEGER ndim, REAL* sl, REAL* su, REAL* xnew,
       REAL* xalt, REAL* d, REAL* vlag, REAL* w);


static void
altmov(const INTEGER n, const INTEGER npt, REAL xpt[],
       REAL* xopt, REAL* bmat, REAL* zmat, const INTEGER ndim,
       REAL* sl, REAL* su, const INTEGER kopt, const INTEGER knew,
       const REAL adelt, REAL* xnew, REAL* xalt, REAL* alpha,
       REAL* cauchy, REAL* glag, REAL* hcol, REAL* w);

static void
prelim(const INTEGER n, const INTEGER npt,
       bobyqa_objfun* objfun, void* data,
       REAL* x, const REAL* xl, const REAL* xu,
       const REAL rhobeg, const INTEGER iprint,
       const INTEGER maxfun, REAL* xbase, REAL* xpt, REAL* fval,
       REAL* gopt, REAL* hq, REAL* pq, REAL* bmat,
       REAL* zmat, const INTEGER ndim, REAL* sl, REAL* su,
       INTEGER* nf, INTEGER* kopt);

static void
trsbox(const INTEGER n, const INTEGER npt, REAL* xpt,
       REAL* xopt, REAL* gopt, REAL* hq, REAL* pq,
       REAL* sl, REAL* su, const REAL delta, REAL* xnew,
       REAL* d, REAL* gnew, REAL* xbdi, REAL* s,
       REAL* hs, REAL* hred, REAL* dsq, REAL* crvmin);

static void
rescue(const INTEGER n, const INTEGER npt,
       bobyqa_objfun* objfun, void* data,
       const REAL* xl, const REAL* xu,
       const INTEGER iprint, const INTEGER maxfun,
       REAL* xbase, REAL* xpt, REAL* fval, REAL* xopt,
       REAL* gopt, REAL* hq, REAL* pq, REAL* bmat, REAL* zmat,
       const INTEGER ndim, REAL* sl, REAL* su, INTEGER* nf,
       const REAL delta, INTEGER* kopt, REAL* vlag, REAL* ptsaux,
       REAL* ptsid, REAL* w);

static void
update(const INTEGER n, const INTEGER npt, REAL* bmat,
       REAL* zmat, const INTEGER ndim, REAL* vlag, const REAL beta,
       const REAL denom, const INTEGER knew, REAL* w);

/*---------------------------------------------------------------------------*/
/* TESTING */

static REAL
objfun_test(const INTEGER n, const REAL* x, void* data);

#ifdef TESTING

int
main(int argc, char* argv[])
{
  bobyqa_test();
  return 0;
}

#ifdef FORTRAN_NAME
int
FORTRAN_NAME(calfun,CALFUN)(const INTEGER* n, REAL* x, REAL* f)
{
  *f = objfun_test(*n, x, NULL);
  return 0;
}
#endif /* FORTRAN_NAME */

#endif

void
bobyqa_test(void)
{
  /* Constants. */
  const REAL twopi = 2.0*M_PI;

  /* Local variables. */
  REAL bdl, bdu, rhobeg, rhoend, temp;
  REAL w[500000], x[100], xl[100], xu[100];
  INTEGER i, iprint, j, jcase, m, maxfun;

  bdl = -1.0;
  bdu =  1.0;
  iprint = 2;
  maxfun = 500000;
  rhobeg = 0.1;
  rhoend = 1e-6;
  for (m = 5; m <= 10; m += m) {
    REAL q = twopi/(REAL)m;
    INTEGER n = 2*m;
    LOOP(i,n) {
      xl[i - 1] = bdl;
      xu[i - 1] = bdu;
    }
    for (jcase = 1; jcase <= 2; ++jcase) {
      INTEGER npt = n + 6;
      if (jcase == 2) {
        npt = 2*n + 1;
      }
      fprintf(OUTPUT,
              "\n\n     2D output with M =%4ld,  N =%4ld  and  NPT =%4ld\n",
              (long)m, (long)n, (long)npt);
      LOOP(j,m) {
        temp = ((REAL)j)*q;
        x[2*j - 2] = COS(temp);
        x[2*j - 1] = SIN(temp);
      }
      bobyqa(n, npt, objfun_test, NULL, x, xl, xu, rhobeg, rhoend,
             iprint, maxfun, w);
    }
  }
}

static REAL
objfun_test(const INTEGER n, const REAL* x, void* data)
{
  INTEGER i;
  REAL f, temp, tempa, tempb;

  f = 0.0;
  for (i = 4; i <= n; i += 2) {
    INTEGER j, j1 = i - 2;
    for (j = 2; j <= j1; j += 2) {
      tempa = x[i - 2] - x[j - 2];
      tempb = x[i - 1] - x[j - 1];
      temp = tempa*tempa + tempb*tempb;
      temp = MAX(temp,1e-6);
      f += 1.0/SQRT(temp);
    }
  }
  return f;
}

/*---------------------------------------------------------------------------*/
/* BOBYQA DRIVER ROUTINES */

int
bobyqa(const INTEGER n, const INTEGER npt,
       bobyqa_objfun* objfun, void* data,
       REAL* x, const REAL* xl, const REAL* xu,
       const REAL rhobeg, const REAL rhoend,
       const INTEGER iprint, const INTEGER maxfun, REAL* w)
{
  /* Constants. */
  const REAL zero = 0.0;

  /* Local variables. */
  REAL temp, tempa, tempb;
  INTEGER ibmat, id, ifv, igo, ihq, ipq, isl, isu, ivl, iw, ixa, ixb, ixn,
    ixo, ixp, izmat, j, jsl, jsu, ndim, np;

  /* Parameter adjustments to comply with FORTRAN indexing. */
  w  -= 1;
  xu -= 1;
  xl -= 1;
  x  -= 1;

  /* Return if the value of NPT is unacceptable. */
  np = n + 1;
  if (npt < n + 2 || npt > (n + 2)*np/2) {
    print_error("NPT is not in the required interval");
    return BOBYQA_BAD_NPT;
  }

  /* Partition the working space array, so that different parts of it can
     be treated separately during the calculation of BOBYQB.  The partition
     requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the
     space that is taken by the last array in the argument list of
     BOBYQB. */
  ndim = npt + n;
  ixb = 1;
  ixp = ixb + n;
  ifv = ixp + n*npt;
  ixo = ifv + npt;
  igo = ixo + n;
  ihq = igo + n;
  ipq = ihq + n*np/2;
  ibmat = ipq + npt;
  izmat = ibmat + ndim*n;
  isl = izmat + npt*(npt - np);
  isu = isl + n;
  ixn = isu + n;
  ixa = ixn + n;
  id = ixa + n;
  ivl = id + n;
  iw = ivl + ndim;

  /* Return if there is insufficient space between the bounds.  Modify the
     initial X if necessary in order to avoid conflicts between the bounds and
     the construction of the first quadratic model.  The lower and upper bounds
     on moves from the updated X are set now, in the ISL and ISU partitions of
     W, in order to provide useful and exact information about components of X
     that become within distance RHOBEG from their bounds. */
  LOOP(j,n) {
    temp = xu[j] - xl[j];
    if (temp < rhobeg + rhobeg) {
      print_error("one of the differences XU(I)-XL(I) is less than 2*RHOBEG");
      return BOBYQA_TOO_CLOSE;
    }
    jsl = isl + j - 1;
    jsu = jsl + n;
    w[jsl] = xl[j] - x[j];
    w[jsu] = xu[j] - x[j];
    if (w[jsl] >= -rhobeg) {
      if (w[jsl] >= zero) {
        x[j] = xl[j];
        w[jsl] = zero;
        w[jsu] = temp;
      } else {
        x[j] = xl[j] + rhobeg;
        w[jsl] = -rhobeg;
        temp = xu[j] - x[j];
        w[jsu] = MAX(temp,rhobeg);
      }
    } else if (w[jsu] <= rhobeg) {
      if (w[jsu] <= zero) {
        x[j] = xu[j];
        w[jsl] = -temp;
        w[jsu] = zero;
      } else {
        x[j] = xu[j] - rhobeg;
        tempa = xl[j] - x[j];
        tempb = -rhobeg;
        w[jsl] = MIN(tempa,tempb);
        w[jsu] = rhobeg;
      }
    }
  }

  /* Make the call of BOBYQB. */
  return bobyqb(n, npt, objfun, data, &x[1], &xl[1], &xu[1],
                rhobeg, rhoend, iprint, maxfun,
                &w[ixb], &w[ixp], &w[ifv], &w[ixo], &w[igo],
                &w[ihq], &w[ipq], &w[ibmat], &w[izmat],
                ndim, &w[isl], &w[isu], &w[ixn], &w[ixa],
                &w[id], &w[ivl], &w[iw]);
} /* bobyqa */

/*---------------------------------------------------------------------------*/
/* FORTRAN SUPPORT */

#ifdef FORTRAN_NAME

REAL
bobyqa_calfun_wrapper(const INTEGER n, const REAL* x, void* data)
{
  REAL f;
  FORTRAN_NAME(calfun,CALFUN)(&n, (REAL*)x, &f);
  return f;
}

int
FORTRAN_NAME(bobyqa,BOBYQA)(const INTEGER* n, const INTEGER* npt,
                            REAL* x, const REAL* xl, const REAL* xu,
                            const REAL* rhobeg, const REAL* rhoend,
                            const INTEGER* iprint, const INTEGER* maxfun,
                            REAL* w)
{
  bobyqa(*n, *npt, bobyqa_calfun_wrapper, NULL, x, xl, xu,
         *rhobeg, *rhoend, *iprint, *maxfun, w);
  return 0;
}

#endif /* FORTRAN_NAME */

/*---------------------------------------------------------------------------*/
/* BOBYQA SUBROUTINES */

static int
bobyqb(const INTEGER n, const INTEGER npt,
       bobyqa_objfun* objfun, void* data,
       REAL* x, const REAL* xl, const REAL* xu,
       const REAL rhobeg, const REAL rhoend,
       const INTEGER iprint, const INTEGER maxfun,
       REAL* xbase, REAL* xpt, REAL* fval, REAL* xopt,
       REAL* gopt, REAL* hq, REAL* pq, REAL* bmat, REAL* zmat,
       const INTEGER ndim, REAL* sl, REAL* su, REAL* xnew,
       REAL* xalt, REAL* d, REAL* vlag, REAL* w)
{
  /* The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT and MAXFUN are
     identical to the corresponding arguments in SUBROUTINE BOBYQA.

     XBASE holds a shift of origin that should reduce the contributions from
     rounding errors to values of the model and Lagrange functions.

     XPT is a two-dimensional array that holds the coordinates of the
     interpolation points relative to XBASE.

     FVAL holds the values of F at the interpolation points.

     XOPT is set to the displacement from XBASE of the trust region centre.

     GOPT holds the gradient of the quadratic model at XBASE+XOPT.

     HQ holds the explicit second derivatives of the quadratic model.

     PQ contains the parameters of the implicit second derivatives of the
     quadratic model.

     BMAT holds the last N columns of H.

     ZMAT holds the factorization of the leading NPT by NPT submatrix of H,
     this factorization being ZMAT times ZMAT^T, which provides both the
     correct rank and positive semi-definiteness.

     NDIM is the first dimension of BMAT and has the value NPT+N.

     SL and SU hold the differences XL-XBASE and XU-XBASE, respectively.  All
     the components of every XOPT are going to satisfy the bounds
     SL(I) <= XOPT(I) <= SU(I), with appropriate equalities when XOPT is on a
     constraint boundary.

     XNEW is chosen by SUBROUTINE TRSBOX or ALTMOV.  Usually XBASE+XNEW is the
     vector of variables for the next call of OBJFUN.  XNEW also satisfies the
     SL and SU constraints in the way that has just been mentioned.

     XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW in
     order to increase the denominator in the updating of UPDATE.

     D is reserved for a trial step from XOPT, which is usually XNEW-XOPT.

     VLAG contains the values of the Lagrange functions at a new point X.  They
     are part of a product that requires VLAG to be of length NDIM.

     W is a one-dimensional array that is used for working space.  Its length
     must be at least 3*NDIM = 3*(NPT+N). */

  /* Constants. */
  const REAL half = 0.5;
  const REAL one = 1.0;
  const REAL ten = 10.0;
  const REAL tenth = 0.1;
  const REAL two = 2.0;
  const REAL zero = 0.0;

  /* Local variables. */
  REAL adelt, alpha, bdtest, bdtol, beta, biglsq, bsum, cauchy, crvmin,
    curv, delsq, delta, den, denom, densav, diff, diffa, diffb, diffc,
    dist, distsq, dnorm, dsq, dx, errbig, f, fopt, fracsq, frhosq, fsave,
    gisq, gqsq, hdiag, pqold, ratio, rho, scaden, sum, suma, sumb, sumpq,
    sumw, sumz, temp, tempa, tempb, vquad, xoptsq;
  INTEGER i, ih, ip, itest, j, jj, jp, k, kbase, knew, kopt, ksav, nf,
    nfsav, nh, np, nptm, nresc, ntrits;
  int status = BOBYQA_SUCCESS;
  const char* reason = NULL;

  /* Parameter adjustments to comply with FORTRAN indexing. */
  x     -= 1;
  xl    -= 1;
  xu    -= 1;
  xbase -= 1;
  xpt   -= 1 + npt;
  fval  -= 1;
  xopt  -= 1;
  gopt  -= 1;
  hq    -= 1;
  pq    -= 1;
  bmat  -= 1 + ndim;
  zmat  -= 1 + npt;
  sl    -= 1;
  su    -= 1;
  xnew  -= 1;
  xalt  -= 1;
  d     -= 1;
  vlag  -= 1;
  w     -= 1;
#define XPT(a1,a2) xpt[(a2)*npt + a1]
#define BMAT(a1,a2) bmat[(a2)*ndim + a1]
#define ZMAT(a1,a2) zmat[(a2)*npt + a1]

  /* Set uninitialized variables to avoid compiler warnings. */
  adelt = zero;
  alpha = zero;
  cauchy = zero;
  denom = zero;
  diff = zero;
  diffc = zero;
  f = zero;
  knew = 0;

  /* Set some constants. */
  np = n + 1;
  nptm = npt - np;
  nh = n*np/2;

  /* The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
     BMAT and ZMAT for the first iteration, with the corresponding values of of
     NF and KOPT, which are the number of calls of OBJFUN so far and the index
     of the interpolation point at the trust region centre.  Then the initial
     XOPT is set too.  The branch to label 720 occurs if MAXFUN is less than
     NPT.  GOPT will be updated if KOPT is different from KBASE. */
  prelim(n, npt, objfun, data, &x[1], &xl[1], &xu[1], rhobeg, iprint, maxfun,
         &xbase[1], &XPT(1,1), &fval[1], &gopt[1], &hq[1], &pq[1], &BMAT(1,1),
         &ZMAT(1,1), ndim, &sl[1], &su[1], &nf, &kopt);
  xoptsq = zero;
  LOOP(i,n) {
    xopt[i] = XPT(kopt,i);
    xoptsq += xopt[i]*xopt[i];
  }
  fsave = fval[1];
  if (nf < npt) {
    goto too_many_evaluations;
  }
  kbase = 1;

  /* Complete the settings that are required for the iterative procedure. */
  rho = rhobeg;
  delta = rho;
  nresc = nf;
  ntrits = 0;
  diffa = zero;
  diffb = zero;
  itest = 0;
  nfsav = nf;

  /* Update GOPT if necessary before the first iteration and after each
     call of RESCUE that makes a call of OBJFUN. */
 L20:
  if (kopt != kbase) {
    ih = 0;
    LOOP(j,n) {
      LOOP(i,j) {
        ++ih;
        if (i < j) {
          gopt[j] += hq[ih]*xopt[i];
        }
        gopt[i] += hq[ih]*xopt[j];
      }
    }
    if (nf > npt) {
      LOOP(k,npt) {
        temp = zero;
        LOOP(j,n) {
          temp += XPT(k,j)*xopt[j];
        }
        temp = pq[k]*temp;
        LOOP(i,n) {
          gopt[i] += temp*XPT(k,i);
        }
      }
    }
  }

  /* Generate the next point in the trust region that provides a small value
     of the quadratic model subject to the constraints on the variables.
     The integer NTRITS is set to the number "trust region" iterations that
     have occurred since the last "alternative" iteration.  If the length
     of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to
     label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW. */
 L60:
  trsbox(n, npt, &XPT(1,1), &xopt[1], &gopt[1], &hq[1], &pq[1],
         &sl[1], &su[1], delta, &xnew[1], &d[1], &w[1], &w[np],
         &w[np + n], &w[np + 2*n], &w[np + n*3], &dsq, &crvmin);
  dnorm = SQRT(dsq);
  dnorm = MIN(dnorm,delta);
  if (dnorm < half*rho) {
    ntrits = -1;
    tempa = ten*rho;
    distsq = tempa*tempa;
    if (nf <= nfsav + 2) {
      goto L650;
    }

    /* The following choice between labels 650 and 680 depends on whether or
       not our work with the current RHO seems to be complete.  Either RHO is
       decreased or termination occurs if the errors in the quadratic model at
       the last three interpolation points compare favourably with predictions
       of likely improvements to the model within distance HALF*RHO of XOPT. */
    errbig = MAX(diffa,diffb);
    errbig = MAX(errbig,diffc);
    frhosq = rho*0.125*rho;
    if (crvmin > zero && errbig > frhosq*crvmin) {
      goto L650;
    }
    bdtol = errbig/rho;
    LOOP(j,n) {
      bdtest = bdtol;
      if (xnew[j] == sl[j]) {
        bdtest = w[j];
      }
      if (xnew[j] == su[j]) {
        bdtest = -w[j];
      }
      if (bdtest < bdtol) {
        curv = hq[(j + j*j)/2];
        LOOP(k,npt) {
          curv += pq[k]*(XPT(k,j)*XPT(k,j));
        }
        bdtest += half*curv*rho;
        if (bdtest < bdtol) {
          goto L650;
        }
      }
    }
    goto L680;
  }
  ++ntrits;

  /* Severe cancellation is likely to occur if XOPT is too far from XBASE.
     If the following test holds, then XBASE is shifted so that XOPT becomes
     zero.  The appropriate changes are made to BMAT and to the second
     derivatives of the current model, beginning with the changes to BMAT
     that do not depend on ZMAT.  VLAG is used temporarily for working space. */
 L90:
  if (dsq <= xoptsq*0.001) {
    fracsq = xoptsq*0.25;
    sumpq = zero;
    LOOP(k,npt) {
      sumpq += pq[k];
      sum = -half*xoptsq;
      LOOP(i,n) {
        sum += XPT(k,i)*xopt[i];
      }
      w[npt + k] = sum;
      temp = fracsq - half*sum;
      LOOP(i,n) {
        w[i] = BMAT(k,i);
        vlag[i] = sum*XPT(k,i) + temp*xopt[i];
        ip = npt + i;
        LOOP(j,i) {
          BMAT(ip,j) = BMAT(ip,j) + w[i]*vlag[j] + vlag[i]*w[j];
        }
      }
    }

    /* Then the revisions of BMAT that depend on ZMAT are calculated. */
    LOOP(jj,nptm) {
      sumz = zero;
      sumw = zero;
      LOOP(k,npt) {
        sumz += ZMAT(k,jj);
        vlag[k] = w[npt + k]*ZMAT(k,jj);
        sumw += vlag[k];
      }
      LOOP(j,n) {
        sum = (fracsq*sumz - half*sumw)*xopt[j];
        LOOP(k,npt) {
          sum += vlag[k]*XPT(k,j);
        }
        w[j] = sum;
        LOOP(k,npt) {
          BMAT(k,j) = BMAT(k,j) + sum*ZMAT(k,jj);
        }
      }
      LOOP(i,n) {
        ip = i + npt;
        temp = w[i];
        LOOP(j,i) {
          BMAT(ip,j) = BMAT(ip,j) + temp*w[j];
        }
      }
    }

    /* The following instructions complete the shift, including the changes to
       the second derivative parameters of the quadratic model. */
    ih = 0;
    LOOP(j,n) {
      w[j] = -half*sumpq*xopt[j];
      LOOP(k,npt) {
        w[j] += pq[k]*XPT(k,j);
        XPT(k,j) = XPT(k,j) - xopt[j];
      }
      LOOP(i,j) {
        ++ih;
        hq[ih] = hq[ih] + w[i]*xopt[j] + xopt[i]*w[j];
        BMAT(npt + i, j) = BMAT(npt + j, i);
      }
    }
    LOOP(i,n) {
      xbase[i] += xopt[i];
      xnew[i] -= xopt[i];
      sl[i] -= xopt[i];
      su[i] -= xopt[i];
      xopt[i] = zero;
    }
    xoptsq = zero;
  }
  if (ntrits == 0) {
    goto L210;
  }
  goto L230;

  /* XBASE is also moved to XOPT by a call of RESCUE.  This calculation is more
     expensive than the previous shift, because new matrices BMAT and ZMAT are
     generated from scratch, which may include the replacement of interpolation
     points whose positions seem to be causing near linear dependence in the
     interpolation conditions.  Therefore RESCUE is called only if rounding
     errors have reduced by at least a factor of two the denominator of the
     formula for updating the H matrix.  It provides a useful safeguard, but is
     not invoked in most applications of BOBYQA. */
 L190:
  nfsav = nf;
  kbase = kopt;
  rescue(n, npt, objfun, data, &xl[1], &xu[1], iprint, maxfun, &xbase[1],
         &XPT(1,1), &fval[1], &xopt[1], &gopt[1], &hq[1],
         &pq[1], &BMAT(1,1), &ZMAT(1,1), ndim, &sl[1], &su[1],
         &nf, delta, &kopt, &vlag[1], &w[1], &w[n + np], &w[ndim + np]);

  /* XOPT is updated now in case the branch below to label 720 is taken.  Any
     updating of GOPT occurs after the branch below to label 20, which leads to
     a trust region iteration as does the branch to label 60. */
  xoptsq = zero;
  if (kopt != kbase) {
    LOOP(i,n) {
      xopt[i] = XPT(kopt,i);
      xoptsq += xopt[i]*xopt[i];
    }
  }
  if (nf < 0) {
    nf = maxfun;
    goto too_many_evaluations;
  }
  nresc = nf;
  if (nfsav < nf) {
    nfsav = nf;
    goto L20;
  }
  if (ntrits > 0) {
    goto L60;
  }

  /* Pick two alternative vectors of variables, relative to XBASE, that are
     suitable as new positions of the KNEW-th interpolation point.  Firstly,
     XNEW is set to the point on a line through XOPT and another interpolation
     point that minimizes the predicted value of the next denominator, subject
     to ||XNEW - XOPT|| <= ADELT and to the SL and SU bounds.  Secondly, XALT
     is set to the best feasible point on a constrained version of the Cauchy
     step of the KNEW-th Lagrange function, the corresponding value of the
     square of this function being returned in CAUCHY.  The choice between
     these alternatives is going to be made when the denominator is
     calculated. */
 L210:
  altmov(n, npt, &XPT(1,1), &xopt[1], &BMAT(1,1), &ZMAT(1,1), ndim,
         &sl[1], &su[1], kopt, knew, adelt, &xnew[1], &xalt[1],
         &alpha, &cauchy, &w[1], &w[np], &w[ndim + 1]);
  LOOP(i,n) {
    d[i] = xnew[i] - xopt[i];
  }

  /* Calculate VLAG and BETA for the current choice of D.  The scalar product
     of D with XPT(K,.) is going to be held in W(NPT+K) for use when VQUAD is
     calculated. */
 L230:
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
    w[npt + k] = suma;
  }
  beta = zero;
  LOOP(jj,nptm) {
    sum = zero;
    LOOP(k,npt) {
      sum += ZMAT(k,jj)*w[k];
    }
    beta -= sum*sum;
    LOOP(k,npt) {
      vlag[k] += sum*ZMAT(k,jj);
    }
  }
  dsq = zero;
  bsum = zero;
  dx = zero;
  LOOP(j,n) {
    dsq += d[j]*d[j];
    sum = zero;
    LOOP(k,npt) {
      sum += w[k]*BMAT(k,j);
    }
    bsum += sum*d[j];
    jp = npt + j;
    LOOP(i,n) {
      sum += BMAT(jp,i)*d[i];
    }
    vlag[jp] = sum;
    bsum += sum*d[j];
    dx += d[j]*xopt[j];
  }
  beta = dx*dx + dsq*(xoptsq + dx + dx + half*dsq) + beta - bsum;
  vlag[kopt] += one;

  /* If NTRITS is zero, the denominator may be increased by replacing the step
     D of ALTMOV by a Cauchy step.  Then RESCUE may be called if rounding
     errors have damaged the chosen denominator. */
  if (ntrits == 0) {
    denom = vlag[knew]*vlag[knew] + alpha*beta;
    if (denom < cauchy && cauchy > zero) {
      LOOP(i,n) {
        xnew[i] = xalt[i];
        d[i] = xnew[i] - xopt[i];
      }
      cauchy = zero;
      goto L230;
    }
    if (denom <= half*(vlag[knew]*vlag[knew])) {
      if (nf > nresc) {
        goto L190;
      }
      goto cancellation_of_denominator;
    }

    /* Alternatively, if NTRITS is positive, then set KNEW to the index of the
       next interpolation point to be deleted to make room for a trust region
       step.  Again RESCUE may be called if rounding errors have damaged the
       chosen denominator, which is the reason for attempting to select KNEW
       before calculating the next value of the objective function. */
  } else {
    delsq = delta*delta;
    scaden = zero;
    biglsq = zero;
    knew = 0;
    LOOP(k,npt) {
      if (k == kopt) continue;
      hdiag = zero;
      LOOP(jj,nptm) {
        hdiag += ZMAT(k,jj)*ZMAT(k,jj);
      }
      den = beta*hdiag + vlag[k]*vlag[k];
      distsq = zero;
      LOOP(j,n) {
        tempa = XPT(k,j) - xopt[j];
        distsq += tempa*tempa;
      }
      temp = distsq/delsq;
      temp = temp*temp;
      temp = MAX(one,temp);
      if (temp*den > scaden) {
        scaden = temp*den;
        knew = k;
        denom = den;
      }
      temp *= vlag[k]*vlag[k];
      biglsq = MAX(biglsq,temp);
    }
    if (scaden <= half*biglsq) {
      if (nf > nresc) {
        goto L190;
      }
      goto cancellation_of_denominator;
    }
  }

  /* Put the variables for the next calculation of the objective function in
     XNEW, with any adjustments for the bounds.

     Calculate the value of the objective function at XBASE+XNEW, unless the
     limit on the number of calculations of F has been reached. */
 L360:
  LOOP(i,n) {
    tempa = xbase[i] + xnew[i];
    tempa = MAX(tempa,xl[i]);
    x[i] = MIN(tempa,xu[i]);
    if (xnew[i] == sl[i]) {
      x[i] = xl[i];
    }
    if (xnew[i] == su[i]) {
      x[i] = xu[i];
    }
  }
  if (nf >= maxfun) {
    goto too_many_evaluations;
  }
  ++nf;
  f = objfun(n, &x[1], data);
  if (iprint == 3) {
    fprintf(OUTPUT,
            "    Function number%6ld    F =%18.10E"
            "    The corresponding X is:\n",
            (long)nf, (double)f);
    print_x(OUTPUT, n, &x[1], NULL);
  }
  if (ntrits == -1) {
    fsave = f;
    goto done;
  }

  /* Use the quadratic model to predict the change in F due to the step D, and
     set DIFF to the error of this prediction. */
  fopt = fval[kopt];
  vquad = zero;
  ih = 0;
  LOOP(j,n) {
    vquad += d[j]*gopt[j];
    LOOP(i,j) {
      ++ih;
      temp = d[i]*d[j];
      if (i == j) {
        temp = half*temp;
      }
      vquad += hq[ih]*temp;
    }
  }
  LOOP(k,npt) {
    vquad += half*pq[k]*(w[npt + k]*w[npt + k]);
  }
  diff = f - fopt - vquad;
  diffc = diffb;
  diffb = diffa;
  diffa = ABS(diff);
  if (dnorm > rho) {
    nfsav = nf;
  }

  /* Pick the next value of DELTA after a trust region step. */
  if (ntrits > 0) {
    if (vquad >= zero) {
      goto step_failed;
    }
    ratio = (f - fopt)/vquad;
    if (ratio <= tenth) {
      delta *= half;
      delta = MIN(delta,dnorm);
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

    /* Recalculate KNEW and DENOM if the new F is less than FOPT. */
    if (f < fopt) {
      ksav = knew;
      densav = denom;
      delsq = delta*delta;
      scaden = zero;
      biglsq = zero;
      knew = 0;
      LOOP(k,npt) {
        hdiag = zero;
        LOOP(jj,nptm) {
          hdiag += ZMAT(k,jj)*ZMAT(k,jj);
        }
        den = beta*hdiag + vlag[k]*vlag[k];
        distsq = zero;
        LOOP(j,n) {
          temp = XPT(k,j) - xnew[j];
          distsq += temp*temp;
        }
        temp = distsq/delsq;
        temp = temp*temp;
        temp = MAX(one,temp);
        if (temp*den > scaden) {
          scaden = temp*den;
          knew = k;
          denom = den;
        }
        temp *= (vlag[k]*vlag[k]);
        biglsq = MAX(biglsq,temp);
      }
      if (scaden <= half*biglsq) {
        knew = ksav;
        denom = densav;
      }
    }
  }

  /* Update BMAT and ZMAT, so that the KNEW-th interpolation point can be
     moved.  Also update the second derivative terms of the model. */
  update(n, npt, &BMAT(1,1), &ZMAT(1,1), ndim, &vlag[1],
         beta, denom, knew, &w[1]);
  ih = 0;
  pqold = pq[knew];
  pq[knew] = zero;
  LOOP(i,n) {
    temp = pqold*XPT(knew,i);
    LOOP(j,i) {
      ++ih;
      hq[ih] += temp*XPT(knew,j);
    }
  }
  LOOP(jj,nptm) {
    temp = diff*ZMAT(knew,jj);
    LOOP(k,npt) {
      pq[k] += temp*ZMAT(k,jj);
    }
  }

  /* Include the new interpolation point, and make the changes to GOPT at the
     old XOPT that are caused by the updating of the quadratic model. */
  fval[knew] = f;
  LOOP(i,n) {
    XPT(knew,i) = xnew[i];
    w[i] = BMAT(knew,i);
  }
  LOOP(k,npt) {
    suma = zero;
    LOOP(jj,nptm) {
      suma += ZMAT(knew,jj)*ZMAT(k,jj);
    }
    sumb = zero;
    LOOP(j,n) {
      sumb += XPT(k,j)*xopt[j];
    }
    temp = suma*sumb;
    LOOP(i,n) {
      w[i] += temp*XPT(k,i);
    }
  }
  LOOP(i,n) {
    gopt[i] += diff*w[i];
  }

  /* Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT. */
  if (f < fopt) {
    kopt = knew;
    xoptsq = zero;
    ih = 0;
    LOOP(j,n) {
      xopt[j] = xnew[j];
      xoptsq += xopt[j]*xopt[j];
      LOOP(i,j) {
        ++ih;
        if (i < j) {
          gopt[j] += hq[ih]*d[i];
        }
        gopt[i] += hq[ih]*d[j];
      }
    }
    LOOP(k,npt) {
      temp = zero;
      LOOP(j,n) {
        temp += XPT(k,j)*d[j];
      }
      temp = pq[k]*temp;
      LOOP(i,n) {
        gopt[i] += temp*XPT(k,i);
      }
    }
  }

  /* Calculate the parameters of the least Frobenius norm interpolant to
     the current data, the gradient of this interpolant at XOPT being put
     into VLAG(NPT+I), I=1,2,...,N. */
  if (ntrits > 0) {
    LOOP(k,npt) {
      vlag[k] = fval[k] - fval[kopt];
      w[k] = zero;
    }
    LOOP(j,nptm) {
      sum = zero;
      LOOP(k,npt) {
        sum += ZMAT(k,j)*vlag[k];
      }
      LOOP(k,npt) {
        w[k] += sum*ZMAT(k,j);
      }
    }
    LOOP(k,npt) {
      sum = zero;
      LOOP(j,n) {
        sum += XPT(k,j)*xopt[j];
      }
      w[k + npt] = w[k];
      w[k] = sum*w[k];
    }
    gqsq = zero;
    gisq = zero;
    LOOP(i,n) {
      sum = zero;
      LOOP(k,npt) {
        sum = sum + BMAT(k,i)*vlag[k] + XPT(k,i)*w[k];
      }
      if (xopt[i] == sl[i]) {
        tempa = MIN(zero,gopt[i]);
        gqsq += tempa*tempa;
        tempa = MIN(zero,sum);
        gisq += tempa*tempa;
      } else if (xopt[i] == su[i]) {
        tempa = MAX(zero,gopt[i]);
        gqsq += tempa*tempa;
        tempa = MAX(zero,sum);
        gisq += tempa*tempa;
      } else {
        gqsq += gopt[i]*gopt[i];
        gisq += sum*sum;
      }
      vlag[npt + i] = sum;
    }

    /* Test whether to replace the new quadratic model by the least Frobenius
       norm interpolant, making the replacement if the test is satisfied. */
    ++itest;
    if (gqsq < ten*gisq) {
      itest = 0;
    }
    if (itest >= 3) {
      INTEGER i1 = MAX(npt,nh);
      LOOP(i,i1) {
        if (i <= n) {
          gopt[i] = vlag[npt + i];
        }
        if (i <= npt) {
          pq[i] = w[npt + i];
        }
        if (i <= nh) {
          hq[i] = zero;
        }
        itest = 0;
      }
    }
  }

  /* If a trust region step has provided a sufficient decrease in F, then
     branch for another trust region calculation.  The case NTRITS=0 occurs
     when the new interpolation point was reached by an alternative step. */
  if (ntrits == 0) {
    goto L60;
  }
  if (f <= fopt + tenth*vquad) {
    goto L60;
  }

  /* Alternatively, find out if the interpolation points are close enough to
     the best point so far. */
  tempa = two*delta;
  tempb = ten*rho;
  tempa = tempa*tempa;
  tempb = tempb*tempb;
  distsq = MAX(tempa,tempb);
 L650:
  knew = 0;
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

  /* If KNEW is positive, then ALTMOV finds alternative new positions for the
     KNEW-th interpolation point within distance ADELT of XOPT.  It is reached
     via label 90.  Otherwise, there is a branch to label 60 for another trust
     region iteration, unless the calculations with the current RHO are
     complete. */
  if (knew > 0) {
    dist = SQRT(distsq);
    if (ntrits == -1) {
      tempa = tenth*delta;
      tempb = half*dist;
      delta = MIN(tempa,tempb);
      if (delta <= rho*1.5) {
        delta = rho;
      }
    }
    ntrits = 0;
    adelt = tenth*dist;
    adelt = MIN(adelt,delta);
    adelt = MAX(adelt,rho);
    dsq = adelt*adelt;
    goto L90;
  }
  if (ntrits == -1) {
    goto L680;
  }
  if (ratio > zero) {
    goto L60;
  }
  if (MAX(delta,dnorm) > rho) {
    goto L60;
  }

  /* The calculations with the current value of RHO are complete.  Pick the
     next values of RHO and DELTA. */
 L680:
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
              (double)rho, (long)nf, (double)fval[kopt]);
      print_x(OUTPUT, n, &xbase[1], &xopt[1]);
    }
    ntrits = 0;
    nfsav = nf;
    goto L60;
  }

  /* Return from the calculation, after another Newton-Raphson step, if it is
     too short to have been tried before. */
  if (ntrits == -1) {
    goto L360;
  }
 done:
  if (fval[kopt] <= fsave) {
    LOOP(i,n) {
      tempa = xbase[i] + xopt[i];
      tempa = MAX(tempa,xl[i]);
      x[i] = MIN(tempa,xu[i]);
      if (xopt[i] == sl[i]) {
        x[i] = xl[i];
      }
      if (xopt[i] == su[i]) {
        x[i] = xu[i];
      }
    }
    f = fval[kopt];
  }
  if (iprint >= 1) {
    fprintf(OUTPUT, "\n"
            "    At the return from BOBYQA "
            "    Number of function values =%6ld\n"
            "    Least value of F =%23.15E     "
            "    The corresponding X is:\n",
            (long)nf, (double)f);
    print_x(OUTPUT, n, &x[1], NULL);
  }

  /* Save function value at the solution and return status. */
  if (status == BOBYQA_SUCCESS) {
    xbase[1] = f;
  }
  return status;

 too_many_evaluations:
  reason = "CALFUN has been called MAXFUN times";
  status = BOBYQA_TOO_MANY_EVALUATIONS;
  goto error;

 cancellation_of_denominator:
  reason = "of much cancellation in a denominator";
  status = BOBYQA_ROUNDING_ERRORS;
  goto error;

 step_failed:
  reason = "a trust region step has failed to reduce Q";
  status = BOBYQA_STEP_FAILED;
  goto error;

 error:
  if (iprint > 0) {
    print_error(reason);
  }
  goto done;

} /* bobyqb */

#undef ZMAT
#undef BMAT
#undef XPT

static void
altmov(const INTEGER n, const INTEGER npt, REAL xpt[],
       REAL* xopt, REAL* bmat, REAL* zmat, const INTEGER ndim,
       REAL* sl, REAL* su, const INTEGER kopt, const INTEGER knew,
       const REAL adelt, REAL* xnew, REAL* xalt, REAL* alpha,
       REAL* cauchy, REAL* glag, REAL* hcol, REAL* w)
{
  /* The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have the
     same meanings as the corresponding arguments of BOBYQB.

     KOPT is the index of the optimal interpolation point.

     KNEW is the index of the interpolation point that is going to be moved.

     ADELT is the current trust region bound.

     XNEW will be set to a suitable new position for the interpolation point
     XPT(KNEW,.).  Specifically, it satisfies the SL, SU and trust region
     bounds and it should provide a large denominator in the next call of
     UPDATE.  The step XNEW-XOPT from XOPT is restricted to moves along the
     straight lines through XOPT and another interpolation point.

     XALT also provides a large value of the modulus of the KNEW-th Lagrange
     function subject to the constraints that have been mentioned, its main
     difference from XNEW being that XALT-XOPT is a constrained version of the
     Cauchy step within the trust region.  An exception is that XALT is not
     calculated if all components of GLAG (see below) are zero.

     ALPHA will be set to the KNEW-th diagonal element of the H matrix.

     CAUCHY will be set to the square of the KNEW-th Lagrange function at the
     step XALT-XOPT from XOPT for the vector XALT that is returned, except that
     CAUCHY is set to zero if XALT is not calculated.

     GLAG is a working space vector of length N for the gradient of the KNEW-th
     Lagrange function at XOPT.

     HCOL is a working space vector of length NPT for the second derivative
     coefficients of the KNEW-th Lagrange function.

     W is a working space vector of length 2N that is going to hold the
     constrained Cauchy step from XOPT of the Lagrange function, followed by
     the downhill version of XALT when the uphill step is calculated. */

  /* Constants. */
  const REAL half = 0.5;
  const REAL one = 1.0;
  const REAL zero = 0.0;
  const REAL one_plus_sqrt2 = one + M_SQRT2;

  /* Local variables. */
  REAL bigstp, csave, curv, dderiv, diff, distsq, ggfree, gw, ha,
    predsq, presav, scale, slbd, step, stpsav, subd, sumin,
    temp, tempa, tempb, tempd, vlag, wfixsq, wsqsav;
  INTEGER i, i1, ibdsav, iflag, ilbd, isbd, iubd, j, k, ksav;

  /* Parameter adjustments to comply with FORTRAN indexing. */
  xpt  -= 1 + npt;
  xopt -= 1;
  bmat -= 1 + ndim;
  zmat -= 1 + npt;
  sl   -= 1;
  su   -= 1;
  xnew -= 1;
  xalt -= 1;
  glag -= 1;
  hcol -= 1;
  w    -= 1;
#define XPT(a1,a2) xpt[(a2)*npt + a1]
#define BMAT(a1,a2) bmat[(a2)*ndim + a1]
#define ZMAT(a1,a2) zmat[(a2)*npt + a1]

  /* FIXME: Set uninitialized variables. */
  csave = zero;
  stpsav = zero;
  step = zero;
  ksav = 0;
  ibdsav = 0;

  /* Set the first NPT components of W to the leading elements of the KNEW-th
     column of the H matrix. */
  LOOP(k,npt) {
    hcol[k] = zero;
  }
  i1 = npt - n - 1;
  LOOP(j,i1) {
    temp = ZMAT(knew, j);
    LOOP(k,npt) {
      hcol[k] += temp*ZMAT(k,j);
    }
  }
  *alpha = hcol[knew];
  ha = half*(*alpha);

  /* Calculate the gradient of the KNEW-th Lagrange function at XOPT. */
  LOOP(i,n) {
    glag[i] = BMAT(knew, i);
  }
  LOOP(k,npt) {
    temp = zero;
    LOOP(j,n) {
      temp += XPT(k,j)*xopt[j];
    }
    temp = hcol[k]*temp;
    LOOP(i,n) {
      glag[i] += temp*XPT(k,i);
    }
  }

  /* Search for a large denominator along the straight lines through XOPT and
     another interpolation point.  SLBD and SUBD will be lower and upper bounds
     on the step along each of these lines in turn.  PREDSQ will be set to the
     square of the predicted denominator for each line.  PRESAV will be set to
     the largest admissible value of PREDSQ that occurs. */
  presav = zero;
  LOOP(k,npt) {
    if (k == kopt) {
      continue;
    }
    dderiv = zero;
    distsq = zero;
    LOOP(i,n) {
      temp = XPT(k,i) - xopt[i];
      dderiv += glag[i]*temp;
      distsq += temp*temp;
    }
    subd = adelt/SQRT(distsq);
    slbd = -subd;
    ilbd = 0;
    iubd = 0;
    sumin = MIN(one,subd);

    /* Revise SLBD and SUBD if necessary because of the bounds in SL and SU. */
    LOOP(i,n) {
      temp = XPT(k,i) - xopt[i];
      if (temp > zero) {
        if (slbd*temp < sl[i] - xopt[i]) {
          slbd = (sl[i] - xopt[i])/temp;
          ilbd = -i;
        }
        if (subd*temp > su[i] - xopt[i]) {
          subd = (su[i] - xopt[i])/temp;
          subd = MAX(subd,sumin);
          iubd = i;
        }
      } else if (temp < zero) {
        if (slbd*temp > su[i] - xopt[i]) {
          slbd = (su[i] - xopt[i])/temp;
          ilbd = i;
        }
        if (subd*temp < sl[i] - xopt[i]) {
          subd = (sl[i] - xopt[i])/temp;
          subd = MAX(subd,sumin);
          iubd = -i;
        }
      }
    }

    /* Seek a large modulus of the KNEW-th Lagrange function when the index of
       the other interpolation point on the line through XOPT is KNEW. */
    if (k == knew) {
      diff = dderiv - one;
      step = slbd;
      vlag = slbd*(dderiv - slbd*diff);
      isbd = ilbd;
      temp = subd*(dderiv - subd*diff);
      if (ABS(temp) > ABS(vlag)) {
        step = subd;
        vlag = temp;
        isbd = iubd;
      }
      tempd = half*dderiv;
      tempa = tempd - diff*slbd;
      tempb = tempd - diff*subd;
      if (tempa*tempb < zero) {
        temp = tempd*tempd/diff;
        if (ABS(temp) > ABS(vlag)) {
          step = tempd/diff;
          vlag = temp;
          isbd = 0;
        }
      }

      /* Search along each of the other lines through XOPT and another
         point. */
    } else {
      step = slbd;
      vlag = slbd*(one - slbd);
      isbd = ilbd;
      temp = subd*(one - subd);
      if (ABS(temp) > ABS(vlag)) {
        step = subd;
        vlag = temp;
        isbd = iubd;
      }
      if (subd > half) {
        if (ABS(vlag) < 0.25) {
          step = half;
          vlag = 0.25;
          isbd = 0;
        }
      }
      vlag *= dderiv;
    }

    /* Calculate PREDSQ for the current line search and maintain PRESAV. */
    temp = step*(one - step)*distsq;
    predsq = vlag*vlag*(vlag*vlag + ha*temp*temp);
    if (predsq > presav) {
      presav = predsq;
      ksav = k;
      stpsav = step;
      ibdsav = isbd;
    }
  }

  /* Construct XNEW in a way that satisfies the bound constraints exactly. */
  LOOP(i,n) {
    temp = xopt[i] + stpsav*(XPT(ksav,i) - xopt[i]);
    temp = MIN(temp,su[i]);
    xnew[i] = MAX(temp,sl[i]);
  }
  if (ibdsav < 0) {
    xnew[-ibdsav] = sl[-ibdsav];
  }
  if (ibdsav > 0) {
    xnew[ibdsav] = su[ibdsav];
  }

  /* Prepare for the iterative method that assembles the constrained Cauchy
     step in W.  The sum of squares of the fixed components of W is formed in
     WFIXSQ, and the free components of W are set to BIGSTP. */
  bigstp = adelt + adelt;
  iflag = 0;
  for (;;) {
    wfixsq = zero;
    ggfree = zero;
    LOOP(i,n) {
      w[i] = zero;
      tempa = xopt[i] - sl[i];
      tempa = MIN(tempa,glag[i]);
      tempb = xopt[i] - su[i];
      tempb = MAX(tempb,glag[i]);
      if (tempa > zero || tempb < zero) {
        w[i] = bigstp;
        ggfree += glag[i]*glag[i];
      }
    }
    if (ggfree == zero) {
      *cauchy = zero;
      return;
    }

    /* Investigate whether more components of W can be fixed. */
    do {
      temp = adelt*adelt - wfixsq;
      if (temp <= zero) break;
      wsqsav = wfixsq;
      step = SQRT(temp/ggfree);
      ggfree = zero;
      LOOP(i,n) {
        if (w[i] == bigstp) {
          temp = xopt[i] - step*glag[i];
          if (temp <= sl[i]) {
            w[i] = sl[i] - xopt[i];
            wfixsq += w[i]*w[i];
          } else if (temp >= su[i]) {
            w[i] = su[i] - xopt[i];
            wfixsq += w[i]*w[i];
          } else {
            ggfree += glag[i]*glag[i];
          }
        }
      }
    } while (wfixsq > wsqsav && ggfree > zero);

    /* Set the remaining free components of W and all components of XALT,
       except that W may be scaled later. */
    gw = zero;
    LOOP(i,n) {
      if (w[i] == bigstp) {
        w[i] = -step*glag[i];
        temp = xopt[i] + w[i];
        temp = MIN(temp,su[i]);
        xalt[i] = MAX(temp,sl[i]);
      } else if (w[i] == zero) {
        xalt[i] = xopt[i];
      } else if (glag[i] > zero) {
        xalt[i] = sl[i];
      } else {
        xalt[i] = su[i];
      }
      gw += glag[i]*w[i];
    }

    /* Set CURV to the curvature of the KNEW-th Lagrange function along W.
       Scale W by a factor less than one if that can reduce the modulus of the
       Lagrange function at XOPT+W.  Set CAUCHY to the final value of the
       square of this function. */
    curv = zero;
    LOOP(k,npt) {
      temp = zero;
      LOOP(j,n) {
        temp += XPT(k,j)*w[j];
      }
      curv += hcol[k]*temp*temp;
    }
    if (iflag == 1) {
      curv = -curv;
    }
    if (curv > -gw && curv < -one_plus_sqrt2*gw) {
      scale = -gw/curv;
      LOOP(i,n) {
        temp = xopt[i] + scale*w[i];
        temp = MIN(temp,su[i]);
        xalt[i] = MAX(temp,sl[i]);
      }
      temp = half*gw*scale;
      *cauchy = temp*temp;
    } else {
      temp = gw + half*curv;
      *cauchy = temp*temp;
    }

    /* If IFLAG is zero, then XALT is calculated as before after reversing the
       sign of GLAG.  Thus two XALT vectors become available.  The one that is
       chosen is the one that gives the larger value of CAUCHY. */
    if (iflag != 0) {
      break;
    }
    LOOP(i,n) {
      glag[i] = -glag[i];
      w[n + i] = xalt[i];
    }
    csave = *cauchy;
    iflag = 1;
  }
  if (csave > *cauchy) {
    LOOP(i,n) {
      xalt[i] = w[n + i];
    }
    *cauchy = csave;
  }
} /* altmov */

#undef ZMAT
#undef BMAT
#undef XPT

static void
prelim(const INTEGER n, const INTEGER npt,
       bobyqa_objfun* objfun, void* data,
       REAL* x, const REAL* xl, const REAL* xu,
       const REAL rhobeg, const INTEGER iprint,
       const INTEGER maxfun, REAL* xbase, REAL* xpt, REAL* fval,
       REAL* gopt, REAL* hq, REAL* pq, REAL* bmat,
       REAL* zmat, const INTEGER ndim, REAL* sl, REAL* su,
       INTEGER* nf, INTEGER* kopt)
{
  /* The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the same as
     the corresponding arguments in SUBROUTINE BOBYQA.  The arguments XBASE,
     XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU are the same as the
     corresponding arguments in BOBYQB, the elements of SL and SU being set in
     BOBYQA.

     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but it
     is set by PRELIM to the gradient of the quadratic model at XBASE.

     If XOPT is nonzero, BOBYQB will change it to its usual value later.

     NF is maintaned as the number of calls of OBJFUN so far.

     KOPT will be such that the least calculated value of F so far is at the
     point XPT(KOPT,.)+XBASE in the space of the variables.

     SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
     BMAT and ZMAT for the first iteration, and it maintains the values of NF
     and KOPT.  The vector X is also changed by PRELIM. */

  /* Constants. */
  const REAL half = 0.5;
  const REAL one = 1.0;
  const REAL two = 2.0;
  const REAL zero = 0.0;

  /* Local variables. */
  REAL diff, f, fbeg, recip, rhosq, stepa, stepb, temp;
  INTEGER i, i1, ih, ipt, itemp, j, jpt, k, nfm, nfx, np;

  /* Parameter adjustments to comply with FORTRAN indexing. */
  x     -= 1;
  xl    -= 1;
  xu    -= 1;
  xbase -= 1;
  xpt   -= 1 + npt;
  fval  -= 1;
  gopt  -= 1;
  hq    -= 1;
  pq    -= 1;
  bmat  -= 1 + ndim;
  zmat  -= 1 + npt;
  sl    -= 1;
  su    -= 1;
#define XPT(a1,a2) xpt[(a2)*npt + a1]
#define BMAT(a1,a2) bmat[(a2)*ndim + a1]
#define ZMAT(a1,a2) zmat[(a2)*npt + a1]

  /* FIXME: Set uninitialized variables. */
  stepa = zero;
  stepb = zero;
  fbeg = zero;
  ipt = 0;
  jpt = 0;

  /* Set some constants. */
  rhosq = rhobeg*rhobeg;
  recip = one/rhosq;
  np = n + 1;

  /* Set XBASE to the initial vector of variables, and set the initial elements
     of XPT, BMAT, HQ, PQ and ZMAT to zero. */
  LOOP(j,n) {
    xbase[j] = x[j];
    LOOP(k,npt) {
      XPT(k,j) = zero;
    }
    LOOP(i,ndim) {
      BMAT(i,j) = zero;
    }
  }
  i1 = n*np/2;
  LOOP(ih,i1) {
    hq[ih] = zero;
  }
  LOOP(k,npt) {
    pq[k] = zero;
    i1 = npt - np;
    LOOP(j,i1) {
      ZMAT(k,j) = zero;
    }
  }

  /* Begin the initialization procedure.  NF becomes one more than the number
     of function values so far.  The coordinates of the displacement of the
     next initial interpolation point from XBASE are set in XPT(NF+1,.). */
  *nf = 0;
  do {
    nfm = *nf;
    nfx = *nf - n;
    ++(*nf);
    if (nfm <= 2*n) {
      if (nfm >= 1 && nfm <= n) {
        stepa = rhobeg;
        if (su[nfm] == zero) {
          stepa = -stepa;
        }
        XPT(*nf, nfm) = stepa;
      } else if (nfm > n) {
        stepa = XPT(*nf - n, nfx);
        stepb = -rhobeg;
        if (sl[nfx] == zero) {
          stepb = two*rhobeg;
          stepb = MIN(stepb,su[nfx]);
        }
        if (su[nfx] == zero) {
          stepb = -two*rhobeg;
          stepb = MAX(stepb,sl[nfx]);
        }
        XPT(*nf, nfx) = stepb;
      }
    } else {
      itemp = (nfm - np)/n;
      jpt = nfm - itemp*n - n;
      ipt = jpt + itemp;
      if (ipt > n) {
        itemp = jpt;
        jpt = ipt - n;
        ipt = itemp;
      }
      XPT(*nf, ipt) = XPT(ipt + 1, ipt);
      XPT(*nf, jpt) = XPT(jpt + 1, jpt);
    }

    /* Calculate the next value of F.  The least function value so far and its
       index are required. */
    LOOP(j,n) {
      temp = xbase[j] + XPT(*nf, j);
      temp = MAX(temp,xl[j]);
      x[j] = MIN(temp,xu[j]);
      if (XPT(*nf, j) == sl[j]) {
        x[j] = xl[j];
      }
      if (XPT(*nf, j) == su[j]) {
        x[j] = xu[j];
      }
    }
    f = objfun(n, &x[1], data);
    if (iprint == 3) {
      fprintf(OUTPUT, "Function number%6ld    F = %.18G"
              "    The corresponding X is: ", (long)*nf, (double)f);
      LOOP(i,n) {
        fprintf(OUTPUT, " %15.6E", x[i]);
      }
      fprintf(OUTPUT, "\n");
    }
    fval[*nf] = f;
    if (*nf == 1) {
      fbeg = f;
      *kopt = 1;
    } else if (f < fval[*kopt]) {
      *kopt = *nf;
    }

    /* Set the nonzero initial elements of BMAT and the quadratic model in the
       cases when NF is at most 2*N+1.  If NF exceeds N+1, then the positions
       of the NF-th and (NF-N)-th interpolation points may be switched, in
       order that the function value at the first of them contributes to the
       off-diagonal second derivative terms of the initial quadratic model. */
    if (*nf <= 2*n + 1) {
      if (*nf >= 2 && *nf <= n + 1) {
        gopt[nfm] = (f - fbeg)/stepa;
        if (npt < *nf + n) {
          BMAT(1,nfm) = -one/stepa;
          BMAT(*nf, nfm) = one/stepa;
          BMAT(npt + nfm, nfm) = -half*rhosq;
        }
      } else if (*nf >= n + 2) {
        ih = nfx*(nfx + 1)/2;
        temp = (f - fbeg)/stepb;
        diff = stepb - stepa;
        hq[ih] = two*(temp - gopt[nfx])/diff;
        gopt[nfx] = (gopt[nfx]*stepb - temp*stepa)/diff;
        if (stepa*stepb < zero) {
          if (f < fval[*nf - n]) {
            fval[*nf] = fval[*nf - n];
            fval[*nf - n] = f;
            if (*kopt == *nf) {
              *kopt = *nf - n;
            }
            XPT(*nf - n, nfx) = stepb;
            XPT(*nf, nfx) = stepa;
          }
        }
        BMAT(1,nfx) = -(stepa + stepb)/(stepa*stepb);
        BMAT(*nf, nfx) = -half/XPT(*nf - n, nfx);
        BMAT(*nf - n, nfx) = -BMAT(1,nfx) - BMAT(*nf, nfx);
        ZMAT(1,nfx) = SQRT(two)/(stepa*stepb);
        ZMAT(*nf, nfx) = SQRT(half)/rhosq;
        ZMAT(*nf - n, nfx) = -ZMAT(1,nfx) - ZMAT(*nf, nfx);
      }
    } else {
      /* Set the off-diagonal second derivatives of the Lagrange functions and
         the initial quadratic model. */
      ih = ipt*(ipt - 1)/2 + jpt;
      ZMAT(1,nfx) = recip;
      ZMAT(*nf, nfx) = recip;
      ZMAT(ipt + 1, nfx) = -recip;
      ZMAT(jpt + 1, nfx) = -recip;
      temp = XPT(*nf, ipt)*XPT(*nf, jpt);
      hq[ih] = (fbeg - fval[ipt + 1] - fval[jpt + 1] + f)/temp;
    }
  } while (*nf < npt && *nf < maxfun);
} /* prelim */

#undef ZMAT
#undef BMAT
#undef XPT

static void
rescue(const INTEGER n, const INTEGER npt,
       bobyqa_objfun* objfun, void* data,
       const REAL* xl, const REAL* xu,
       const INTEGER iprint, const INTEGER maxfun,
       REAL* xbase, REAL* xpt, REAL* fval, REAL* xopt,
       REAL* gopt, REAL* hq, REAL* pq, REAL* bmat, REAL* zmat,
       const INTEGER ndim, REAL* sl, REAL* su, INTEGER* nf,
       const REAL delta, INTEGER* kopt, REAL* vlag, REAL* ptsaux,
       REAL* ptsid, REAL* w)
{
  /* The arguments N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT,
     GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU have the same meanings as the
     corresponding arguments of BOBYQB on the entry to RESCUE.

     NF is maintained as the number of calls of OBJFUN so far, except that NF
     is set to -1 if the value of MAXFUN prevents further progress.

     KOPT is maintained so that FVAL(KOPT) is the least calculated function
     value.  Its correct value must be given on entry.  It is updated if a new
     least function value is found, but the corresponding changes to XOPT and
     GOPT have to be made later by the calling program.

     DELTA is the current trust region radius.

     VLAG is a working space vector that will be used for the values of the
     provisional Lagrange functions at each of the interpolation points.  They
     are part of a product that requires VLAG to be of length NDIM.

     PTSAUX is also a working space array.  For J=1,2,...,N, PTSAUX(1,J) and
     PTSAUX(2,J) specify the two positions of provisional interpolation points
     when a nonzero step is taken along e_J (the J-th coordinate direction)
     through XBASE+XOPT, as specified below.  Usually these steps have length
     DELTA, but other lengths are chosen if necessary in order to satisfy the
     given bounds on the variables.

     PTSID is also a working space array.  It has NPT components that denote
     provisional new positions of the original interpolation points, in case
     changes are needed to restore the linear independence of the interpolation
     conditions.  The K-th point is a candidate for change if and only if
     PTSID(K) is nonzero.  In this case let p and q be the integer parts of
     PTSID(K) and (PTSID(K)-p) multiplied by N+1.  If p and q are both positive,
     the step from XBASE+XOPT to the new K-th interpolation point is
     PTSAUX(1,p)*e_p + PTSAUX(1,q)*e_q.  Otherwise the step is PTSAUX(1,p)*e_p
     or PTSAUX(2,q)*e_q in the cases q=0 or p=0, respectively.

     The first NDIM+NPT elements of the array W are used for working space.
     The final elements of BMAT and ZMAT are set in a well-conditioned way to
     the values that are appropriate for the new interpolation points.  The
     elements of GOPT, HQ and PQ are also revised to the values that are
     appropriate to the final quadratic model. */

  /* Constants. */
  const REAL half = 0.5;
  const REAL one = 1.0;
  const REAL zero = 0.0;

  /* Local variables. */
  REAL beta, bsum, den, denom, diff, distsq, dsqmin, f, fbase, hdiag,
    sfrac, sum, sumpq, temp, vlmxsq, vquad, winc, xp, xq;
  INTEGER i, ih, ihp, ihq, ip, iq, iw, j, jp, jpn, k, knew, kold, kpt,
    np, nptm, nrem;

  /* Parameter adjustments to comply with FORTRAN indexing. */
  zmat   -= 1 + npt;
  xpt    -= 1 + npt;
  xl     -= 1;
  xu     -= 1;
  xbase  -= 1;
  fval   -= 1;
  xopt   -= 1;
  gopt   -= 1;
  hq     -= 1;
  pq     -= 1;
  bmat   -= 1 + ndim;
  sl     -= 1;
  su     -= 1;
  vlag   -= 1;
  ptsaux -= 3;
  ptsid  -= 1;
  w      -= 1;
#define XPT(a1,a2) xpt[(a2)*npt + a1]
#define BMAT(a1,a2) bmat[(a2)*ndim + a1]
#define ZMAT(a1,a2) zmat[(a2)*npt + a1]
#define PTSAUX(a1,a2) ptsaux[(a2)*2 + a1]

  /* FIXME: Set uninitialized variables. */
  beta = zero;
  denom = zero;

  /* Set some constants. */
  np = n + 1;
  sfrac = half/(REAL)np;
  nptm = npt - np;

  /* Shift the interpolation points so that XOPT becomes the origin, and set
     the elements of ZMAT to zero.  The value of SUMPQ is required in the
     updating of HQ below.  The squares of the distances from XOPT to the other
     interpolation points are set at the end of W.  Increments of WINC may be
     added later to these squares to balance the consideration of the choice of
     point that is going to become current. */
  sumpq = zero;
  winc = zero;
  LOOP(k,npt) {
    distsq = zero;
    LOOP(j,n) {
      XPT(k,j) = XPT(k,j) - xopt[j];
      distsq += XPT(k,j)*XPT(k,j);
    }
    sumpq += pq[k];
    w[ndim + k] = distsq;
    winc = MAX(winc,distsq);
    LOOP(j,nptm) {
      ZMAT(k,j) = zero;
    }
  }

  /* Update HQ so that HQ and PQ define the second derivatives of the model
     after XBASE has been shifted to the trust region centre. */
  ih = 0;
  LOOP(j,n) {
    w[j] = half*sumpq*xopt[j];
    LOOP(k,npt) {
      w[j] += pq[k]*XPT(k,j);
    }
    LOOP(i,j) {
      ++ih;
      hq[ih] = hq[ih] + w[i]*xopt[j] + w[j]*xopt[i];
    }
  }

  /* Shift XBASE, SL, SU and XOPT.  Set the elements of BMAT to zero, and also
     set the elements of PTSAUX. */
  LOOP(j,n) {
    xbase[j] += xopt[j];
    sl[j] -= xopt[j];
    su[j] -= xopt[j];
    xopt[j] = zero;
    PTSAUX(1,j) = MIN( delta, su[j]);
    PTSAUX(2,j) = MAX(-delta, sl[j]);
    if (PTSAUX(1,j) + PTSAUX(2,j) < zero) {
      temp = PTSAUX(1,j);
      PTSAUX(1,j) = PTSAUX(2,j);
      PTSAUX(2,j) = temp;
    }
    if (ABS(PTSAUX(2,j)) < half*ABS(PTSAUX(1,j))) {
      PTSAUX(2,j) = half*PTSAUX(1,j);
    }
    LOOP(i,ndim) {
      BMAT(i,j) = zero;
    }
  }
  fbase = fval[*kopt];

  /* Set the identifiers of the artificial interpolation points that are along
     a coordinate direction from XOPT, and set the corresponding nonzero
     elements of BMAT and ZMAT. */
  ptsid[1] = sfrac;
  LOOP(j,n) {
    jp = j + 1;
    jpn = jp + n;
    ptsid[jp] = (REAL)j + sfrac;
    if (jpn <= npt) {
      ptsid[jpn] = (REAL)j/(REAL)np + sfrac;
      temp = one/(PTSAUX(1,j) - PTSAUX(2,j));
      BMAT(jp,j) = -temp + one/PTSAUX(1,j);
      BMAT(jpn,j) = temp + one/PTSAUX(2,j);
      BMAT(1,j) = -BMAT(jp,j) - BMAT(jpn,j);
      ZMAT(1,j) = SQRT(2.0)/ABS(PTSAUX(1,j)*PTSAUX(2,j));
      ZMAT(jp,j) = ZMAT(1,j)*PTSAUX(2,j)*temp;
      ZMAT(jpn,j) = -ZMAT(1,j)*PTSAUX(1,j)*temp;
    } else {
      BMAT(1,j) = -one/PTSAUX(1,j);
      BMAT(jp,j) = one/PTSAUX(1,j);
      BMAT(j + npt, j) = -half*(PTSAUX(1,j)*PTSAUX(1,j));
    }
  }

  /* Set any remaining identifiers with their nonzero elements of ZMAT. */
  if (npt >= n + np) {
    for (k = 2*np; k <= npt; ++k) {
      iw = (INTEGER) (((REAL)(k - np) - half)/(REAL)n);
      ip = k - np - iw*n;
      iq = ip + iw;
      if (iq > n) {
        iq -= n;
      }
      ptsid[k] = (REAL)ip + (REAL)iq/(REAL)np + sfrac;
      temp = one/(PTSAUX(1,ip)*PTSAUX(1,iq));
      ZMAT(1, k - np) = temp;
      ZMAT(ip + 1, k - np) = -temp;
      ZMAT(iq + 1, k - np) = -temp;
      ZMAT(k, k - np) = temp;
    }
  }
  nrem = npt;
  kold = 1;
  knew = *kopt;

  /* Reorder the provisional points in the way that exchanges PTSID(KOLD) with
     PTSID(KNEW). */
 L80:
  LOOP(j,n) {
    temp = BMAT(kold,j);
    BMAT(kold,j) = BMAT(knew,j);
    BMAT(knew,j) = temp;
  }
  LOOP(j,nptm) {
    temp = ZMAT(kold,j);
    ZMAT(kold,j) = ZMAT(knew,j);
    ZMAT(knew,j) = temp;
  }
  ptsid[kold] = ptsid[knew];
  ptsid[knew] = zero;
  w[ndim + knew] = zero;
  --nrem;
  if (knew != *kopt) {
    temp = vlag[kold];
    vlag[kold] = vlag[knew];
    vlag[knew] = temp;

    /* Update the BMAT and ZMAT matrices so that the status of the KNEW-th
       interpolation point can be changed from provisional to original.  The
       branch to label 350 occurs if all the original points are reinstated.
       The nonnegative values of W(NDIM+K) are required in the search below. */
    update(n, npt, &BMAT(1,1), &ZMAT(1,1), ndim, &vlag[1],
           beta, denom, knew, &w[1]);
    if (nrem == 0) {
      goto L350;
    }
    LOOP(k,npt) {
      w[ndim + k] = ABS(w[ndim + k]);
    }
  }

  /* Pick the index KNEW of an original interpolation point that has not yet
     replaced one of the provisional interpolation points, giving attention to
     the closeness to XOPT and to previous tries with KNEW. */
 L120:
  dsqmin = zero;
  LOOP(k,npt) {
    if (w[ndim + k] > zero) {
      if (dsqmin == zero || w[ndim + k] < dsqmin) {
        knew = k;
        dsqmin = w[ndim + k];
      }
    }
  }
  if (dsqmin == zero) {
    goto L260;
  }

  /* Form the W-vector of the chosen original interpolation point. */
  LOOP(j,n) {
    w[npt + j] = XPT(knew,j);
  }
  LOOP(k,npt) {
    sum = zero;
    if (k == *kopt) {
    } else if (ptsid[k] == zero) {
      LOOP(j,n) {
        sum += w[npt + j]*XPT(k,j);
      }
    } else {
      ip = (INTEGER) ptsid[k];
      if (ip > 0) {
        sum = w[npt + ip]*PTSAUX(1,ip);
      }
      iq = (INTEGER) ((REAL)np*ptsid[k] - (REAL)(ip*np));
      if (iq > 0) {
        iw = 1;
        if (ip == 0) {
          iw = 2;
        }
        sum += w[npt + iq]*PTSAUX(iw,iq);
      }
    }
    w[k] = half*sum*sum;
  }

  /* Calculate VLAG and BETA for the required updating of the H matrix if
     XPT(KNEW,.) is reinstated in the set of interpolation points. */
  LOOP(k,npt) {
    sum = zero;
    LOOP(j,n) {
      sum += BMAT(k,j)*w[npt + j];
    }
    vlag[k] = sum;
  }
  beta = zero;
  LOOP(j,nptm) {
    sum = zero;
    LOOP(k,npt) {
      sum += ZMAT(k,j)*w[k];
    }
    beta -= sum*sum;
    LOOP(k,npt) {
      vlag[k] += sum*ZMAT(k,j);
    }
  }
  bsum = zero;
  distsq = zero;
  LOOP(j,n) {
    sum = zero;
    LOOP(k,npt) {
      sum += BMAT(k,j)*w[k];
    }
    jp = j + npt;
    bsum += sum*w[jp];
    for (ip = npt + 1; ip <= ndim; ++ip) {
      sum += BMAT(ip,j)*w[ip];
    }
    bsum += sum*w[jp];
    vlag[jp] = sum;
    distsq += XPT(knew,j)*XPT(knew,j);
  }
  beta = half*distsq*distsq + beta - bsum;
  vlag[*kopt] += one;

  /* KOLD is set to the index of the provisional interpolation point that is
     going to be deleted to make way for the KNEW-th original interpolation
     point.  The choice of KOLD is governed by the avoidance of a small value
     of the denominator in the updating calculation of UPDATE. */
  denom = zero;
  vlmxsq = zero;
  LOOP(k,npt) {
    if (ptsid[k] != zero) {
      hdiag = zero;
      LOOP(j,nptm) {
        hdiag += ZMAT(k,j)*ZMAT(k,j);
      }
      den = beta*hdiag + vlag[k]*vlag[k];
      if (den > denom) {
        kold = k;
        denom = den;
      }
    }
    temp = vlag[k]*vlag[k];
    vlmxsq = MAX(vlmxsq,temp);
  }
  if (denom <= vlmxsq*0.01) {
    w[ndim + knew] = -w[ndim + knew] - winc;
    goto L120;
  }
  goto L80;

  /* When label 260 is reached, all the final positions of the interpolation
     points have been chosen although any changes have not been included yet in
     XPT.  Also the final BMAT and ZMAT matrices are complete, but, apart from
     the shift of XBASE, the updating of the quadratic model remains to be
     done.  The following cycle through the new interpolation points begins by
     putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero, except
     that a RETURN occurs if MAXFUN prohibits another value of F. */
 L260:
  LOOP(kpt,npt) {
    if (ptsid[kpt] == zero) {
      continue;
    }
    if (*nf >= maxfun) {
      *nf = -1;
      goto L350;
    }
    ih = 0;
    LOOP(j,n) {
      w[j] = XPT(kpt,j);
      XPT(kpt,j) = zero;
      temp = pq[kpt]*w[j];
      LOOP(i,j) {
        ++ih;
        hq[ih] += temp*w[i];
      }
    }
    pq[kpt] = zero;
    ip = (INTEGER) ptsid[kpt];
    iq = (INTEGER) ((REAL)np*ptsid[kpt] - (REAL)(ip*np))
      ;
    if (ip > 0) {
      xp = PTSAUX(1,ip);
      XPT(kpt,ip) = xp;
    }
    if (iq > 0) {
      xq = PTSAUX(1,iq);
      if (ip == 0) {
        xq = PTSAUX(2,iq);
      }
      XPT(kpt,iq) = xq;
    }

    /* Set VQUAD to the value of the current model at the new point. */
    vquad = fbase;
    if (ip > 0) {
      ihp = (ip + ip*ip)/2;
      vquad += xp*(gopt[ip] + half*xp*hq[ihp]);
    }
    if (iq > 0) {
      ihq = (iq + iq*iq)/2;
      vquad += xq*(gopt[iq] + half*xq*hq[ihq]);
      if (ip > 0) {
        iw = MAX(ihp,ihq) - (ip >= iq ? ip - iq : iq - ip);
        vquad += xp*xq*hq[iw];
      }
    }
    LOOP(k,npt) {
      temp = zero;
      if (ip > 0) {
        temp += xp*XPT(k,ip);
      }
      if (iq > 0) {
        temp += xq*XPT(k,iq);
      }
      vquad += half*pq[k]*temp*temp;
    }

    /* Calculate F at the new interpolation point, and set DIFF to the factor
       that is going to multiply the KPT-th Lagrange function when the model is
       updated to provide interpolation to the new function value. */
    LOOP(i,n) {
      temp = xbase[i] + XPT(kpt,i);
      temp = MAX(temp,xl[i]);
      w[i] = MIN(temp,xu[i]);
      if (XPT(kpt,i) == sl[i]) {
        w[i] = xl[i];
      }
      if (XPT(kpt,i) == su[i]) {
        w[i] = xu[i];
      }
    }
    ++(*nf);
    f = objfun(n, &w[1], data);
    if (iprint == 3) {
      fprintf(OUTPUT,
              "    Function number%6ld"
              "    F =%18.10E"
              "    The corresponding X is:\n",
              (long)*nf, (double)f);
      print_x(OUTPUT, n, &w[1], NULL);
    }
    fval[kpt] = f;
    if (f < fval[*kopt]) {
      *kopt = kpt;
    }
    diff = f - vquad;

    /* Update the quadratic model.  The RETURN from the subroutine occurs when
       all the new interpolation points are included in the model. */
    LOOP(i,n) {
      gopt[i] += diff*BMAT(kpt,i);
    }
    LOOP(k,npt) {
      sum = zero;
      LOOP(j,nptm) {
        sum += ZMAT(k,j)*ZMAT(kpt,j);
      }
      temp = diff*sum;
      if (ptsid[k] == zero) {
        pq[k] += temp;
      } else {
        ip = (INTEGER) ptsid[k];
        iq = (INTEGER) ((REAL)np*ptsid[k] - (REAL)(ip*np));
        ihq = (iq*iq + iq)/2;
        if (ip == 0) {
          hq[ihq] += temp*(PTSAUX(2,iq)*PTSAUX(2,iq));
        } else {
          ihp = (ip*ip + ip)/2;
          hq[ihp] += temp*(PTSAUX(1,ip)*PTSAUX(1,ip));
          if (iq > 0) {
            hq[ihq] += temp*(PTSAUX(1,iq)*PTSAUX(1,iq));
            iw = MAX(ihp,ihq) - (ip >= iq ? ip - iq : iq - ip);
            hq[iw] += temp*PTSAUX(1,ip)*PTSAUX(1,iq);
          }
        }
      }
    }
    ptsid[kpt] = zero;
  }
 L350:
  return;
} /* rescue */

#undef PTSAUX
#undef ZMAT
#undef BMAT
#undef XPT

static void
trsbox(const INTEGER n, const INTEGER npt, REAL* xpt,
       REAL* xopt, REAL* gopt, REAL* hq, REAL* pq,
       REAL* sl, REAL* su, const REAL delta, REAL* xnew,
       REAL* d, REAL* gnew, REAL* xbdi, REAL* s,
       REAL* hs, REAL* hred, REAL* dsq, REAL* crvmin)
{
  /* The arguments N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same
     meanings as the corresponding arguments of BOBYQB.

     DELTA is the trust region radius for the present calculation, which seeks
     a small value of the quadratic model within distance DELTA of

     XOPT subject to the bounds on the variables.

     XNEW will be set to a new vector of variables that is approximately the
     one that minimizes the quadratic model within the trust region subject to
     the SL and SU constraints on the variables.  It satisfies as equations the
     bounds that become active during the calculation.

     D is the calculated trial step from XOPT, generated iteratively from an
     initial value of zero.  Thus XNEW is XOPT+D after the final iteration.

     GNEW holds the gradient of the quadratic model at XOPT+D.  It is updated
     when D is updated.

     XBDI is a working space vector.  For I=1,2,...,N, the element XBDI(I) is
     set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the I-th
     variable has become fixed at a bound, the bound being SL(I) or SU(I) in
     the case XBDI(I)=-1.0 or XBDI(I)=1.0, respectively.  This information is
     accumulated during the construction of XNEW.  The arrays S, HS and HRED
     are also used for working space.  They hold the current search direction,
     and the changes in the gradient of Q along S and the reduced D,
     respectively, where the reduced D is the same as D, except that the
     components of the fixed variables are zero.

     DSQ will be set to the square of the length of XNEW-XOPT.

     CRVMIN is set to zero if D reaches the trust region boundary.  Otherwise it
     is set to the least curvature of H that occurs in the conjugate gradient
     searches that are not restricted by any constraints.  The value
     CRVMIN=-1.0D0 is set, however, if all of these searches are constrained.

     A version of the truncated conjugate gradient is applied.  If a line search
     is restricted by a constraint, then the procedure is restarted, the values
     of the variables that are at their bounds being fixed.  If the trust region
     boundary is reached, then further changes may be made to D, each one being
     in the two dimensional space that is spanned by the current D and the
     gradient of Q at XOPT+D, staying on the trust region boundary.  Termination
     occurs when the reduction in Q seems to be close to the greatest reduction
     that can be achieved. */

  /* Constants. */
  const REAL half = 0.5;
  const REAL one = 1.0;
  const REAL onemin = -1.0;
  const REAL zero = 0.0;

  /* Local variables. */
  REAL angbd, angt, beta, blen, cth, delsq, dhd, dhs, dredg, dredsq, ds,
    ggsav, gredsq, qred, rdnext, rdprev, redmax, rednew, redsav, resid,
    sdec, shs, sredg, ssq, stepsq, sth, stplen, temp, tempa, tempb, xsav, xsum;
  INTEGER i, iact, ih, isav, itcsav, iterc, itermax, iu, j, k, nact;

  /* Parameter adjustments to comply with FORTRAN indexing. */
  xpt  -= 1 + npt;
  xopt -= 1;
  gopt -= 1;
  hq   -= 1;
  pq   -= 1;
  sl   -= 1;
  su   -= 1;
  xnew -= 1;
  d    -= 1;
  gnew -= 1;
  xbdi -= 1;
  s    -= 1;
  hs   -= 1;
  hred -= 1;
#define XPT(a1,a2) xpt[(a2)*npt + a1]

  /* FIXME: Set uninitialized variables. */
  angbd = zero;
  dredg = zero;
  dredsq = zero;
  ggsav = zero;
  gredsq = zero;
  rdnext = zero;
  sredg = zero;
  xsav = zero;
  iact = 0;
  itcsav = 0;
  itermax = 0;

  /* The sign of GOPT(I) gives the sign of the change to the I-th variable that
     will reduce Q from its value at XOPT.  Thus XBDI(I) shows whether or not
     to fix the I-th variable at one of its bounds initially, with NACT being
     set to the number of fixed variables.  D and GNEW are also set for the
     first iteration.  DELSQ is the upper bound on the sum of squares of the
     free variables.  QRED is the reduction in Q so far. */
  iterc = 0;
  nact = 0;
  LOOP(i,n) {
    xbdi[i] = zero;
    if (xopt[i] <= sl[i]) {
      if (gopt[i] >= zero) {
        xbdi[i] = onemin;
      }
    } else if (xopt[i] >= su[i]) {
      if (gopt[i] <= zero) {
        xbdi[i] = one;
      }
    }
    if (xbdi[i] != zero) {
      ++nact;
    }
    d[i] = zero;
    gnew[i] = gopt[i];
  }
  delsq = delta*delta;
  qred = zero;
  *crvmin = onemin;

  /* Set the next search direction of the conjugate gradient method.  It is the
     steepest descent direction initially and when the iterations are restarted
     because a variable has just been fixed by a bound, and of course the
     components of the fixed variables are zero.  ITERMAX is an upper bound on
     the indices of the conjugate gradient iterations. */
 L20:
  beta = zero;
 L30:
  stepsq = zero;
  LOOP(i,n) {
    if (xbdi[i] != zero) {
      s[i] = zero;
    } else if (beta == zero) {
      s[i] = -gnew[i];
    } else {
      s[i] = beta*s[i] - gnew[i];
    }
    stepsq += s[i]*s[i];
  }
  if (stepsq == zero) {
    goto L190;
  }
  if (beta == zero) {
    gredsq = stepsq;
    itermax = iterc + n - nact;
  }
  if (gredsq*delsq <= qred*1e-4*qred) {
    goto L190;
  }

  /* Multiply the search direction by the second derivative matrix of Q and
     calculate some scalars for the choice of steplength.  Then set BLEN to the
     length of the the step to the trust region boundary and STPLEN to the
     steplength, ignoring the simple bounds. */
  goto L210;
 L50:
  resid = delsq;
  ds = zero;
  shs = zero;
  LOOP(i,n) {
    if (xbdi[i] == zero) {
      resid -= d[i]*d[i];
      ds += s[i]*d[i];
      shs += s[i]*hs[i];
    }
  }
  if (resid <= zero) {
    goto L90;
  }
  temp = SQRT(stepsq*resid + ds*ds);
  if (ds < zero) {
    blen = (temp - ds)/stepsq;
  } else {
    blen = resid/(temp + ds);
  }
  if (shs > zero) {
    stplen = gredsq/shs;
    stplen = MIN(blen,stplen);
  } else {
    stplen = blen;
  }

  /* Reduce STPLEN if necessary in order to preserve the simple bounds, letting
     IACT be the index of the new constrained variable. */
  iact = 0;
  LOOP(i,n) {
    if (s[i] != zero) {
      xsum = xopt[i] + d[i];
      if (s[i] > zero) {
        temp = (su[i] - xsum)/s[i];
      } else {
        temp = (sl[i] - xsum)/s[i];
      }
      if (temp < stplen) {
        stplen = temp;
        iact = i;
      }
    }
  }

  /* Update CRVMIN, GNEW and D.  Set SDEC to the decrease that occurs in Q. */
  sdec = zero;
  if (stplen > zero) {
    ++iterc;
    temp = shs/stepsq;
    if (iact == 0 && temp > zero) {
      *crvmin = MIN(*crvmin,temp);
      if (*crvmin == onemin) {
        *crvmin = temp;
      }
    }
    ggsav = gredsq;
    gredsq = zero;
    LOOP(i,n) {
      gnew[i] += stplen*hs[i];
      if (xbdi[i] == zero) {
        gredsq += gnew[i]*gnew[i];
      }
      d[i] += stplen*s[i];
    }
    sdec = stplen*(ggsav - half*stplen*shs);
    sdec = MAX(sdec,zero);
    qred += sdec;
  }

  /* Restart the conjugate gradient method if it has hit a new bound. */
  if (iact > 0) {
    ++nact;
    xbdi[iact] = one;
    if (s[iact] < zero) {
      xbdi[iact] = onemin;
    }
    delsq -= d[iact]*d[iact];
    if (delsq <= zero) {
      goto L90;
    }
    goto L20;
  }

  /* If STPLEN is less than BLEN, then either apply another conjugate gradient
     iteration or RETURN. */
  if (stplen < blen) {
    if (iterc == itermax) {
      goto L190;
    }
    if (sdec <= qred*0.01) {
      goto L190;
    }
    beta = gredsq/ggsav;
    goto L30;
  }
 L90:
  *crvmin = zero;

  /* Prepare for the alternative iteration by calculating some scalars and by
     multiplying the reduced D by the second derivative matrix of Q. */
 L100:
  if (nact >= n - 1) {
    goto L190;
  }
  dredsq = zero;
  dredg = zero;
  gredsq = zero;
  LOOP(i,n) {
    if (xbdi[i] == zero) {
      dredsq += d[i]*d[i];
      dredg += d[i]*gnew[i];
      gredsq += gnew[i]*gnew[i];
      s[i] = d[i];
    } else {
      s[i] = zero;
    }
  }
  itcsav = iterc;
  goto L210;

  /* Let the search direction S be a linear combination of the reduced D and
     the reduced G that is orthogonal to the reduced D. */
 L120:
  ++iterc;
  temp = gredsq*dredsq - dredg*dredg;
  if (temp <= qred*1e-4*qred) {
    goto L190;
  }
  temp = SQRT(temp);
  LOOP(i,n) {
    if (xbdi[i] == zero) {
      s[i] = (dredg*d[i] - dredsq*gnew[i])/temp;
    } else {
      s[i] = zero;
    }
  }
  sredg = -temp;

  /* By considering the simple bounds on the variables, calculate an upper
     bound on the tangent of half the angle of the alternative iteration,
     namely ANGBD, except that, if already a free variable has reached a bound,
     there is a branch back to label 100 after fixing that variable. */
  angbd = one;
  iact = 0;
  LOOP(i,n) {
    if (xbdi[i] == zero) {
      tempa = xopt[i] + d[i] - sl[i];
      tempb = su[i] - xopt[i] - d[i];
      if (tempa <= zero) {
        ++nact;
        xbdi[i] = onemin;
        goto L100;
      } else if (tempb <= zero) {
        ++nact;
        xbdi[i] = one;
        goto L100;
      }
      ssq = d[i]*d[i] + s[i]*s[i];
      temp = xopt[i] - sl[i];
      temp = ssq - temp*temp;
      if (temp > zero) {
        temp = SQRT(temp) - s[i];
        if (angbd*temp > tempa) {
          angbd = tempa/temp;
          iact = i;
          xsav = onemin;
        }
      }
      temp = su[i] - xopt[i];
      temp = ssq - temp*temp;
      if (temp > zero) {
        temp = SQRT(temp) + s[i];
        if (angbd*temp > tempb) {
          angbd = tempb/temp;
          iact = i;
          xsav = one;
        }
      }
    }
  }
  goto L210;

  /* Calculate HHD and some curvatures for the alternative iteration. */
 L150:
  shs = zero;
  dhs = zero;
  dhd = zero;
  LOOP(i,n) {
    if (xbdi[i] == zero) {
      shs += s[i]*hs[i];
      dhs += d[i]*hs[i];
      dhd += d[i]*hred[i];
    }
  }

  /* Seek the greatest reduction in Q for a range of equally spaced values of
     ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of the
     alternative iteration. */
  redmax = zero;
  isav = 0;
  redsav = zero;
  iu = (INTEGER) (angbd*17.0 + 3.1);
  LOOP(i,iu) {
    angt = angbd*(REAL)i/(REAL)iu;
    sth = (angt + angt)/(one + angt*angt);
    temp = shs + angt*(angt*dhd - dhs - dhs);
    rednew = sth*(angt*dredg - sredg - half*sth*temp);
    if (rednew > redmax) {
      redmax = rednew;
      isav = i;
      rdprev = redsav;
    } else if (i == isav + 1) {
      rdnext = rednew;
    }
    redsav = rednew;
  }

  /* Return if the reduction is zero.  Otherwise, set the sine and cosine of
     the angle of the alternative iteration, and calculate SDEC. */
  if (isav == 0) {
    goto L190;
  }
  if (isav < iu) {
    temp = (rdnext - rdprev)/(redmax + redmax - rdprev - rdnext);
    angt = angbd*((REAL)isav + half*temp)/(REAL)iu;
  }
  cth = (one - angt*angt)/(one + angt*angt);
  sth = (angt + angt)/(one + angt*angt);
  temp = shs + angt*(angt*dhd - dhs - dhs);
  sdec = sth*(angt*dredg - sredg - half*sth*temp);
  if (sdec <= zero) {
    goto L190;
  }

  /* Update GNEW, D and HRED.  If the angle of the alternative iteration is
     restricted by a bound on a free variable, that variable is fixed at the
     bound. */
  dredg = zero;
  gredsq = zero;
  LOOP(i,n) {
    gnew[i] = gnew[i] + (cth - one)*hred[i] + sth*hs[i];
    if (xbdi[i] == zero) {
      d[i] = cth*d[i] + sth*s[i];
      dredg += d[i]*gnew[i];
      gredsq += gnew[i]*gnew[i];
    }
    hred[i] = cth*hred[i] + sth*hs[i];
  }
  qred += sdec;
  if (iact > 0 && isav == iu) {
    ++nact;
    xbdi[iact] = xsav;
    goto L100;
  }

  /* If SDEC is sufficiently small, then RETURN after setting XNEW to XOPT+D,
     giving careful attention to the bounds. */
  if (sdec > qred*0.01) {
    goto L120;
  }
 L190:
  *dsq = zero;
  LOOP(i,n) {
    temp = xopt[i] + d[i];
    temp = MIN(temp,su[i]);
    xnew[i] = MAX(temp,sl[i]);
    if (xbdi[i] == onemin) {
      xnew[i] = sl[i];
    }
    if (xbdi[i] == one) {
      xnew[i] = su[i];
    }
    d[i] = xnew[i] - xopt[i];
    *dsq += d[i]*d[i];
  }
  return;

  /* The following instructions multiply the current S-vector by the second
     derivative matrix of the quadratic model, putting the product in HS.  They
     are reached from three different parts of the software above and they can
     be regarded as an external subroutine. */
 L210:
  ih = 0;
  LOOP(j,n) {
    hs[j] = zero;
    LOOP(i,j) {
      ++ih;
      if (i < j) {
        hs[j] += hq[ih]*s[i];
      }
      hs[i] += hq[ih]*s[j];
    }
  }
  LOOP(k,npt) {
    if (pq[k] != zero) {
      temp = zero;
      LOOP(j,n) {
        temp += XPT(k,j)*s[j];
      }
      temp *= pq[k];
      LOOP(i,n) {
        hs[i] += temp*XPT(k,i);
      }
    }
  }
  if (*crvmin != zero) {
    goto L50;
  }
  if (iterc > itcsav) {
    goto L150;
  }
  LOOP(i,n) {
    hred[i] = hs[i];
  }
  goto L120;
} /* trsbox */

#undef XPT

static void
update(const INTEGER n, const INTEGER npt, REAL* bmat,
       REAL* zmat, const INTEGER ndim, REAL* vlag, const REAL beta,
       const REAL denom, const INTEGER knew, REAL* w)
{
  /* The arrays BMAT and ZMAT are updated, as required by the new position of
     the interpolation point that has the index KNEW.  The vector VLAG has
     N+NPT components, set on entry to the first NPT and last N components of
     the product Hw in equation (4.11) of the Powell (2006) paper on NEWUOA.
     Further, BETA is set on entry to the value of the parameter with that
     name, and DENOM is set to the denominator of the updating formula.
     Elements of ZMAT may be treated as zero if their moduli are at most ZTEST.
     The first NDIM elements of W are used for working space. */

  /* Constants. */
  const REAL one = 1.0;
  const REAL zero = 0.0;

  /* Local variables. */
  REAL alpha, tau, temp, tempa, tempb, ztest;
  INTEGER i, j, jp, k, nptm;

  /* Parameter adjustments to comply with FORTRAN indexing. */
  zmat -= 1 + npt;
  bmat -= 1 + ndim;
  vlag -= 1;
  w    -= 1;
#define BMAT(a1,a2) bmat[(a2)*ndim + a1]
#define ZMAT(a1,a2) zmat[(a2)*npt + a1]

  /* Initialization. */
  nptm = npt - n - 1;
  ztest = zero;
  LOOP(k,npt) {
    LOOP(j,nptm) {
      temp = ABS(ZMAT(k,j));
      ztest = MAX(ztest, temp);
    }
  }
  ztest *= 1e-20;

  /* Apply the rotations that put zeros in the KNEW-th row of ZMAT. */
  for (j = 2; j <= nptm; ++j) {
    if (ABS(ZMAT(knew,j)) > ztest) {
      tempa = ZMAT(knew,1);
      tempb = ZMAT(knew,j);
      temp = SQRT(tempa*tempa + tempb*tempb);
      tempa /= temp;
      tempb /= temp;
      LOOP(i,npt) {
        temp = tempa*ZMAT(i,1) + tempb*ZMAT(i,j);
        ZMAT(i,j) = tempa*ZMAT(i,j) - tempb*ZMAT(i,1);
        ZMAT(i,1) = temp;
      }
    }
    ZMAT(knew, j) = zero;
  }

  /* Put the first NPT components of the KNEW-th column of HLAG into W, and
     calculate the parameters of the updating formula. */
  LOOP(i,npt) {
    w[i] = ZMAT(knew,1)*ZMAT(i,1);
  }
  alpha = w[knew];
  tau = vlag[knew];
  vlag[knew] -= one;

  /* Complete the updating of ZMAT. */
  temp = SQRT(denom);
  tempb = ZMAT(knew,1)/temp;
  tempa = tau/temp;
  LOOP(i,npt) {
    ZMAT(i,1) = tempa*ZMAT(i,1) - tempb*vlag[i];
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
} /* update */

#undef ZMAT
#undef BMAT

static void
print_error(const char* reason)
{
  fprintf(stderr, "\n    Return from BOBYQA because %s.\n", reason);
}

static void
print_x(FILE* output, const INTEGER n, const REAL x[], const REAL dx[])
{
  INTEGER i, k;

  if (output == NULL) {
    output = stdout;
  }
  for (i = 0; i < n; ++i) {
    k = i%5;
    if (k == 0) {
      fprintf(output, "  ");
    }
    fprintf(output, "%15.6E", (dx == NULL ? x[i] : x[i] + dx[i]));
    if (i == n - 1 || k == 4) {
      fprintf(output, "\n");
    }
  }
}

/*---------------------------------------------------------------------------*/
