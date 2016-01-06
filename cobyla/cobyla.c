/*
 * cobyla.c -
 *
 * Implementation (in C) of Mike Powell's COBYLA algorithm for minimizing a
 * function of a few variables.  The method is "derivatives free" (only the
 * function values are needed) and accounts for constraints on the variables.
 * The algorithm is described in:
 *
 *   M.J.D. Powell, "A direct search optimization method that models the
 *   objective and constraint functions by linear interpolation," in Advances
 *   in Optimization and Numerical Analysis Mathematics and Its Applications,
 *   vol. 275 (eds. Susana Gomez and Jean-Pierre Hennart), Kluwer Academic
 *   Publishers, pp. 51-67 (1994).
 *
 * The present code is based on the original FORTRAN version written by Mike
 * Powell who kindly provides his code on demand and has been converted to C by
 * É. Thiébaut.
 *
 * Copyright (c) 1992, Mike Powell (FORTRAN version).
 * Copyright (c) 2015, Éric Thiébaut (C version).
 */

/* To re-use as much as the code for the reverse communication routines, we use
   a trick which consists in "self-including" this file with different macros
   (_COBYLA_PART1, _COBYLA_PART2, etc.) defined so as to skip or modify certain
   parts of the source file. */
#ifndef _COBYLA_PART1
#define _COBYLA_PART1 1

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>

#include "cobyla.h"

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

#define MAX(a, b) ((a) >= (b) ? (a) : (b))
#define MIN(a, b) ((a) <= (b) ? (a) : (b))

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

static int
cobylb(const INTEGER n, const INTEGER m,
       cobyla_calcfc* calcfc, void* calcfc_data, REAL x[],
       const REAL rhobeg, const REAL rhoend, const INTEGER iprint,
       INTEGER* maxfun, REAL con[], REAL sim[], REAL simi[],
       REAL datmat[], REAL a[], REAL vsig[], REAL veta[],
       REAL sigbar[], REAL dx[], REAL w[], INTEGER iact[]);

static void
trstlp(const INTEGER n, const INTEGER m, const REAL a[], const REAL b[],
       const REAL rho, REAL dx[], INTEGER* ifull, INTEGER iact[], REAL z[],
       REAL zdota[], REAL vmultc[], REAL sdirn[], REAL dxnew[], REAL vmultd[]);

static void
print_calcfc(FILE* output, INTEGER n, INTEGER nfvals,
             REAL f, REAL maxcv, const REAL x[]);


int
cobyla(INTEGER n, INTEGER m, cobyla_calcfc* calcfc, void* calcfc_data,
       REAL x[], REAL rhobeg, REAL rhoend, INTEGER iprint,
       INTEGER* maxfun, REAL w[], INTEGER iact[])
{
  /* Partition the working space array W to provide the storage that is needed
    for the main calculation. */
  INTEGER mpp = m + 2;
  REAL* con = w;
  REAL* sim = con + mpp;
  REAL* simi = sim + n*n + n;
  REAL* datmat = simi + n*n;
  REAL* a = datmat + n*mpp + mpp;
  REAL* vsig = a + m*n + n;
  REAL* veta = vsig + n;
  REAL* sigbar = veta + n;
  REAL* dx = sigbar + n;
  REAL* work = dx + n;
  return cobylb(n, m, calcfc, calcfc_data, x, rhobeg, rhoend, iprint, maxfun,
                con, sim, simi, datmat, a, vsig, veta, sigbar, dx, work, iact);
}

/*---------------------------------------------------------------------------*/
/* FORTRAN support. */

#ifdef FORTRAN_NAME

/* CALCFC_WRAPPER is a wrapper function, it assumes a CALCFC subroutine written
   in FORTRAN. */
static REAL
calcfc_wrapper(INTEGER n, INTEGER m, const REAL x[], REAL con[], void* data)
{
# define fc FORTRAN_NAME(calcfc,CALCFC)
  REAL f;
  extern void fc(INTEGER* n, INTEGER* m, const REAL x[], REAL* f, REAL con[]);
  fc(&n, &m, x, &f, con);
  return f;
# undef fc
}

int
FORTRAN_NAME(cobyla,COBYLA)(INTEGER* n, INTEGER* m, REAL x[], REAL* rhobeg,
                            REAL* rhoend, INTEGER* iprint, INTEGER* maxfun,
                            REAL w[], INTEGER iact[])
{
  (void)cobyla(*n, *m, calcfc_wrapper, NULL, x, *rhobeg, *rhoend,
               *iprint, maxfun, w, iact);
  return 0;
}

#endif /* FORTRAN_NAME */

/*---------------------------------------------------------------------------*/
/* Reverse communication version. */

struct _cobyla_context {
  INTEGER n;      /* number of variables */
  INTEGER m;      /* number of constraints */
  INTEGER iprint;
  INTEGER maxfun;
  INTEGER nfvals;
  REAL rhobeg;
  REAL rhoend;

  /* Workspace addresses. */
  INTEGER* iact;
  REAL* con;
  REAL* sim;
  REAL* simi;
  REAL* datmat;
  REAL* a;
  REAL* vsig;
  REAL* veta;
  REAL* sigbar;
  REAL* dx;
  REAL* w;

  REAL parmu, parsig, prerec, prerem, rho, f;
  INTEGER ibrnch, iflag, ifull, jdrop;

  int status;
};

cobyla_context_t*
cobyla_create(INTEGER n, INTEGER m, REAL rhobeg, REAL rhoend,
              INTEGER iprint, INTEGER maxfun)
{
  cobyla_context_t* ctx;
  long size, offset1, offset2;
  INTEGER mpp;

  /* Check arguments. */
  if (n < 1 || m < 0 || rhobeg < rhoend || rhoend <= 0 || maxfun < 1) {
    errno = EINVAL;
    return NULL;
  }

  /* Allocate memory. */
  size = sizeof(cobyla_context_t);
  offset1 = ROUND_UP(size, sizeof(INTEGER));
  size = offset1 + (m + 1)*sizeof(INTEGER);
  offset2 = ROUND_UP(size, sizeof(REAL));
  size = offset2 + (n*(3*n + 2*m + 11) + 4*m + 6)*sizeof(REAL);
  ctx = (cobyla_context_t*)malloc(size);
  if (ctx == NULL) {
    return NULL;
  }
  memset(ctx, 0, size);
  ctx->n = n;
  ctx->m = m;
  ctx->nfvals = 0; /* indicate a fresh start */
  ctx->status = COBYLA_ITERATE;
  ctx->iprint = iprint;
  ctx->maxfun = maxfun;
  ctx->rhobeg = rhobeg;
  ctx->rhoend = rhoend;

  /* Partition the working space array W to provide the storage that is needed
    for the main calculation. */
  mpp = m + 2;
  ctx->iact   = ADDRESS(INTEGER, ctx, offset1);
  ctx->con    = ADDRESS(REAL, ctx, offset2);
  ctx->sim    = ctx->con + mpp;
  ctx->simi   = ctx->sim + n*n + n;
  ctx->datmat = ctx->simi + n*n;
  ctx->a      = ctx->datmat + n*mpp + mpp;
  ctx->vsig   = ctx->a + m*n + n;
  ctx->veta   = ctx->vsig + n;
  ctx->sigbar = ctx->veta + n;
  ctx->dx     = ctx->sigbar + n;
  ctx->w      = ctx->dx + n;

  return ctx;
}

void
cobyla_delete(cobyla_context_t* ctx)
{
  if (ctx != NULL) {
    free((void*)ctx);
  }
}

int
cobyla_restart(cobyla_context_t* ctx)
{
  if (ctx == NULL) {
    errno = EFAULT;
    return COBYLA_BAD_ADDRESS;
  }
  ctx->nfvals = 0;
  ctx->status = COBYLA_ITERATE;
  return ctx->status;
}

int
cobyla_get_status(const cobyla_context_t* ctx)
{
  if (ctx == NULL) {
    errno = EFAULT;
    return COBYLA_BAD_ADDRESS;
  }
  return ctx->status;
}

INTEGER
cobyla_get_nevals(const cobyla_context_t* ctx)
{
  if (ctx == NULL) {
    errno = EFAULT;
    return -1;
  }
  return ctx->nfvals;
}

REAL
cobyla_get_rho(const cobyla_context_t* ctx)
{
  if (ctx == NULL) {
    errno = EFAULT;
    return -1;
  }
  return ctx->rho;
}

REAL
cobyla_get_last_f(const cobyla_context_t* ctx)
{
  if (ctx == NULL) {
    errno = EFAULT;
    return -1;
  }
  return ctx->f;
}

/* Include this file with the macro _COBYLA_REVCOM defined to
   generate the code `newuoa_iterate` in place of `newuob`. */
#define _COBYLA_REVCOM 1
#include __FILE__
#undef _COBYLA_REVCOM

#endif /* _COBYLA_PART1 */

/* Define macros mimicking FORTRAN indexing. */
#define      A(a1,a2)      a[a1 - 1 + n*(a2 - 1)]
#define    CON(a1)       con[a1 - 1]
#define DATMAT(a1,a2) datmat[a1 - 1 + mpp*(a2 - 1)]
#define     DX(a1)        dx[a1 - 1]
#define   IACT(a1)      iact[a1 - 1]
#define SIGBAR(a1)    sigbar[a1 - 1]
#define    SIM(a1,a2)    sim[a1 - 1 + n*(a2 - 1)]
#define   SIMI(a1,a2)   simi[a1 - 1 + n*(a2 - 1)]
#define   VETA(a1)      veta[a1 - 1]
#define   VSIG(a1)      vsig[a1 - 1]
#define      W(a1)         w[a1 - 1]
#define      X(a1)         x[a1 - 1]

#ifdef _COBYLA_REVCOM
#  define RESTORE(var)  var = ctx->var
#  define SAVE(var)     ctx->var = var
#endif

#ifdef _COBYLA_REVCOM
int
cobyla_iterate(cobyla_context_t* ctx, REAL f, REAL x[], REAL c[])
#else
static int
cobylb(const INTEGER n, const INTEGER m,
       cobyla_calcfc* calcfc, void* calcfc_data, REAL x[],
       const REAL rhobeg, const REAL rhoend, const INTEGER iprint,
       INTEGER* _maxfun, REAL con[], REAL sim[], REAL simi[],
       REAL datmat[], REAL a[], REAL vsig[], REAL veta[],
       REAL sigbar[], REAL dx[], REAL w[], INTEGER iact[])
#endif
{
  /* Constants. */
  const REAL zero  = FLT(0.0);
  const REAL one   = FLT(1.0);
  const REAL alpha = FLT(0.25);
  const REAL beta  = FLT(2.1);
  const REAL gamma = FLT(0.5);
  const REAL delta = FLT(1.1);

  /* Local variables (for each type, real or integer, the first group of
     variables are those which must be preserved between calls in a reverse
     communication version of the algorithm). */
  REAL parmu, parsig, prerec, prerem, rho;
  REAL barmu, cvmaxm, cvmaxp, dxsign, edgmax, error, pareta, phi, phimin,
    ratio, resmax, resnew, sum, temp, tempa, trured, vmnew, vmold;
  INTEGER ibrnch, iflag, ifull, jdrop, nfvals, maxfun;
  INTEGER i, j, k, l, mp, mpp, np, nbest;
#ifdef _COBYLA_REVCOM
  INTEGER n, m,  iprint;
  REAL rhobeg, rhoend;
  REAL *con, *sim, *simi, *datmat, *a, *vsig, *veta, *sigbar, *dx, *w;
  INTEGER *iact;
#else
  REAL f;
#endif
  int status;

#ifdef _COBYLA_REVCOM
  /* Minimal checking and restore initial set of variables. */
  if (ctx == NULL) {
    errno = EFAULT;
    return COBYLA_BAD_ADDRESS;
  }
  RESTORE(n);
  RESTORE(m);
  RESTORE(iprint);
  RESTORE(maxfun);
  RESTORE(nfvals);
  RESTORE(rhobeg);
  RESTORE(rhoend);
  RESTORE(iact);
  RESTORE(con);
  RESTORE(sim);
  RESTORE(simi);
  RESTORE(datmat);
  RESTORE(a);
  RESTORE(vsig);
  RESTORE(veta);
  RESTORE(sigbar);
  RESTORE(dx);
  RESTORE(w);
  RESTORE(status);
  if (x == NULL || (c == NULL && m > 0)) {
    errno = EFAULT;
    ctx->status = COBYLA_BAD_ADDRESS;
    return COBYLA_BAD_ADDRESS;
  }
  if (status != COBYLA_ITERATE || nfvals < 0) {
    errno = EINVAL;
    ctx->status = COBYLA_CORRUPTED;
    return ctx->status;
  }
#endif

  /* Set the initial values of some parameters.  The last column of SIM holds
     the optimal vertex of the current simplex, and the preceding N columns
     hold the displacements from the optimal vertex to the other vertices.
     Further, SIMI holds the inverse of the matrix that is contained in the
     first N columns of SIM. */
  np = n + 1;
  mp = m + 1;
  mpp = m + 2;
#ifdef _COBYLA_REVCOM
  if (nfvals == 0) {
    /* This is the first function evaluation.  Proceed with initialization. */
    status = COBYLA_INITIAL_ITERATE;
    prerec = zero; /* avoids warnings */
    prerem = zero; /* avoids warnings */
    parsig = zero; /* avoids warnings */
    iflag  = 0;    /* avoids warnings */
#else
    maxfun = *_maxfun;
    nfvals = 0;
#endif
    rho = rhobeg;
    parmu = zero;
    jdrop = np;
    ibrnch = 0;
    temp = one/rho;
    LOOP(i,n) {
      SIM(i,np) = X(i);
      LOOP(j,n) {
        SIM(i,j) = zero;
        SIMI(i,j) = zero;
      }
      SIM(i,i) = rho;
      SIMI(i,i) = temp;
    }
    if (iprint >= 2) {
      fprintf(stdout, "\n   The initial value of RHO is%13.6E"
              "  and PARMU is set to zero.\n", (double)rho);
    }
#ifdef _COBYLA_REVCOM
  } else {
    /* This is not the first function evaluation.  Restore other local
       variables and jump to the place where function value is expected. */
    RESTORE(parmu);
    RESTORE(parsig);
    RESTORE(prerec);
    RESTORE(prerem);
    RESTORE(rho);
    RESTORE(ibrnch);
    RESTORE(iflag);
    RESTORE(ifull);
    RESTORE(jdrop);
    goto new_eval;
  }
#endif

  /* Make the next call of the user-supplied subroutine CALCFC.  These
     instructions are also used for calling CALCFC during the iterations
     of the algorithm. */
 call_calcfc:
  if (nfvals >= maxfun && nfvals > 0) {
    status = COBYLA_TOO_MANY_EVALUATIONS;
    fprintf(stderr, "Return from subroutine COBYLA because %s.\n",
            cobyla_reason(status));
    goto L_600;
  }
#ifdef _COBYLA_REVCOM
  if (status == COBYLA_INITIAL_ITERATE) {
    /* We already know the functiuon value. */
    status = COBYLA_ITERATE;
  } else {
    /* Request a new function evaluation from the caller and branch to save
       local variables in the context structure. */
    status = COBYLA_ITERATE;
    goto save;
  }
  /* We arrive here when caller has computed a new function value. */
 new_eval:
#else
  f = calcfc(n, m, x, con, calcfc_data);
#endif
  ++nfvals;

  /* Estimate the worst constraint RESMAX. */
  resmax = zero;
#ifdef _COBYLA_REVCOM
  for (k = 0; k < m; ++k) {
    temp = c[k];
    con[k] = temp;
    if (resmax < -temp) {
      resmax = -temp;
    }
  }
#else
  for (k = 0; k < m; ++k) {
    if (resmax < (temp = -con[k])) {
      resmax = temp;
    }
  }
#endif
  if (nfvals == iprint - 1 || iprint == 3) {
    print_calcfc(stdout, n, nfvals, f, resmax, &X(1));
  }
  CON(mp) = f;
  CON(mpp) = resmax;
  if (ibrnch == 1) goto L_440;

  /* Set the recently calculated function values in a column of DATMAT.  This
     array has a column for each vertex of the current simplex, the entries of
     each column being the values of the constraint functions (if any) followed
     by the objective function and the greatest constraint violation at the
     vertex. */
  LOOP(k,mpp) {
    DATMAT(k,jdrop) = CON(k);
  }
  if (nfvals > np) goto L_130;

  /* Exchange the new vertex of the initial simplex with the optimal vertex if
     necessary.  Then, if the initial simplex is not complete, pick its next
     vertex and calculate the function values there. */
  if (jdrop <= n) {
    if (DATMAT(mp,np) <= f) {
      X(jdrop) = SIM(jdrop,np);
    } else {
      SIM(jdrop,np) = X(jdrop);
      LOOP(k,mpp) {
        DATMAT(k,jdrop) = DATMAT(k,np);
        DATMAT(k,np) = CON(k);
      }
      LOOP(k,jdrop) {
        temp = zero;
        SIM(jdrop,k) = -rho;
        for (i = k; i <= jdrop; ++i) {
          temp -= SIMI(i,k);
        }
        SIMI(jdrop,k) = temp;
      }
    }
  }
  if (nfvals <= n) {
    jdrop = nfvals;
    X(jdrop) += rho;
    goto call_calcfc;
  }
 L_130:
  ibrnch = 1;

  /* Identify the optimal vertex of the current simplex. */
 L_140:
  phimin = DATMAT(mp,np) + parmu*DATMAT(mpp,np);
  nbest = np;
  LOOP(j,n) {
    temp = DATMAT(mp,j) + parmu*DATMAT(mpp,j);
    if (temp < phimin) {
      nbest = j;
      phimin = temp;
    } else if (temp == phimin && parmu == zero) {
      if (DATMAT(mpp,j) < DATMAT(mpp,nbest)) nbest = j;
    }
  }

  /* Switch the best vertex into pole position if it is not there already,
     and also update SIM, SIMI and DATMAT. */
  if (nbest <= n) {
    LOOP(i,mpp) {
      temp = DATMAT(i,np);
      DATMAT(i,np) = DATMAT(i,nbest);
      DATMAT(i,nbest) = temp;
    }
    LOOP(i,n) {
      temp = SIM(i,nbest);
      SIM(i,nbest) = zero;
      SIM(i,np) += temp;
      tempa = zero;
      LOOP(k,n) {
        SIM(i,k) -= temp;
        tempa -= SIMI(k,i);
      }
      SIMI(nbest,i) = tempa;
    }
  }

  /* Make an error return if SIGI is a poor approximation to the inverse of
     the leading N by N submatrix of SIG. */
  error = zero;
  LOOP(i,n) {
    LOOP(j,n) {
      temp = (i == j ? -one : zero);
      LOOP(k,n) {
        temp += SIMI(i,k)*SIM(k,j);
      }
      temp = ABS(temp);
      error = MAX(error, temp);
    }
  }
  if (error > FLT(0.1)) {
    status = COBYLA_ROUNDING_ERRORS;
    if (iprint >= 1) {
      fprintf(stderr, "Return from subroutine COBYLA because %s.\n",
              cobyla_reason(status));
    }
    goto L_600;
  }

  /* Calculate the coefficients of the linear approximations to the objective
     and constraint functions, placing minus the objective function gradient
     after the constraint gradients in the array A.  The vector W is used for
     working space. */
  LOOP(k,mp) {
    CON(k) = -DATMAT(k,np);
    LOOP(j,n) {
      W(j) = DATMAT(k,j) + CON(k);
    }
    LOOP(i,n) {
      temp = zero;
      LOOP(j,n) {
        temp += W(j)*SIMI(j,i);
      }
      A(i,k) = (k == mp ? -temp : temp);
    }
  }

  /* Calculate the values of SIGMA and ETA, and set IFLAG=0 if the current
     simplex is not acceptable. */
  iflag = 1;
  parsig = alpha*rho;
  pareta = beta*rho;
  LOOP(j,n) {
    REAL wsig = zero;
    REAL weta = zero;
    LOOP(i,n) {
      wsig += SIMI(j,i)*SIMI(j,i);
      weta += SIM(i,j)*SIM(i,j);
    }
    VSIG(j) = one/SQRT(wsig);
    VETA(j) = SQRT(weta);
    if (VSIG(j) < parsig || VETA(j) > pareta) iflag = 0;
  }

  /* If a new vertex is needed to improve acceptability, then decide which
     vertex to drop from the simplex. */
  if (ibrnch == 1 || iflag == 1) goto L_370;
  jdrop = 0;
  temp = pareta;
  LOOP(j,n) {
    if (VETA(j) > temp) {
      jdrop = j;
      temp = VETA(j);
    }
  }
  if (jdrop == 0) {
    LOOP(j,n) {
      if (VSIG(j) < temp) {
        jdrop = j;
        temp = VSIG(j);
      }
    }
  }

  /* Calculate the step to the new vertex and its sign. */
  temp = gamma*rho*VSIG(jdrop);
  LOOP(i,n) {
    DX(i) = temp*SIMI(jdrop,i);
  }
  cvmaxp = zero;
  cvmaxm = zero;
  sum = zero; /* FIXME: `sum` was uninitialized */
  LOOP(k,mp) {
    sum = zero;
    LOOP(i,n) {
      sum = sum + A(i,k)*DX(i);
    }
    if (k < mp) {
      temp = DATMAT(k,np);
      cvmaxp = MAX(cvmaxp, -sum - temp);
      cvmaxm = MAX(cvmaxm,  sum - temp);
    }
  }
  if (parmu*(cvmaxp - cvmaxm) > sum + sum) {
    dxsign = -one;
  } else {
    dxsign = one;
  }

  /* Update the elements of SIM and SIMI, and set the next X. */
  temp = zero;
  LOOP(i,n) {
    DX(i) *= dxsign;
    SIM(i,jdrop) = DX(i);
    temp += SIMI(jdrop,i)*DX(i);
  }
  LOOP(i,n) {
    SIMI(jdrop,i) /= temp;
  }
  LOOP(j,n) {
    if (j != jdrop) {
      temp = zero;
      LOOP(i,n) {
        temp += SIMI(j,i)*DX(i);
      }
      LOOP(i,n) {
        SIMI(j,i) -= temp*SIMI(jdrop,i);
      }
    }
    X(j) = SIM(j,np) + DX(j);
  }
  goto call_calcfc;

  /* Calculate DX = X(*) - X(0).  Branch if the length of DX is less than
     0.5*RHO. */
 L_370:
  {
    REAL* z = w;
    REAL* zdota = z + n*n;
    REAL* vmc = zdota + n;
    REAL* sdirn = vmc + mp;
    REAL* dxnew = sdirn + n;
    REAL* vmd = dxnew + n;
    trstlp(n, m, a, con, rho, dx, &ifull, iact, z, zdota,
           vmc, sdirn, dxnew, vmd);
  }
  if (ifull == 0) {
    temp = zero;
    LOOP(i,n) {
      temp += DX(i)*DX(i);
    }
    if (temp < FLT(0.25)*rho*rho) {
      ibrnch = 1;
      goto L_550;
    }
  }

  /* Predict the change to F and the new maximum constraint violation if the
     variables are altered from X(0) to X(0)+DX. */
  resnew = zero;
  CON(mp) = zero;
  sum = zero; /* FIXME: `sum` was uninitialized */
  LOOP(k,mp) {
    sum = CON(k);
    LOOP(i,n) {
      sum -= A(i,k)*DX(i);
    }
    if (k < mp) resnew = MAX(resnew, sum);
  }

  /* Increase PARMU if necessary and branch back if this change alters the
     optimal vertex.  Otherwise PREREM and PREREC will be set to the predicted
     reductions in the merit function and the maximum constraint violation
     respectively. */
  prerec = DATMAT(mpp,np) - resnew;
  barmu = (prerec > zero ? sum/prerec : zero);
  if (parmu < FLT(1.5)*barmu) {
    parmu = FLT(2.0)*barmu;
    if (iprint >= 2) {
      fprintf(stdout, "\n   Increase in PARMU to%13.6E\n", (double)parmu);
    }
    phi = DATMAT(mp,np) + parmu*DATMAT(mpp,np);
    LOOP(j,n) {
      temp = DATMAT(mp,j) + parmu*DATMAT(mpp,j);
      if (temp < phi) goto L_140;
      if (temp == phi && parmu == zero) {
        if (DATMAT(mpp,j) < DATMAT(mpp,np)) goto L_140;
      }
    }
  }
  prerem = parmu*prerec - sum;

  /* Calculate the constraint and objective functions at X(*).  Then find the
     actual reduction in the merit function. */
  LOOP(i,n) {
    X(i) = SIM(i,np) + DX(i);
  }
  ibrnch = 1;
  goto call_calcfc;
 L_440:
  vmold = DATMAT(mp,np) + parmu*DATMAT(mpp,np);
  vmnew = f + parmu*resmax;
  trured = vmold - vmnew;
  if (parmu == zero && f == DATMAT(mp,np)) {
    prerem = prerec;
    trured = DATMAT(mpp,np) - resmax;
  }

  /* Begin the operations that decide whether X(*) should replace one of the
     vertices of the current simplex, the change being mandatory if TRURED is
     positive.  Firstly, JDROP is set to the index of the vertex that is to be
     replaced. */
  ratio = (trured <= zero ? one : zero);
  jdrop = 0;
  LOOP(j,n) {
    temp = zero;
    LOOP(i,n) {
      temp += SIMI(j,i)*DX(i);
    }
    temp = ABS(temp);
    if (temp > ratio) {
      jdrop = j;
      ratio = temp;
    }
    SIGBAR(j) = temp*VSIG(j);
  }

  /* Calculate the value of ell. */
  edgmax = delta*rho;
  l = 0;
  LOOP(j,n) {
    if (SIGBAR(j) >= parsig || SIGBAR(j) >= VSIG(j)) {
      temp = VETA(j);
      if (trured > zero) {
        temp = zero;
        LOOP(i,n) {
          REAL tempb = DX(i) - SIM(i,j);
          temp += tempb*tempb;
        }
        temp = SQRT(temp);
      }
      if (temp > edgmax) {
        l = j;
        edgmax = temp;
      }
    }
  }
  if (l > 0) jdrop = l;
  if (jdrop == 0) goto L_550;

  /* Revise the simplex by updating the elements of SIM, SIMI and DATMAT. */
  temp = zero;
  LOOP(i,n) {
    SIM(i,jdrop) = DX(i);
    temp += SIMI(jdrop,i)*DX(i);
  }
  LOOP(i,n) {
    SIMI(jdrop,i) = SIMI(jdrop,i)/temp;
  }
  LOOP(j,n) {
    if (j != jdrop) {
      temp = zero;
      LOOP(i,n) {
        temp += SIMI(j,i)*DX(i);
      }
      LOOP(i,n) {
        SIMI(j,i) -= temp*SIMI(jdrop,i);
      }
    }
  }
  LOOP(k,mpp) {
    DATMAT(k,jdrop) = CON(k);
  }

  /* Branch back for further iterations with the current RHO. */
  if (trured > zero && trured >= FLT(0.1)*prerem) goto L_140;
 L_550:
  if (iflag == 0) {
    ibrnch = 0;
    goto L_140;
  }

  /* Otherwise reduce RHO if it is not at its least value and reset PARMU. */
  if (rho > rhoend) {
    rho = FLT(0.5)*rho;
    if (rho <= FLT(1.5)*rhoend) rho = rhoend;
    if (parmu > zero) {
      REAL cmin, cmax, denom;
      denom = zero;
      LOOP(k,mp) {
        cmax = cmin = DATMAT(k,np);
        LOOP(i,n) {
          temp = DATMAT(k,i);
          cmin = MIN(cmin, temp);
          cmax = MAX(cmax, temp);
        }
        if (k <= m && cmin < FLT(0.5)*cmax) {
          temp = MAX(cmax, zero) - cmin;
          if (denom <= zero) {
            denom = temp;
          } else {
            denom = MIN(denom, temp);
          }
        }
      }
      if (denom == zero) {
        parmu = zero;
      } else if (cmax - cmin < parmu*denom) {
        parmu = (cmax - cmin)/denom;
      }
    }
    if (iprint >= 2) {
      fprintf(stdout, "\n   Reduction in RHO to%13.6E  and PARMU =%13.6E\n",
              (double)rho, (double)parmu);
      if (iprint == 2) {
        print_calcfc(stdout, n, nfvals, DATMAT(mp,np), DATMAT(mpp,np),
                     &SIM(1,np));
      }
    }
    goto L_140;
  }

  /* Return the best calculated values of the variables. */
  if (iprint >= 1) {
    fprintf(stdout, "\n   Normal return from subroutine COBYLA\n");
  }
  status = COBYLA_SUCCESS;
  if (ifull == 1) goto L_620;
 L_600:
  LOOP(i,n) {
    X(i) = SIM(i,np);
  }
  f = DATMAT(mp,np);
  resmax = DATMAT(mpp,np);
 L_620:
  if (iprint >= 1) {
    print_calcfc(stdout, n, nfvals, f, resmax, &X(1));
  }

#ifdef _COBYLA_REVCOM
  /* Save local variables and return status. */
 save:
  SAVE(nfvals);
  SAVE(parmu);
  SAVE(parsig);
  SAVE(prerec);
  SAVE(prerem);
  SAVE(rho);
  SAVE(f);
  SAVE(ibrnch);
  SAVE(iflag);
  SAVE(ifull);
  SAVE(jdrop);
  SAVE(status);
#else
  /* Save number of function evaluations and return status. */
  *_maxfun = nfvals;
#endif
  return status;
}

/* Undefine macros mimicking FORTRAN indexing. */
#undef A
#undef CON
#undef DATMAT
#undef DX
#undef IACT
#undef SIGBAR
#undef SIM
#undef SIMI
#undef VETA
#undef VSIG
#undef W
#undef X

#ifdef _COBYLA_REVCOM
#  undef RESTORE
#  undef SAVE
#endif

#ifndef _COBYLA_PART2
#define _COBYLA_PART2 1

/*
 * This  subroutine  calculates  an  N-component  vector  DX  by  applying  the
 * following two stages.  In the first stage,  DX is set to the shortest vector
 * that minimizes the greatest violation of the constraints:
 *
 *   A(1,K)*DX(1)+A(2,K)*DX(2)+...+A(N,K)*DX(N)  >=  B(K), K=2,3,...,M,
 *
 * subject to the Euclidean  length of DX being at most RHO.   If its length is
 * strictly less than RHO, then we use  the resultant freedom in DX to minimize
 * the objective function:
 *
 *            -A(1,M+1)*DX(1)-A(2,M+1)*DX(2)-...-A(N,M+1)*DX(N)
 *
 * subject to no increase in  any greatest constraint violation.  This notation
 * allows the gradient of the objective function to be regarded as the gradient
 * of a constraint.  Therefore the two stages are distinguished by MCON = M and
 * MCON > M respectively.  It is possible that a degeneracy may prevent DX from
 * attaining the target  length RHO.  Then the value IFULL=0  would be set, but
 * usually IFULL=1 on return.
 *
 * In  general  NACT  is the  number  of  constraints  in  the active  set  and
 * IACT(1),...,IACT(NACT)  are  their  indices,  while the  remainder  of  IACT
 * contains a permutation  of the remaining constraint indices.   Further, Z is
 * an orthogonal matrix whose first NACT  columns can be regarded as the result
 * of   Gram-Schmidt  applied   to  the   active  constraint   gradients.   For
 * J=1,2,...,NACT, the number ZDOTA(J) is the scalar product of the J-th column
 * of Z  with the gradient  of the J-th active  constraint.  DX is  the current
 * vector of variables and here the  residuals of the active constraints should
 * be  zero.   Further,  the   active  constraints  have  nonnegative  Lagrange
 * multipliers that are held at the beginning of VMULTC.  The remainder of this
 * vector holds the  residuals of the inactive constraints at  DX, the ordering
 * of the components  of VMULTC being in agreement with  the permutation of the
 * indices  of the  constraints  that  is in  IACT.   All  these residuals  are
 * nonnegative, which  is achieved  by the  shift RESMAX  that makes  the least
 * residual zero.
 *
 * Initialize  Z  and some  other  variables.   The  value  of RESMAX  will  be
 * appropriate  to DX=0,  while  ICON will  be  the index  of  a most  violated
 * constraint if RESMAX is positive.  Usually during the first stage the vector
 * SDIRN  gives a  search  direction  that reduces  all  the active  constraint
 * violations by one simultaneously.
 */

/* Define macros mimicking FORTRAN indexing. */
#define      A(a1,a2)      a[a1 - 1 + n*(a2 - 1)]
#define      B(a1)         b[a1 - 1]
#define     DX(a1)        dx[a1 - 1]
#define  DXNEW(a1)     dxnew[a1 - 1]
#define   IACT(a1)      iact[a1 - 1]
#define  SDIRN(a1)     sdirn[a1 - 1]
#define VMULTC(a1)    vmultc[a1 - 1]
#define VMULTD(a1)    vmultd[a1 - 1]
#define      Z(a1,a2)      z[a1 - 1 + n*(a2 - 1)]
#define  ZDOTA(a1)     zdota[a1 - 1]

static void
trstlp(const INTEGER n, const INTEGER m, const REAL a[], const REAL b[],
       const REAL rho, REAL dx[], INTEGER* ifull, INTEGER iact[], REAL z[],
       REAL zdota[], REAL vmultc[], REAL sdirn[], REAL dxnew[], REAL vmultd[])
{
  /* Constants. */
  const REAL zero  = FLT(0.0);
  const REAL one   = FLT(1.0);
  const REAL tiny  = FLT(1.0e-6);
  const REAL Op1   = FLT(0.1);
  const REAL Op2   = FLT(0.2);

  /* Local variables (FIXME: `iout` was unused). */
  REAL acca, accb, alpha, beta, dd, optnew, optold, ratio, resmax, resold, sd,
    sp, spabs, ss, step, stpful, sum, sumabs, temp, tempa, tempb, tot, vsave,
    zdotv, zdotw, zdvabs, zdwabs;
  INTEGER i, icon, icount, isave, j, k, kk, kl, kp, kw, mcon, nact, nactx;

  *ifull = 1;
  icon = 0; /* FIXME: `icon` was uninitialized */
  mcon = m;
  nact = 0;
  resmax = zero;
  resold = zero; /* FIXME: `resold` was uninitialized */
  LOOP(i,n) {
    LOOP(j,n) {
      Z(i,j) = zero;
    }
    Z(i,i) = one;
    DX(i) = zero;
  }
  if (m >= 1) {
    LOOP(k,m) {
      if (B(k) > resmax) {
        resmax = B(k);
        icon = k;
      }
    }
    LOOP(k,m) {
      IACT(k) = k;
      VMULTC(k) = resmax - B(k);
    }
  }
  if (resmax == zero) goto L_480;
  LOOP(i,n) {
    SDIRN(i) = zero;
  }

  /* End the current stage of the calculation if 3 consecutive iterations have
     either failed to reduce the best calculated value of the objective
     function or to increase the number of active constraints since the best
     value was calculated.  This strategy prevents cycling, but there is a
     remote possibility that it will cause premature termination. */
 L_60:
  optold = zero;
  icount = 0;
 L_70:
  if (mcon == m) {
    optnew = resmax;
  } else {
    optnew = zero;
    LOOP(i,n) {
      optnew -= DX(i)*A(i,mcon);
    }
  }
  if (icount == 0 || optnew < optold) {
    optold = optnew;
    nactx = nact;
    icount = 3;
  } else if (nact > nactx) {
    nactx = nact;
    icount = 3;
  } else {
    --icount;
    if (icount == 0) goto L_490;
  }

  /* If ICON exceeds NACT, then we add the constraint with index IACT(ICON) to
     the active set.  Apply Givens rotations so that the last N-NACT-1 columns
     of Z are orthogonal to the gradient of the new constraint, a scalar
     product being set to zero if its nonzero value could be due to computer
     rounding errors.  The array DXNEW is used for working space. */
  if (icon <= nact) goto L_260;
  kk = IACT(icon);
  LOOP(i,n) {
    DXNEW(i) = A(i,kk);
  }
  tot = zero;
  for (k = n; k > nact; --k) {
    sp = zero;
    spabs = zero;
    LOOP(i,n) {
      temp = Z(i,k)*DXNEW(i);
      sp += temp;
      spabs += ABS(temp);
    }
    acca = spabs + Op1*ABS(sp);
    accb = spabs + Op2*ABS(sp);
    if (spabs >= acca || acca >= accb) sp = zero;
    if (tot == zero) {
      tot = sp;
    } else {
      kp = k + 1;
      temp = SQRT(sp*sp + tot*tot);
      alpha = sp/temp;
      beta = tot/temp;
      tot = temp;
      LOOP(i,n) {
        temp = alpha*Z(i,k) + beta*Z(i,kp);
        Z(i,kp) = alpha*Z(i,kp) - beta*Z(i,k);
        Z(i,k) = temp;
      }
    }
  }

  /* Add the new constraint if this can be done without a deletion from the
     active set. */
  if (tot != zero) {
    ++nact;
    ZDOTA(nact) = tot;
    VMULTC(icon) = VMULTC(nact);
    VMULTC(nact) = zero;
    goto L_210;
  }

  /* The next instruction is reached if a deletion has to be made from the
     active set in order to make room for the new active constraint, because
     the new constraint gradient is a linear combination of the gradients of
     the old active constraints.  Set the elements of VMULTD to the multipliers
     of the linear combination.  Further, set IOUT to the index of the
     constraint to be deleted, but branch if no suitable index can be found. */
  ratio = -one;
  k = nact;
 L_130:
  zdotv = zero;
  zdvabs = zero;
  LOOP(i,n) {
    temp = Z(i,k)*DXNEW(i);
    zdotv += temp;
    zdvabs += ABS(temp);
  }
  acca = zdvabs + Op1*ABS(zdotv);
  accb = zdvabs + Op2*ABS(zdotv);
  if (zdvabs < acca && acca < accb) {
    temp = zdotv/ZDOTA(k);
    if (temp > zero && IACT(k) <= m) {
      tempa = VMULTC(k)/temp;
      if (ratio < zero || tempa < ratio) {
        ratio = tempa;
        /* iout = k; (FIXME: `iout` unused) */
      }
    }
    if (k >= 2) {
      kw = IACT(k);
      LOOP(i,n) {
        DXNEW(i) -= temp*A(i,kw);
      }
    }
    VMULTD(k) = temp;
  } else {
    VMULTD(k) = zero;
  }
  --k;
  if (k > 0) goto L_130;
  if (ratio < zero) goto L_490;

  /* Revise the Lagrange multipliers and reorder the active constraints so that
     the one to be replaced is at the end of the list.  Also calculate the new
     value of ZDOTA(NACT) and branch if it is not acceptable. */
  LOOP(k,nact) {
    tempb = VMULTC(k) - ratio*VMULTD(k);
    VMULTC(k) = MAX(zero, tempb);
  }
  if (icon < nact) {
    isave = IACT(icon);
    vsave = VMULTC(icon);
    k = icon;
    do {
      kp = k + 1;
      kw = IACT(kp);
      sp = zero;
      LOOP(i,n) {
        sp += Z(i,k)*A(i,kw);
      }
      temp = SQRT(sp*sp + ZDOTA(kp)*ZDOTA(kp));
      alpha = ZDOTA(kp)/temp;
      beta = sp/temp;
      ZDOTA(kp) = alpha*ZDOTA(k);
      ZDOTA(k) = temp;
      LOOP(i,n) {
        temp = alpha*Z(i,kp) + beta*Z(i,k);
        Z(i,kp) = alpha*Z(i,k) - beta*Z(i,kp);
        Z(i,k) = temp;
      }
      IACT(k) = kw;
      VMULTC(k) = VMULTC(kp);
      k = kp;
    } while (k < nact);
    IACT(k) = isave;
    VMULTC(k) = vsave;
  }
  temp = zero;
  LOOP(i,n) {
    temp += Z(i,nact)*A(i,kk);
  }
  if (temp == zero) goto L_490;
  ZDOTA(nact) = temp;
  VMULTC(icon) = zero;
  VMULTC(nact) = ratio;

  /* Update IACT and ensure that the objective function continues to be treated
     as the last active constraint when MCON > M. */
 L_210:
  IACT(icon) = IACT(nact);
  IACT(nact) = kk;
  if (mcon > m && kk != mcon) {
    k = nact - 1;
    sp = zero;
    LOOP(i,n) {
      sp += Z(i,k)*A(i,kk);
    }
    temp = SQRT(sp*sp + ZDOTA(nact)*ZDOTA(nact));
    alpha = ZDOTA(nact)/temp;
    beta = sp/temp;
    ZDOTA(nact) = alpha*ZDOTA(k);
    ZDOTA(k) = temp;
    LOOP(i,n) {
      temp = alpha*Z(i,nact) + beta*Z(i,k);
      Z(i,nact) = alpha*Z(i,k) - beta*Z(i,nact);
      Z(i,k) = temp;
    }
    IACT(nact) = IACT(k);
    IACT(k) = kk;
    temp = VMULTC(k);
    VMULTC(k) = VMULTC(nact);
    VMULTC(nact) = temp;
  }

  /* If stage one is in progress, then set SDIRN to the direction of the next
     change to the current vector of variables. */
  if (mcon > m) goto L_320;
  kk = IACT(nact);
  temp = zero;
  LOOP(i,n) {
    temp += SDIRN(i)*A(i,kk);
  }
  temp = (temp - one)/ZDOTA(nact);
  LOOP(i,n) {
    SDIRN(i) -= temp*Z(i,nact);
  }
  goto L_340;

  /* Delete the constraint that has the index IACT(ICON) from the active
     set. */
 L_260:
  if (icon < nact) {
    isave = IACT(icon);
    vsave = VMULTC(icon);
    k = icon;
    do {
      kp = k + 1;
      kk = IACT(kp);
      sp = zero;
      LOOP(i,n) {
        sp += Z(i,k)*A(i,kk);
      }
      temp = SQRT(sp*sp + ZDOTA(kp)*ZDOTA(kp));
      alpha = ZDOTA(kp)/temp;
      beta = sp/temp;
      ZDOTA(kp) = alpha*ZDOTA(k);
      ZDOTA(k) = temp;
      LOOP(i,n) {
        temp = alpha*Z(i,kp) + beta*Z(i,k);
        Z(i,kp) = alpha*Z(i,k) - beta*Z(i,kp);
        Z(i,k) = temp;
      }
      IACT(k) = kk;
      VMULTC(k) = VMULTC(kp);
      k = kp;
    } while (k < nact);
    IACT(k) = isave;
    VMULTC(k) = vsave;
  }
  --nact;

  /* If stage one is in progress, then set SDIRN to the direction of the next
     change to the current vector of variables. */
  if (mcon > m) goto L_320;
  temp = zero;
  LOOP(i,n) {
    temp += SDIRN(i)*Z(i,nact + 1);
  }
  LOOP(i,n) {
    SDIRN(i) -= temp*Z(i,nact + 1);
  }
  goto L_340;

  /* Pick the next search direction of stage two. */
 L_320:
  temp = one/ZDOTA(nact);
  LOOP(i,n) {
    SDIRN(i) = temp*Z(i,nact);
  }

  /* Calculate the step to the boundary of the trust region or take the step
     that reduces RESMAX to zero.  The two statements below that include the
     factor 1.0E-6 prevent some harmless underflows that occurred in a test
     calculation.  Further, we skip the step if it could be zero within a
     reasonable tolerance for computer rounding errors. */
 L_340:
  dd = rho*rho;
  sd = zero;
  ss = zero;
  LOOP(i,n) {
    if (ABS(DX(i)) >= tiny*rho) {
      dd -= DX(i)*DX(i);
    }
    sd += SDIRN(i)*DX(i);
    ss += SDIRN(i)*SDIRN(i);
  }
  if (dd <= zero) goto L_490;
  temp = SQRT(ss*dd);
  if (ABS(sd) >= tiny*temp) {
    temp = SQRT(ss*dd + sd*sd);
  }
  stpful = dd/(temp + sd);
  step = stpful;
  if (mcon == m) {
    acca = step + Op1*resmax;
    accb = step + Op2*resmax;
    if (step >= acca || acca >= accb) goto L_480;
    step = MIN(step, resmax);
  }

  /* Set DXNEW to the new variables if STEP is the steplength, and reduce
     RESMAX to the corresponding maximum residual if stage one is being done.
     Because DXNEW will be changed during the calculation of some Lagrange
     multipliers, it will be restored to the following value later. */
  LOOP(i,n) {
    DXNEW(i) = DX(i) + step*SDIRN(i);
  }
  if (mcon == m) {
    resold = resmax;
    resmax = zero;
    LOOP(k,nact) {
      kk = IACT(k);
      temp = B(kk);
      LOOP(i,n) {
        temp -= A(i,kk)*DXNEW(i);
      }
      resmax = MAX(resmax, temp);
    }
  }

  /* Set VMULTD to the VMULTC vector that would occur if DX became DXNEW.  A
     device is included to force VMULTD(K) = 0 if deviations from this value
     can be attributed to computer rounding errors.  First calculate the new
     Lagrange multipliers. */
  k = nact;
  for (;;) {
    zdotw = zero;
    zdwabs = zero;
    LOOP(i,n) {
      temp = Z(i,k)*DXNEW(i);
      zdotw += temp;
      zdwabs += ABS(temp);
    }
    acca = zdwabs + Op1*ABS(zdotw);
    accb = zdwabs + Op2*ABS(zdotw);
    if (zdwabs >= acca || acca >= accb) zdotw = zero;
    VMULTD(k) = zdotw/ZDOTA(k);
    if (k < 2) break;
    kk = IACT(k);
    LOOP(i,n) {
      DXNEW(i) -= VMULTD(k)*A(i,kk);
    }
    --k;
  }
  if (mcon > m && VMULTD(nact) < zero) VMULTD(nact) = zero;

  /* Complete VMULTC by finding the new constraint residuals. */
  LOOP(i,n) {
    DXNEW(i) = DX(i) + step*SDIRN(i);
  }
  if (mcon > nact) {
    kl = nact + 1;
    for (k = kl; k <= mcon; ++k) {
      kk = IACT(k);
      sum = resmax - B(kk);
      sumabs = resmax + ABS(B(kk));
      LOOP(i,n) {
        temp = A(i,kk)*DXNEW(i);
        sum += temp;
        sumabs += ABS(temp);
      }
      acca = sumabs + Op1*ABS(sum);
      accb = sumabs + Op2*ABS(sum);
      if (sumabs >= acca || acca >= accb) sum = zero;
      VMULTD(k) = sum;
    }
  }

  /* Calculate the fraction of the step from DX to DXNEW that will be taken. */
  ratio = one;
  icon = 0;
  LOOP(k,mcon) {
    if (VMULTD(k) < zero) {
      temp = VMULTC(k)/(VMULTC(k) - VMULTD(k));
      if (temp < ratio) {
        ratio = temp;
        icon = k;
      }
    }
  }

  /* Update DX, VMULTC and RESMAX. */
  temp = one - ratio;
  LOOP(i,n) {
    DX(i) = temp*DX(i) + ratio*DXNEW(i);
  }
  LOOP(k,mcon) {
    tempb = temp*VMULTC(k) + ratio*VMULTD(k);
    VMULTC(k) = MAX(zero, tempb);
  }
  if (mcon == m) resmax = resold + ratio*(resmax - resold);

  /* If the full step is not acceptable then begin another iteration.
     Otherwise switch to stage two or end the calculation. */
  if (icon > 0) goto L_70;
  if (step == stpful) goto L_500;
 L_480:
  mcon = m + 1;
  icon = mcon;
  IACT(mcon) = mcon;
  VMULTC(mcon) = zero;
  goto L_60;

  /* We employ any freedom that may be available to reduce the objective
     function before returning a DX whose length is less than RHO. */
 L_490:
  if (mcon == m) goto L_480;
  *ifull = 0;
 L_500:
  return;
}

/* Undefine macros mimicking FORTRAN indexing. */
#undef A
#undef B
#undef DX
#undef DXNEW
#undef IACT
#undef SDIRN
#undef VMULTC
#undef VMULTD
#undef Z
#undef ZDOTA

/*---------------------------------------------------------------------------*/

const char*
cobyla_reason(int status)
{
  switch (status) {
  case COBYLA_ITERATE:
    return "user requested to compute F(X) and C(X)";
  case COBYLA_SUCCESS:
    return "algorithm was successful";
  case COBYLA_ROUNDING_ERRORS:
    return "rounding errors are becoming damaging";
  case COBYLA_TOO_MANY_EVALUATIONS:
    return "MAXFUN limit has been reached";
  case COBYLA_BAD_ADDRESS:
    return "illegal NULL address";
  case COBYLA_CORRUPTED:
    return "unexpected parameter or corrupted workspace";
  default:
    return "unknown status";
  }
}

static void
print_calcfc(FILE* output, INTEGER n, INTEGER nfvals,
             REAL f, REAL maxcv, const REAL x[])
{
  INTEGER i;
  fprintf(output,
          "\n   NFVALS =%5d   F =%13.6E    MAXCV =%13.6E"
          "\n   X =%13.6E", (int)nfvals, f, maxcv, x[0]);
  for (i = 1; i < n; ++i) {
    fprintf(output, (i%5 == 0 ? "\n%19.6E" : "%15.6E"), x[i]);
  }
  fprintf(output, "\n");
}

/*---------------------------------------------------------------------------*/

#if defined(TESTING_REVCOM)
#  if defined(TESTING_FWRAP)
#    error only one of TESTING_REVCOM or TESTING_FWRAP can be defined
#  endif
#  define TESTING
static void
testing_revcom(INTEGER n, INTEGER m, REAL rhobeg, REAL rhoend,
               INTEGER iprint, INTEGER maxfun, REAL x[]);
#elif defined(TESTING_FWRAP)
#  define TESTING
#endif

#ifdef TESTING

/***********************************************************/
/* Main program of test problems in Report DAMTP 1992/NA5. */
/***********************************************************/

static int nprob = 1; /* Problem number. */
static cobyla_calcfc calcfc;

int
main(void)
{
  REAL r1, rhobeg, rhoend, temp, tempa, tempb, tempc, tempd;
  REAL x[10], xopt[10];
  INTEGER i, m, n, icase, maxfun, iprint;
#ifndef TESTING_REVCOM
  REAL w[3000];
  INTEGER iact[51];
#endif

  for (nprob = 1; nprob <= 10; ++nprob) {

#   define PRT(s) fprintf(stdout, "\n       " s "\n")

    switch(nprob) {

    case 1:
      /* Minimization of a simple quadratic function of two variables. */
      PRT("Output from test problem 1 (Simple quadratic)");
      n = 2;
      m = 0;
      xopt[0] = -1.0;
      xopt[1] = 0.0;
      break;

    case 2:
      /* Easy two dimensional minimization in unit circle. */
      PRT("Output from test problem 2 (2D unit circle calculation)");
      n = 2;
      m = 1;
      xopt[0] = SQRT(0.5);
      xopt[1] = -xopt[0];
      break;

    case 3:
      /* Easy three dimensional minimization in ellipsoid. */
      PRT("Output from test problem 3 (3D ellipsoid calculation)");
      n = 3;
      m = 1;
      xopt[0] = 1.0/SQRT(3.0);
      xopt[1] = 1.0/SQRT(6.0);
      xopt[2] = -0.33333333333333331;
      break;

    case 4:
      /* Weak version of Rosenbrock's problem. */
      PRT("Output from test problem 4 (Weak Rosenbrock)");
      n = 2;
      m = 0;
      xopt[0] = -1.0;
      xopt[1] = 1.0;
      break;

    case 5:
      /* Intermediate version of Rosenbrock's problem. */
      PRT("Output from test problem 5 (Intermediate Rosenbrock)");
      n = 2;
      m = 0;
      xopt[0] = -1.0;
      xopt[1] = 1.0;
      break;

    case 6:
      /* This problem is taken from Fletcher's book Practical Methods of
         Optimization and has the equation number (9.1.15). */
      PRT("Output from test problem 6 (Equation (9.1.15) in Fletcher)");
      n = 2;
      m = 2;
      xopt[0] = SQRT(0.5);
      xopt[1] = xopt[0];
      break;

    case 7:
      /* This problem is taken from Fletcher's book Practical Methods of
         Optimization and has the equation number (14.4.2). */
      PRT("Output from test problem 7 (Equation (14.4.2) in Fletcher)");
      n = 3;
      m = 3;
      xopt[0] = 0.0;
      xopt[1] = -3.0;
      xopt[2] = -3.0;
      break;

    case 8:
      /* This problem is taken from page 66 of Hock and Schittkowski's book
         Test Examples for Nonlinear Programming Codes. It is their test
         problem Number 43, and has the name Rosen-Suzuki. */
      PRT("Output from test problem 8 (Rosen-Suzuki)");
      n = 4;
      m = 3;
      xopt[0] = 0.0;
      xopt[1] = 1.0;
      xopt[2] = 2.0;
      xopt[3] = -1.0;
      break;

    case 9:
      /* This problem is taken from page 111 of Hock and Schittkowski's book
         Test Examples for Nonlinear Programming Codes. It is their test
         problem Number 100. */
      PRT("Output from test problem 9 (Hock and Schittkowski 100)");
      n = 7;
      m = 4;
      xopt[0] =  2.330499;
      xopt[1] =  1.951372;
      xopt[2] = -0.4775414;
      xopt[3] =  4.365726;
      xopt[4] = -0.624487;
      xopt[5] =  1.038131;
      xopt[6] =  1.594227;
      break;

    case 10:
    default:
      /* This problem is taken from page 415 of Luenberger's book Applied
         Nonlinear Programming. It is to maximize the area of a hexagon of unit
         diameter. */
      PRT("Output from test problem 10 (Hexagon area)");
      n = 9;
      m = 14;
      break;

    }

#  undef PRT

    for (icase = 1; icase <= 2; ++icase) {
      for (i = 1; i <= n; ++i) {
        x[i - 1] = 1.0;
      }
      rhobeg = 0.5;
      rhoend = ((icase == 2) ? 1e-4 : 0.001);
      iprint = 1;
      maxfun = 2000;
#if defined(TESTING_REVCOM)
      testing_revcom(n, m, rhobeg, rhoend, iprint, maxfun, x);
#elif defined(TESTING_FWRAP)
      cobyla_(&n, &m, x, &rhobeg, &rhoend, &iprint, &maxfun, w, iact);
#else
      cobyla(n, m, calcfc, NULL, x, rhobeg, rhoend, iprint, &maxfun, w, iact);
#endif
      if (nprob == 10) {
        tempa = x[0] + x[2] + x[4] + x[6];
        tempb = x[1] + x[3] + x[5] + x[7];
        tempc = 0.5/SQRT(tempa*tempa + tempb*tempb);
        tempd = tempc*SQRT(3.0);
        xopt[0] = tempd*tempa + tempc*tempb;
        xopt[1] = tempd*tempb - tempc*tempa;
        xopt[2] = tempd*tempa - tempc*tempb;
        xopt[3] = tempd*tempb + tempc*tempa;
        for (i = 1; i <= 4; ++i) {
          xopt[i + 3] = xopt[i - 1];
        }
      }
      temp = 0.0;
      for (i = 1; i <= n; ++i) {
        r1 = x[i - 1] - xopt[i - 1];
        temp += r1*r1;
      }
      fprintf(stdout, "\n     Least squares error in variables =%16.6E\n",
              (double)SQRT(temp));
    }
    fprintf(stdout, "  ----------------------------------------------"
            "--------------------\n");
  }
  return 0;
}

void
calcfc_(INTEGER* n, INTEGER* m, const REAL x[], REAL* f, REAL con[])
{
  *f = calcfc(*n, *m, x, con, NULL);
}

static REAL
calcfc(INTEGER n, INTEGER m, const REAL x[], REAL con[], void* data)
{
  REAL r1, r2, r3, r4, r5, r6, r7, fc;

#define POW2(x) ((x)*(x))
#define C(x)    FLT(x)
#define CON(a1) con[a1 - 1]
#define X(a1)   x[a1 - 1]


  /* Beware that order of operations may affect the result (whithin rounding
     errors).  I have tried to keep the same ordering as F2C which takes care
     of that, in particular when converting expressions involving powers. */

  switch (nprob) {

  case 1: /* Test problem 1 (Simple quadratic) */
    r1 = X(1) + C(1.0);
    r2 = X(2);
    fc = C(10.0)*(r1*r1) + (r2*r2);
    break;

  case 2: /* Test problem 2 (2D unit circle calculation) */
    fc = X(1)*X(2);
    CON(1) = C(1.0) - X(1)*X(1) - X(2)*X(2);
    break;

  case 3: /* Test problem 3 (3D ellipsoid calculation) */
    fc = X(1)*X(2)*X(3);
    CON(1) = C(1.0) - (X(1)*X(1)) - C(2.0)*(X(2)*X(2))
      - C(3.0)*(X(3)*X(3));
    break;

  case 4: /* Test problem 4 (Weak Rosenbrock) */
    r2 = X(1);
    r1 = r2*r2 - X(2);
    r3 = X(1) + C(1.0);
    fc = r1*r1 + r3*r3;
    break;

  case 5: /* Test problem 5 (Intermediate Rosenbrock) */
    r2 = X(1);
    r1 = r2*r2 - X(2);
    r3 = X(1) + C(1.0);
    fc = r1*r1*C(10.0) + r3*r3;
    break;

  case 6: /* Test problem 6 (Equation (9.1.15) in Fletcher's book) */
    fc = -X(1) - X(2);
    r1 = X(1);
    CON(1) = X(2) - r1*r1;
    r1 = X(1);
    r2 = X(2);
    CON(2) = C(1.0) - r1*r1 - r2*r2;
    break;

  case 7: /* Test problem 7 (Equation (14.4.2) in Fletcher's book) */
    fc = X(3);
    CON(1) = X(1)*C(5.0) - X(2) + X(3);
    r1 = X(1);
    r2 = X(2);
    CON(2) = X(3) - r1*r1 - r2*r2 - X(2)*C(4.0);
    CON(3) = X(3) - X(1)*C(5.0) - X(2);
    break;

  case 8: /* Test problem 8 (Rosen-Suzuki) */
    r1 = X(1);
    r2 = X(2);
    r3 = X(3);
    r4 = X(4);
    fc = r1*r1 + r2*r2 + r3*r3*C(2.0) + r4*r4 - X(1)*C(5.0) - X(2)*C(5.0)
      - X(3)*C(21.0) + X(4)*C(7.0);
    r1 = X(1);
    r2 = X(2);
    r3 = X(3);
    r4 = X(4);
    CON(1) = C(8.0) - r1*r1 - r2*r2 - r3*r3 - r4*r4 - X(1) + X(2) - X(3)
      + X(4);
    r1 = X(1);
    r2 = X(2);
    r3 = X(3);
    r4 = X(4);
    CON(2) = C(10.0) - r1*r1 - r2*r2*C(2.0) - r3*r3 - r4*r4*C(2.0) + X(1)
      + X(4);
    r1 = X(1);
    r2 = X(2);
    r3 = X(3);
    CON(3) = C(5.0) - r1*r1*C(2.0) - r2*r2 - r3*r3 - X(1)*C(2.0) + X(2) + X(4);
    break;

  case 9: /* Test problem 9 (Hock and Schittkowski 100) */
    r1 = X(1) - C(10.0);
    r2 = X(2) - C(12.0);
    r3 = X(3);
    r3 *= r3;
    r4 = X(4) - C(11.0);
    r5 = X(5);
    r5 *= r5;
    r6 = X(6);
    r7 = X(7);
    r7 *= r7;
    fc = r1*r1 + r2*r2*C(5.0) + r3*r3 + r4*r4*C(3.0) + r5*(r5*r5)*C(10.0)
      + r6*r6*C(7.0) + r7*r7 - X(6)*C(4.0)*X(7) - X(6)*C(10.0) - X(7)*C(8.0);
    r1 = X(1);
    r2 = X(2);
    r2 *= r2;
    r3 = X(4);
    CON(1) = C(127.0) - r1*r1*C(2.0) - r2*r2*C(3.0) - X(3) - r3*r3*C(4.0)
      - X(5)*C(5.0);
    r1 = X(3);
    CON(2) = C(282.0) - X(1)*C(7.0) - X(2)*C(3.0) - r1*r1*C(10.0)
      - X(4) + X(5);
    r1 = X(2);
    r2 = X(6);
    CON(3) = C(196.0) - X(1)*C(23.0) - r1*r1 - r2*r2*C(6.0) + X(7)*C(8.0);
    r1 = X(1);
    r2 = X(2);
    r3 = X(3);
    CON(4) = r1*r1*-C(4.0) - r2*r2 + X(1)*C(3.0)*X(2) - r3*r3*C(2.0)
      - X(6)*C(5.0) + X(7)*C(11.0);
    break;

  case 10: /* Test problem 10 (Hexagon area) */
  default:
    fc = -C(0.5)*(X(1)*X(4) - X(2)*X(3) + X(3)*X(9) - X(5)*X(9) + X(5)*X(8)
               - X(6)*X(7));
    r1 = X(3);
    r2 = X(4);
    CON(1) = C(1.0) - r1*r1 - r2*r2;
    r1 = X(9);
    CON(2) = C(1.0) - r1*r1;
    r1 = X(5);
    r2 = X(6);
    CON(3) = C(1.0) - r1*r1 - r2*r2;
    r1 = X(1);
    r2 = X(2) - X(9);
    CON(4) = C(1.0) - r1*r1 - r2*r2;
    r1 = X(1) - X(5);
    r2 = X(2) - X(6);
    CON(5) = C(1.0) - r1*r1 - r2*r2;
    r1 = X(1) - X(7);
    r2 = X(2) - X(8);
    CON(6) = C(1.0) - r1*r1 - r2*r2;
    r1 = X(3) - X(5);
    r2 = X(4) - X(6);
    CON(7) = C(1.0) - r1*r1 - r2*r2;
    r1 = X(3) - X(7);
    r2 = X(4) - X(8);
    CON(8) = C(1.0) - r1*r1 - r2*r2;
    r1 = X(7);
    r2 = X(8) - X(9);
    CON(9) = C(1.0) - r1*r1 - r2*r2;
    CON(10) = X(1)*X(4) - X(2)*X(3);
    CON(11) = X(3)*X(9);
    CON(12) = -X(5)*X(9);
    CON(13) = X(5)*X(8) - X(6)*X(7);
    CON(14) = X(9);
  }
#undef POW2
#undef C
#undef CON
#undef X

  return fc;
}

#ifdef TESTING_REVCOM
static void
testing_revcom(INTEGER n, INTEGER m, REAL rhobeg, REAL rhoend,
               INTEGER iprint, INTEGER maxfun, REAL x[])
{
  REAL f;
  REAL* c;
  cobyla_context_t* ctx;
  int status;
  const char* reason;

  if (m > 0) {
    c = (REAL*)malloc(m*sizeof(REAL));
    if (c == NULL) {
      goto enomem;
    }
  } else {
    c = NULL;
  }
  ctx = cobyla_create(n, m, rhobeg, rhoend, iprint, maxfun);
  if (ctx == NULL) {
    if (errno == ENOMEM) {
      goto enomem;
    } else {
      goto einval;
    }
  }
  status = cobyla_get_status(ctx);
  while (status == COBYLA_ITERATE) {
    f = calcfc(n, m, x, c, NULL);
    status = cobyla_iterate(ctx, f, x, c);
  }
  cobyla_delete(ctx);
#if 0
  if (status == COBYLA_SUCCESS) {
    return;
  }
  reason = cobyla_reason(status);
#else
  return;
#endif

 error:
  fprintf(stderr, "Something work occured in COBYLA: %s\n", reason);
  return;

 enomem:
   reason = "insufficient memory";
   goto error;

 einval:
   reason = "invalid parameters";
   goto error;
}
#endif /* TESTING_REVCOM */

#endif /* TESTING */

#endif /* _COBYLA_PART2 */

/*---------------------------------------------------------------------------*/
