/*
 * tolmin.c --
 *
 * Implementation of Mike Powell's TOLMIN algorithm for minimizing a
 * differentiable function of many variables subject to linear and bound
 * constraints.  The algorithm is described in:
 *
 *    M.J.D. Powell, "A tolerant algorithm for linearly constrained
 *    optimization calculations", Math. Programming B, Vol. 45, pp. 547-566
 *    (1989).
 *
 * The present code is based on the original FORTRAN version written by Mike
 * Powell who released his code under the GNU Lesser General Public License.
 * His original code is available at CCPForge
 * <https://ccpforge.cse.rl.ac.uk/gf/project/powell/>
 *
 * ----------------------------------------------------------------------------
 *
 * Copyright (c) 1989, Mike J. D. Powell (FORTRAN version released under the
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
 *  - Pass objective function as argument.
 *  - Rename GETMIN as TOLMIN.
 *
 * Things to do:
 *  - Use 2 workspaces, one for integers, one for reals.
 *  - Define status constants and related messages.
 *  - Write FORTRAN wrapper with the same syntax as the original code.
 *  - Cleanup goto statements.
 */

#include <math.h>
#include <stdio.h>

#include "tolmin.h"

#define OUTPUT stdout


/* Fast macros to emulate FORTRAN intrinsics. */
#define sign(x,y)         _sign(typeof(x), x, y)
#define abs(x)            _abs(typeof(x), x)
#define pow2(x)           _pow2(typeof(x), x)
#define min(a,b)          _min(typeof((a) + (b)), a, b)
#define max(a,b)          _max(typeof((a) + (b)), a, b)
#define clamp(x,lo,hi)    _clamp(typeof((x) + (lo) + (hi)), x, lo, hi)
#define _abs(T,x)         ({ T _x = (x);  _x >= (T)0 ? _x : -_x; })
#define _pow2(T,x)        ({ T _x = (x);  _x*_x; })
#define _min(T,a,b)       ({ T _a = (a), _b = (b); _a < _b ? _a : _b; })
#define _max(T,a,b)       ({ T _a = (a), _b = (b); _a > _b ? _a : _b; })
#define _clamp(T,x,lo,hi) _min(T, _max(T, x, lo), hi)
#define _sign(T,x,y)      ({ T _x = (x); (_x < 0) != ((y) < 0) ? -_x : _x; })

/*
 * The following macros are defined to mimic FORTRAN indexing.  All recent
 * compilers will notice the constant offsets and make the necessary
 * adjustments.  So the readability is not to the detriment of efficiency.
 */
#define A(i1,i2)     a[(i2)*ia + (i1) - (ia+1)]
#define AM(i1)       am[(i1) - 1]
#define B(i1)        b[(i1) - 1]
#define BRES(i1)     bres[(i1) - 1]
#define CGRAD(i1)    cgrad[(i1) - 1]
#define D(i1)        d[(i1) - 1]
#define G(i1)        g[(i1) - 1]
#define GM(i1)       gm[(i1) - 1]
#define GMNEW(i1)    gmnew[(i1) - 1]
#define GOPT(i1)     gopt[(i1) - 1]
#define GS(i1)       gs[(i1) - 1]
#define IACT(i1)     iact[(i1) - 1]
#define PAR(i1)      par[(i1) - 1]
#define PARNEW(i1)   parnew[(i1) - 1]
#define PARW(i1)     parw[(i1) - 1]
#define RESKT(i1)    reskt[(i1) - 1]
#define RESKTW(i1)   resktw[(i1) - 1]
#define U(i1)        u[(i1) - 1]
#define W(i1)        w[(i1) - 1]
#define X(i1)        x[(i1) - 1]
#define XBIG(i1)     xbig[(i1) - 1]
#define XL(i1)       xl[(i1) - 1]
#define XS(i1)       xs[(i1) - 1]
#define XU(i1)       xu[(i1) - 1]
#define Z(i1)        z[(i1) - 1]
#define ZTC(i1)      ztc[(i1) - 1]
#define ZTG(i1)      ztg[(i1) - 1]
#define ZZDIAG(i1)   zzdiag[(i1) - 1]


static void
prtrvec(FILE* out, const char* prefix, const INTEGER ncols,
        const char* format, const INTEGER n, const REAL x[])
{
  INTEGER i;

  if (prefix != NULL) {
    fputs(prefix, out);
  }
  for (i = 0; i < n; ++i) {
    if (i > 0 && i%ncols == 0) {
      fputs("\n", out);
      if (prefix != NULL) {
        fputs(prefix, out);
      }
    }
    fprintf(out, format, (double)x[i]);
  }
  fputs("\n", out);
}


static void
prtivec(FILE* out, const char* prefix, const INTEGER ncols,
        const char* format, const INTEGER n, const INTEGER x[])
{
  INTEGER i;

  if (prefix != NULL) {
    fputs(prefix, out);
  }
  for (i = 0; i < n; ++i) {
    if (i > 0 && i%ncols == 0) {
      fputs("\n", out);
      if (prefix != NULL) {
        fputs(prefix, out);
      }
    }
    fprintf(out, format, (long)x[i]);
  }
  fputs("\n", out);
}


static void
addcon(const INTEGER n, const INTEGER m, const REAL a[], const INTEGER ia,
       INTEGER iact[], INTEGER* nact, REAL z[], REAL u[], const REAL relacc,
       const INTEGER indxbd, REAL ztc[], REAL cgrad[])
{
  INTEGER i, j, jp, np, iz, icon, inewbd, ipiv, iznbd;
  REAL sum, sumabs, temp, tempa, tempb, wcos, wsin, wpiv;

  iznbd = 0; /* FIXME: Fix uninitialized variables. */
  np = *nact + 1;
  icon = IACT(indxbd);
  IACT(indxbd) = IACT(np);
  IACT(np) = icon;
  if (icon > m) {
    /* Form ZTC when the new constraint is a bound. */
    inewbd = icon - m;
    if (inewbd <= n) {
      temp = -1.0;
    } else {
      inewbd -= n;
      temp = 1.0;
    }
    iznbd = inewbd*n - n;
    for (j = 1; j <= n; ++j) {
      ZTC(j) = temp*Z(iznbd+j);
    }
  } else {
    /* Else form ZTC for an ordinary constraint. */
    for (i = 1; i <= n; ++i) {
      CGRAD(i) = A(i,icon);
    }
    for (j = 1; j <= n; ++j) {
      ZTC(j) = 0.0;
      iz = j;
      for (i = 1; i <= n; ++i) {
        ZTC(j) += Z(iz)*CGRAD(i);
        iz += n;
      }
    }
  }

  /* Find any Givens rotations to apply to the last columns of Z. */
  j = n;
 L40:
  jp = j;
  --j;
  if (j > *nact) {
    if (ZTC(jp) == 0.0) {
      goto L40;
    }
    if (abs(ZTC(jp)) <= relacc*abs(ZTC(j))) {
      temp = abs(ZTC(j));
    } else if (abs(ZTC(j)) <= relacc*abs(ZTC(jp))) {
      temp = abs(ZTC(jp));
    } else {
      temp = abs(ZTC(jp))*sqrt(pow2(ZTC(j)/ZTC(jp)) + 1.0);
    }
    wcos = ZTC(j)/temp;
    wsin = ZTC(jp)/temp;
    ZTC(j) = temp;

    iz = j;
    if (icon > m) {
      /* Apply the rotation when the new constraint is a bound. */
      for (i = 1; i <= n; ++i) {
        temp = wcos*Z(iz+1) - wsin*Z(iz);
        Z(iz) = wcos*Z(iz) + wsin*Z(iz+1);
        Z(iz+1) = temp;
        iz += n;
      }
      Z(iznbd+jp) = 0.0;
    } else {
      /* Else apply the rotation for an ordinary constraint. */
      wpiv = 0.0;
      ipiv = 0;
      for (i = 1; i <= n; ++i) {
        tempa = wcos*Z(iz+1);
        tempb = wsin*Z(iz);
        temp = abs(CGRAD(i))*(abs(tempa) + abs(tempb));
        if (temp > wpiv) {
          wpiv = temp;
          ipiv = i;
        }
        Z(iz) = wcos*Z(iz) + wsin*Z(iz+1);
        Z(iz+1) = tempa - tempb;
        iz += n;
      }

      /* Ensure orthogonality of Z(.,JP) to CGRAD. */
      sum = 0.0;
      iz = jp;
      for (i = 1; i <= n; ++i) {
        sum += Z(iz)*CGRAD(i);
        iz += n;
      }
      if (sum != 0.0) {
        iz = ipiv*n - n + jp;
        Z(iz) -= sum/CGRAD(ipiv);
      }
    }
    goto L40;
  }

  /* Test for linear independence in the proposed new active set. */
  if (ZTC(np) == 0.0) {
    goto L90;
  }
  if (icon <= m) {
    sum = 0.0;
    sumabs = 0.0;
    iz = np;
    for (i = 1; i <= n; ++i) {
      temp = Z(iz)*CGRAD(i);
      sum += temp;
      sumabs += abs(temp);
      iz += n;
    }
    if (abs(sum) <= relacc*sumabs) {
      goto L90;
    }
  }

  /* Set the new diagonal element of U and return. */
  U(np) = 1.0/ZTC(np);
  *nact = np;
 L90:
  return;
}


static void
newcon(const INTEGER n, const INTEGER m, const REAL a[], const INTEGER ia,
       INTEGER iact[], INTEGER* nact, REAL z[], REAL u[], REAL d[],
       const REAL relacc, const INTEGER mdeg, REAL zzdiag[], REAL gmnew[],
       REAL cgrad[])
{
  INTEGER i, j, k, jm, np, iz, jmv, iadd, khigh;
  REAL cviol, cvmax, savabs, sumabs, savsum, sum, temp, sumd;

  /* FIXME: Fix uninitialized variables. */
  iadd = 0;
  savsum = 0.0;
  savabs = 0.0;

  /* Initialization. */
  np = *nact + 1;
  khigh = mdeg;
  iz = 0;
  for (i = 1; i <= n; ++i) {
    ZZDIAG(i) = 0.0;
    for (j = np; j <= n; ++j) {
      ZZDIAG(i) += pow2(Z(iz+j));
    }
    iz += n;
  }

  /* Calculate the scalar products of D with its constraints. */
 L30:
  cvmax = 0.0;
  for (k = np; k <= khigh; ++k) {
    j = IACT(k);
    if (j <= m) {
      sum = 0.0;
      sumabs = 0.0;
      sumd = 0.0;
      for (i = 1; i <= n; ++i) {
        temp = D(i)*A(i,j);
        sum += temp;
        sumabs += abs(temp);
        sumd += ZZDIAG(i)*pow2(A(i,j));
      }
    } else {
      jm = j - m;
      if (jm <= n) {
        sum = -D(jm);
      } else {
        jm -= n;
        sum = D(jm);
      }
      sumabs = abs(sum);
      sumd = ZZDIAG(jm);
    }

    /* Pick out the most violated constraint, or return if the
     *   violation is negligible. */
    if (sum > relacc*sumabs) {
      cviol = sum*sum/sumd;
      if (cviol > cvmax) {
        cvmax = cviol;
        iadd = k;
        savsum = sum;
        savabs = sumabs;
      }
    }
  }
  if (cvmax <= 0.0) {
    goto L140;
  }
  if (*nact == 0) {
    goto L120;
  }

  /* Set GMNEW to the gradient of the most violated constraint. */
  j = IACT(iadd);
  if (j <= m) {
    jmv = 0;
    for (i = 1; i <= n; ++i) {
      GMNEW(i) = A(i,j);
    }
  } else {
    jmv = j - m;
    for (i = 1; i <= n; ++i) {
      GMNEW(i) = 0.0;
    }
    if (jmv <= n) {
      GMNEW(jmv) = -1.0;
    } else {
      jmv -= n;
      GMNEW(jmv) = 1.0;
    }
  }

  /* Modify GMNEW for the next active constraint. */
  k = *nact;
 L80:
  temp = 0.0;
  iz = k;
  for (i = 1; i <= n; ++i) {
    temp += Z(iz)*GMNEW(i);
    iz += n;
  }
  temp *= U(k);
  j = IACT(k);
  if (j <= m) {
    for (i = 1; i <= n; ++i) {
      GMNEW(i) -= temp*A(i,j);
    }
  } else {
    jm = j - m;
    if (jm <= n) {
      GMNEW(jm) += temp;
    } else {
      GMNEW(jm-n) -= temp;
    }
  }

  /* Revise the values of SAVSUM and SAVABS. */
  sum = 0.0;
  sumabs = 0.0;
  for (i = 1; i <= n; ++i) {
    temp = D(i)*GMNEW(i);
    sum += temp;
    sumabs += abs(temp);
  }
  savsum = min(savsum, sum);
  savabs = max(savabs, sumabs);
  --k;
  if (k >= 1) {
    goto L80;
  }

  /* Add the new constraint to the active set if the constraint
   *   violation is still significant. */
  if (jmv > 0) {
    D(jmv) = 0.0;
  }
  if (savsum <= relacc*savabs) {
    goto L130;
  }
 L120:
  k = *nact;
  addcon(n, m, a, ia, iact, nact, z, u, relacc, iadd, gmnew, cgrad);
  if (*nact > k) {
    goto L140;
  }

  /* Seek another constraint violation. */
  iadd = np;
 L130:
  if (np < khigh) {
    k = IACT(khigh);
    IACT(khigh) = IACT(iadd);
    IACT(iadd) = k;
    --khigh;
    goto L30;
  }
 L140:
  return;
}


static void
delcon(const INTEGER n, const INTEGER m, const REAL a[], const INTEGER ia,
       INTEGER iact[], INTEGER* nact, REAL z[], REAL u[], const REAL relacc,
       const INTEGER idrop)
{
  INTEGER i, j, jp, nm, iz, ibd, icon, izbd, ipiv, isave;
  REAL rjjp, temp, tempb, ujp, sum, wcos, wsin, wpiv, denom, tempa;


  /* FIXME: Fix uninitialized variables. */
  ipiv = 0;
  izbd = 0;

  /* Function Body */
  nm = *nact - 1;
  if (idrop == *nact) {
    goto L60;
  }
  isave = IACT(idrop);

  /* Cycle through the constraint exchanges that are needed. */
  for (j = idrop; j <= nm; ++j) {
    jp = j + 1;
    icon = IACT(jp);
    IACT(j) = icon;

    /* Calculate the (J,JP) element of R. */
    if (icon <= m) {
      rjjp = 0.0;
      iz = j;
      for (i = 1; i <= n; ++i) {
        rjjp += Z(iz)*A(i,icon);
        iz += n;
      }
    } else {
      ibd = icon - m;
      if (ibd <= n) {
        izbd = ibd*n - n;
        rjjp = -Z(izbd+j);
      } else {
        ibd -= n;
        izbd = ibd*n - n;
        rjjp = Z(izbd+j);
      }
    }

    /* Calculate the parameters of the next rotation. */
    ujp = U(jp);
    temp = rjjp*ujp;
    denom = abs(temp);
    if (denom*relacc < 1.0) {
      denom = sqrt(denom*denom + 1.0);
    }
    wcos = temp/denom;
    wsin = 1.0/denom;

    /* Rotate Z when a bound constraint is promoted. */
    iz = j;
    if (icon > m) {
      for (i = 1; i <= n; ++i) {
        temp = wcos*Z(iz+1) - wsin*Z(iz);
        Z(iz) = wcos*Z(iz) + wsin*Z(iz+1);
        Z(iz+1) = temp;
        iz += n;
      }
      Z(izbd+jp) = 0.0;

      /* Rotate Z when an ordinary constraint is promoted. */
    } else {
      wpiv = 0.0;
      for (i = 1; i <= n; ++i) {
        tempa = wcos*Z(iz+1);
        tempb = wsin*Z(iz);
        temp = abs(A(i,icon))*(abs(tempa) + abs(tempb));
        if (temp > wpiv) {
          wpiv = temp;
          ipiv = i;
        }
        Z(iz) = wcos*Z(iz) + wsin*Z(iz+1);
        Z(iz+1) = tempa - tempb;
        iz += n;
      }

      /* Ensure orthogonality to promoted constraint. */
      sum = 0.0;
      iz = jp;
      for (i = 1; i <= n; ++i) {
        sum += Z(iz)*A(i,icon);
        iz += n;
      }
      if (sum != 0.0) {
        iz = ipiv*n - n + jp;
        Z(iz) -= sum/A(ipiv,icon);
      }
    }

    /* Set the new diagonal elements of U. */
    U(jp) = -denom*U(j);
    U(j) = ujp/denom;
  }

  /* Return. */
  IACT(*nact) = isave;
 L60:
  *nact = nm;
  return;
}


static void
sdirn(const INTEGER n, const INTEGER nact, const REAL z[], REAL d[],
      REAL ztg[], const REAL gm[], const REAL relacc, REAL* ddotgm)
{
  INTEGER i, j, np, iz;
  REAL sum, temp, sumabs;

  /* Function Body */
  *ddotgm = 0.0;
  if (nact >= n) {
    goto L60;
  }

  /* Premultiply GM by the transpose of Z. */
  np = nact + 1;
  for (j = np; j <= n; ++j) {
    sum = 0.0;
    sumabs = 0.0;
    iz = j;
    for (i = 1; i <= n; ++i) {
      temp = Z(iz)*GM(i);
      sum += temp;
      sumabs += abs(temp);
      iz += n;
    }
    if (abs(sum) <= relacc*sumabs) {
      sum = 0.0;
    }
    ZTG(j) = sum;
  }

  /* Form D by premultiplying ZTG by -Z. */
  iz = 0;
  for (i = 1; i <= n; ++i) {
    sum = 0.0;
    sumabs = 0.0;
    for (j = np; j <= n; ++j) {
      temp = Z(iz+j)*ZTG(j);
      sum -= temp;
      sumabs += abs(temp);
    }
    if (abs(sum) <= relacc*sumabs) {
      sum = 0.0;
    }
    D(i) = sum;
    iz += n;
  }

  /* Test that the search direction is downhill. */
  sumabs = 0.0;
  for (i = 1; i <= n; ++i) {
    temp = D(i)*GM(i);
    *ddotgm += temp;
    sumabs += abs(temp);
  }
  if (*ddotgm + relacc*sumabs >= 0.0) {
    *ddotgm = 0.0;
  }
 L60:
  return;
}

static void
sdegen(const INTEGER n, const INTEGER m, const REAL a[], const INTEGER ia,
       INTEGER iact[], INTEGER* nact, REAL par[], REAL z[], REAL u[],
       REAL d[], REAL ztg[], REAL gm[], const REAL relacc, REAL* ddotgm,
       const INTEGER meql, const INTEGER mdeg, REAL gmnew[], REAL parnew[],
       REAL cgrad[])
{
  INTEGER i, j, k, jm, mp, np, ku, iz, idrop, itest, n1;
  REAL dtest, sum, temp, theta, ratio;


  /* FIXME: Fix uninitialized variables. */
  idrop = 0;
  itest = 0;

  /* Function Body */
  mp = meql + 1;
  dtest = 0.0;

  /* Calculate the search direction and branch if it is not downhill. */
 L10:
  sdirn(n, *nact, z, d, ztg, gm, relacc, ddotgm);
  if (*ddotgm == 0.0) {
    goto L120;
  }

  /* Branch if there is no need to consider any degenerate constraints.
   *   The test gives termination if two consecutive additions to the
   *   active set fail to increase the predicted new value of F. */
  if (*nact == mdeg) {
    goto L120;
  }
  np = *nact + 1;
  sum = 0.0;
  for (j = np; j <= n; ++j) {
    sum += pow2(ZTG(j));
  }
  if (dtest > 0.0 && sum >= dtest) {
    if (itest == 1) {
      goto L120;
    }
    itest = 1;
  } else {
    dtest = sum;
    itest = 0;
  }

  /* Add a constraint to the active set if there are any significant
   *   violations of degenerate constraints. */
  k = *nact;
  newcon(n, m, a, ia, iact, nact, z, u, d, relacc, mdeg, gmnew, parnew, cgrad);
  if (*nact == k) {
    goto L120;
  }
  PAR(*nact) = 0.0;

  /* Calculate the new reduced gradient and Lagrange parameters. */
 L30:
  for (i = 1; i <= n; ++i) {
    GMNEW(i) = GM(i);
  }
  k = *nact;
 L50:
  temp = 0.0;
  iz = k;
  for (i = 1; i <= n; ++i) {
    temp += Z(iz)*GMNEW(i);
    iz += n;
  }
  temp *= U(k);
  PARNEW(k) = PAR(k) + temp;
  if (k == *nact) {
    PARNEW(k) = min(PARNEW(k), 0.0);
  }
  j = IACT(k);
  if (j <= m) {
    for (i = 1; i <= n; ++i) {
      GMNEW(i) -= temp*A(i,j);
    }
  } else {
    jm = j - m;
    if (jm <= n) {
      GMNEW(jm) += temp;
    } else {
      GMNEW(jm-n) -= temp;
    }
  }
  --k;
  if (k > meql) {
    goto L50;
  }

  /* Set RATIO for linear interpolation between PAR and PARNEW. */
  ratio = 0.0;
  if (mp < *nact) {
    ku = *nact - 1;
    for (k = mp; k <= ku; ++k) {
      if (PARNEW(k) > 0.0) {
        ratio = PARNEW(k)/(PARNEW(k) - PAR(k));
        idrop = k;
      }
    }
  }

  /* Apply the linear interpolation. */
  theta = 1.0 - ratio;
  n1 = *nact;
  for (k = mp; k <= n1; ++k) {
    PAR(k) = min(theta*PARNEW(k) + ratio*PAR(k) ,0.0);
  }
  for (i = 1; i <= n; ++i) {
    GM(i) = theta*GMNEW(i) + ratio*GM(i);
  }

  /* Drop a constraint if RATIO is positive. */
  if (ratio > 0.0) {
    delcon(n, m, a, ia, iact, nact, z, u, relacc, idrop);
    n1 = *nact;
    for (k = idrop; k <= n1; ++k) {
      PAR(k) = PAR(k+1);
    }
    goto L30;
  }

  /* Return if there is no freedom for a new search direction. */
  if (*nact < n) {
    goto L10;
  }
  *ddotgm = 0.0;
 L120:
  return;
}


static void
getd(const INTEGER n, const INTEGER m, const REAL a[], const INTEGER ia,
     INTEGER iact[], INTEGER* nact, REAL par[], const REAL g[], REAL z[],
     REAL u[], REAL d[], REAL ztg[], const REAL relacc, REAL* ddotg,
     const INTEGER meql, const INTEGER mdeg, REAL gm[], REAL gmnew[],
     REAL parnew[], REAL cgrad[])
{
  INTEGER i, j, k, jm, iz;
  REAL temp, ddotgm;

  /* Initialize GM and cycle backwards through the active set. */
 L10:
  for (i = 1; i <= n; ++i) {
    GM(i) = G(i);
  }
  k = *nact;
 L30:
  if (k > 0) {
    /* Set TEMP to the next multiplier, but reduce the active set if
     *   TEMP has an unacceptable sign. */
    temp = 0.0;
    iz = k;
    for (i = 1; i <= n; ++i) {
      temp += Z(iz)*GM(i);
      iz += n;
    }
    temp *= U(k);
    if (k > meql && temp > 0.0) {
      delcon(n, m, a, ia, iact, nact, z, u, relacc, k);
      goto L10;
    }

    /* Update GM using the multiplier that has just been calculated. */
    j = IACT(k);
    if (j <= m) {
      for (i = 1; i <= n; ++i) {
        GM(i) -= temp*A(i,j);
      }
    } else {
      jm = j - m;
      if (jm <= n) {
        GM(jm) += temp;
      } else {
        GM(jm-n) -= temp;
      }
    }
    PAR(k) = temp;
    --k;
    goto L30;
  }

  /* Calculate the search direction and DDOTG. */
  *ddotg = 0.0;
  if (*nact < n) {
    sdegen(n, m, a, ia, iact, nact, par, z, u, d, ztg, gm, relacc,
           &ddotgm, meql, mdeg, gmnew, parnew, cgrad);
    if (ddotgm < 0.0) {
      for (i = 1; i <= n; ++i) {
        *ddotg += D(i)*G(i);
      }
    }
  }
}


static void
stepbd(const INTEGER n, const INTEGER m, const REAL a[], const INTEGER ia,
       INTEGER iact[], REAL bres[], const REAL d[], REAL* stepcb, REAL* ddotg,
       const INTEGER mdeg, INTEGER* msat, const INTEGER mtot, INTEGER* indxbd)
{
  INTEGER i, j, k, jm, kl;
  int iflag;
  REAL sp, temp;


  /* Set steps to constraint boundaries and find the least positive one. */
  iflag = 0;
  *stepcb = 0.0;
  *indxbd = 0;
  k = mdeg;
 L10:
  ++k;
  if (k > mtot) {
    goto L40;
  }

  /* Form the scalar product of D with the current constraint normal. */
 L20:
  j = IACT(k);
  if (j <= m) {
    sp = 0.0;
    for (i = 1; i <= n; ++i) {
      sp += D(i)*A(i,j);
    }
  } else {
    jm = j - m;
    if (jm <= n) {
      sp = -D(jm);
    } else {
      sp = D(jm-n);
    }
  }

  /* The next branch is taken if label 20 was reached via label 50.0 */
  if (iflag == 1) {
    goto L60;
  }

  /* Set BRES(J) to indicate the status of the j-th constraint. */
  if (sp*BRES(j) <= 0.0) {
    BRES(j) = 0.0;
  } else {
    BRES(j) /= sp;
    if (*stepcb == 0.0 || BRES(j) < *stepcb) {
      *stepcb = BRES(j);
      *indxbd = k;
    }
  }
  goto L10;
 L40:

  /* Try to pass through the boundary of a violated constraint. */
 L50:
  if (*indxbd <= *msat) {
    goto L80;
  }
  iflag = 1;
  k = *indxbd;
  goto L20;
 L60:
  ++(*msat);
  IACT(*indxbd) = IACT(*msat);
  IACT(*msat) = j;
  BRES(j) = 0.0;
  *indxbd = *msat;
  *ddotg -= sp;
  if (*ddotg < 0.0 && *msat < mtot) {

    /* Seek the next constraint boundary along the search direction. */
    temp = 0.0;
    kl = mdeg + 1;
    for (k = kl; k <= mtot; ++k) {
      j = IACT(k);
      if (BRES(j) > 0.0) {
        if (temp == 0.0 || BRES(j) < temp) {
          temp = BRES(j);
          *indxbd = k;
        }
      }
    }
    if (temp > 0.0) {
      *stepcb = temp;
      goto L50;
    }
  }
 L80:
  return;
}


static void
conres(const INTEGER n, const INTEGER m, const REAL a[], const INTEGER ia,
       const REAL b[], const REAL xl[], const REAL xu[], const REAL x[],
       INTEGER iact[], INTEGER* nact, REAL par[], REAL g[], REAL z[], REAL u[],
       const REAL xbig[], REAL bres[], REAL d[], REAL ztg[], const REAL relacc,
       const REAL tol, REAL* stepcb, REAL* sumres, const INTEGER meql,
       INTEGER* msat, const INTEGER mtot, INTEGER* indxbd, REAL gm[],
       REAL gmnew[], REAL parnew[], REAL cgrad[])
{
  INTEGER i, j, k, kl, jm, idiff, mdeg, msatk;
  REAL ddotg, res, sum, resabs, temp;

  /* Function Body */
  idiff = mtot - *msat;

  /* Calculate and partition the residuals of the inactive constraints,
   *   and set the gradient vector when seeking feasibility. */
  if (idiff > 0) {
    for (i = 1; i <= n; ++i) {
      G(i) = 0.0;
    }
    *sumres = 0.0;
  }
  msatk = *msat;
  mdeg = *nact;
  *msat = *nact;
  kl = meql + 1;
  for (k = kl; k <= mtot; ++k) {
    j = IACT(k);

    /* Calculate the residual of the current constraint. */
    if (j <= m) {
      res = B(j);
      resabs = abs(B(j));
      for (i = 1; i <= n; ++i) {
        res -= X(i)*A(i,j);
        resabs += abs(XBIG(i)*A(i,j));
      }
    } else {
      jm = j - m;
      if (jm <= n) {
        res = X(jm) - XL(jm);
        resabs = abs(XBIG(jm)) + abs(XL(jm));
      } else {
        jm -= n;
        res = XU(jm) - X(jm);
        resabs = abs(XBIG(jm)) + abs(XU(jm));
      }
    }
    BRES(j) = res;

    /* Set TEMP to the relative residual. */
    temp = 0.0;
    if (resabs != 0.0) {
      temp = res/resabs;
    }
    if (k > msatk && temp < 0.0) {
      if (temp + relacc >= 0.0) {
        if (j <= m) {
          sum = abs(B(j));
          for (i = 1; i <= n; ++i) {
            sum += abs(X(i)*A(i,j));
          }
        } else {
          jm = j - m;
          if (jm <= n) {
            sum = abs(X(jm)) + abs(XL(jm));
          } else {
            sum = abs(X(jm-n)) + abs(XU(jm-n));
          }
        }
        if (abs(res) <= sum*relacc) {
          temp = 0.0;
        }
      }
    }

    /* Place the residual in the appropriate position. */
    if (k <= *nact) {
      goto L50;
    }
    if (k <= msatk || temp >= 0.0) {
      ++(*msat);
      if (*msat < k) {
        IACT(k) = IACT(*msat);
      }
      if (temp > tol) {
        IACT(*msat) = j;
      } else {
        ++mdeg;
        IACT(*msat) = IACT(mdeg);
        IACT(mdeg) = j;
      }

      /* Update the gradient and SUMRES if the constraint is violated when
       *   seeking feasibility. */
    } else {
      if (j <= m) {
        for (i = 1; i <= n; ++i) {
          G(i) += A(i,j);
        }
      } else {
        j -= m;
        if (j <= n) {
          G(j) += -1.0;
        } else {
          G(j-n) += 1.0;
        }
      }
      *sumres += abs(res);
    }
  L50:
    ;
  }

  /* Seek the next search direction unless CONRES was called from GETFES
   *   and feasibility has been achieved. */
  *stepcb = 0.0;
  if (idiff > 0 && *msat == mtot) {
    goto L60;
  }
  getd(n, m, a, ia, iact, nact, par, g, z, u, d, ztg, relacc, &ddotg,
       meql, mdeg, gm, gmnew, parnew, cgrad);

  /* Calculate the (bound on the) step-length due to the constraints. */
  if (ddotg < 0.0) {
    stepbd(n, m, a, ia, iact, bres, d, stepcb, &ddotg, mdeg, msat,
           mtot, indxbd);
  }
  if (idiff == 0) {
    *sumres = ddotg;
  }
 L60:
  return;
}


static void
ktvec(const INTEGER n, const INTEGER m, const REAL a[], const INTEGER ia,
      const INTEGER iact[], const INTEGER nact, REAL par[], const REAL g[],
      REAL reskt[], const REAL z[], const REAL u[], const REAL bres[],
      REAL* relaxf, const INTEGER meql, REAL* ssqkt, REAL parw[],
      REAL resktw[])
{
  INTEGER i, j, k, kk, jm, kl, iz;
  int icase;
  REAL ssqktw, temp;

  /* Calculate the Lagrange parameters and the residual vector. */
  /* FIXME: Fix uninitialized variables. */
  ssqktw = 0.0;

  /* Function Body */
  for (i = 1; i <= n; ++i) {
    RESKT(i) = G(i);
  }
  if (nact > 0) {
    icase = 0;
  L20:
    for (kk = 1; kk <= nact; ++kk) {
      k = nact + 1 - kk;
      j = IACT(k);
      temp = 0.0;
      iz = k;
      for (i = 1; i <= n; ++i) {
        temp += Z(iz)*RESKT(i);
        iz += n;
      }
      temp *= U(k);
      if (icase == 0) {
        PAR(k) = 0.0;
      }
      if (k <= meql || PAR(k) + temp < 0.0) {
        PAR(k) += temp;
      } else {
        temp = -PAR(k);
        PAR(k) = 0.0;
      }
      if (temp != 0.0) {
        if (j <= m) {
          for (i = 1; i <= n; ++i) {
            RESKT(i) -= temp*A(i,j);
          }
        } else {
          jm = j - m;
          if (jm <= n) {
            RESKT(jm) += temp;
          } else {
            RESKT(jm-n) -= temp;
          }
        }
      }
    }

    /* Calculate the sum of squares of the KT residual vector. */
    *ssqkt = 0.0;
    if (nact == n) {
      goto L130;
    }
    for (i = 1; i <= n; ++i) {
      *ssqkt += pow2(RESKT(i));
    }

    /* Apply iterative refinement to the residual vector. */
    if (icase == 0) {
      icase = 1;
      for (k = 1; k <= nact; ++k) {
        PARW(k) = PAR(k);
      }
      for (i = 1; i <= n; ++i) {
        RESKTW(i) = RESKT(i);
      }
      ssqktw = *ssqkt;
      goto L20;
    }

    /* Undo the iterative refinement if it does not reduce SSQKT. */
    if (ssqktw < *ssqkt) {
      for (k = 1; k <= nact; ++k) {
        PAR(k) = PARW(k);
      }
      for (i = 1; i <= n; ++i) {
        RESKT(i) = RESKTW(i);
      }
      *ssqkt = ssqktw;
    }

    /* Calculate SSQKT when there are no active constraints. */
  } else {
    *ssqkt = 0.0;
    for (i = 1; i <= n; ++i) {
      *ssqkt += pow2(G(i));
    }
  }

  /* Predict the reduction in F if one corrects any positive residuals
   *   of active inequality constraints. */
  *relaxf = 0.0;
  if (meql < nact) {
    kl = meql + 1;
    for (k = kl; k <= nact; ++k) {
      j = IACT(k);
      if (BRES(j) > 0.0) {
        *relaxf -= PAR(k)*BRES(j);
      }
    }
  }
 L130:
  return;
}


static void
lsrch(tolmin_objective fg, void* ctx, const INTEGER n, REAL x[], REAL g[],
      const REAL d[], REAL xs[], REAL gs[], const REAL relacc,
      const REAL stepcb, const REAL ddotg, REAL* f, REAL* step,
      INTEGER* nfvals, const INTEGER nfmax, REAL gopt[])
{
  INTEGER i, icount;
  REAL ddotgb, dgknot, relint, stphgh, stpmin, stplow, stpopt, ratio;
  REAL fhgh, temp, flow, fopt, fbase, dghgh, dgmid, sbase, dglow, dgopt;

  /* FIXME: Fix uninitialized variables. */
  dghgh = 0.0;
  fhgh = 0.0;

  /* Initialization. */
  relint = 0.9;
  icount = 0;
  ratio = -1.0;
  for (i = 1; i <= n; ++i) {
    XS(i) = X(i);
    GS(i) = G(i);
    GOPT(i) = G(i);
    if (D(i) != 0.0) {
      temp = abs(X(i)/D(i));
      if (ratio < 0.0 || temp < ratio) {
        ratio = temp;
      }
    }
  }
  *step = min(1.0,stepcb);

  /* The following number 1.0D-12 is independent of the working
   *   accuracy of the computer arithmetic. */
  stpmin = max(relacc*ratio, (*step)*1e-12);
  *step = max(stpmin, *step);
  sbase = 0.0;
  fbase = *f;
  ddotgb = ddotg;
  stplow = 0.0;
  flow = *f;
  dglow = ddotg;
  stphgh = 0.0;
  stpopt = 0.0;
  fopt = *f;
  dgopt = abs(ddotg);

  /* Calculate another function and gradient value. */
 L20:
  for (i = 1; i <= n; ++i) {
    X(i) = XS(i) + *step * D(i);
  }
  *f = fg(ctx, x, g);
  ++icount;
  dgmid = 0.0;
  for (i = 1; i <= n; ++i) {
    dgmid += D(i) * G(i);
  }
  if (*f <= fopt) {
    if (*f < fopt || abs(dgmid) < dgopt) {
      stpopt = *step;
      fopt = *f;
      for (i = 1; i <= n; ++i) {
        GOPT(i) = G(i);
      }
      dgopt = abs(dgmid);
    }
  }
  if (*nfvals + icount == nfmax) {
    goto L70;
  }

  /* Modify the bounds on the steplength or convergence. */
  if (*f >= fbase + (*step - sbase) * 0.1 * ddotgb) {
    if (stphgh > 0.0 || *f > fbase || dgmid > ddotg * 0.5) {
      stphgh = *step;
      fhgh = *f;
      dghgh = dgmid;
      goto L60;
    }
    sbase = *step;
    fbase = *f;
    ddotgb = dgmid;
  }
  if (dgmid >= ddotgb * 0.7) {
    goto L70;
  }
  stplow = *step;
  flow = *f;
  dglow = dgmid;
 L60:
  if (stphgh > 0.0 && stplow >= relint * stphgh) {
    goto L70;
  }

  /* Calculate the next step length or end the iterations. */
  if (stphgh == 0.0) {
    if (*step == stepcb) {
      goto L70;
    }
    temp = 10.0;
    if (dgmid > ddotg*0.9) {
      temp = ddotg/(ddotg - dgmid);
    }
    *step = min(temp*(*step), stepcb);
    goto L20;
  } else if (icount == 1 || stplow > 0.0) {
    dgknot = (fhgh - flow)*2.0/(stphgh - stplow) - (dglow + dghgh)*0.5;
    if (dgknot >= 0.0) {
      ratio = max(0.1, dglow*0.5/(dglow - dgknot));
    } else {
      ratio = (dghgh * 0.5 - dgknot)/(dghgh - dgknot);
    }
    *step = stplow + ratio*(stphgh - stplow);
    goto L20;
  } else {
    *step *= 0.1;
    if (*step >= stpmin) {
      goto L20;
    }
  }

  /* Return from subroutine. */
 L70:
  if (*step != stpopt) {
    *step = stpopt;
    *f = fopt;
    for (i = 1; i <= n; ++i) {
      X(i) = XS(i) + *step * D(i);
      G(i) = GOPT(i);
    }
  }
  *nfvals += icount;
  return;
}


static void
zbfgs(const INTEGER n, const REAL x[], const INTEGER nact, const REAL g[],
      REAL z[], REAL ztg[], REAL xs[], REAL gs[], REAL* zznorm)
{
  INTEGER i, k, km, kp, np, iz;
  REAL dd, dg, sum, temp, wcos, wsin;

  /* Test if there is sufficient convexity for the update. */
  dd = 0.0;
  dg = 0.0;
  temp = 0.0;
  for (i = 1; i <= n; ++i) {
    XS(i) = X(i) - XS(i);
    dd += pow2(XS(i));
    temp += GS(i)*XS(i);
    GS(i) = G(i) - GS(i);
    dg += GS(i)*XS(i);
  }
  if (dg < abs(temp)*0.1) {
    goto L90;
  }

  /* Transform the Z matrix. */
  k = n;
 L20:
  kp = k;
  --k;
  if (k > nact) {
    if (ZTG(kp) == 0.0) {
      goto L20;
    }
    temp = abs(ZTG(kp))*sqrt(pow2(ZTG(k)/ZTG(kp)) + 1.0);
    wcos = ZTG(k)/temp;
    wsin = ZTG(kp)/temp;
    ZTG(k) = temp;
    iz = k;
    for (i = 1; i <= n; ++i) {
      temp = wcos * Z(iz+1) - wsin * Z(iz);
      Z(iz) = wcos * Z(iz) + wsin * Z(iz+1);
      Z(iz+1) = temp;
      iz += n;
    }
    goto L20;
  }

  /* Update the value of ZZNORM. */
  if (*zznorm < 0.0) {
    *zznorm = dd/dg;
  } else {
    temp = sqrt(*zznorm * dd/dg);
    *zznorm = min(*zznorm, temp);
    *zznorm = max(*zznorm, temp*0.1);
  }

  /* Complete the updating of Z. */
  np = nact + 1;
  temp = sqrt(dg);
  iz = np;
  for (i = 1; i <= n; ++i) {
    Z(iz) = XS(i)/temp;
    iz += n;
  }
  if (np < n) {
    km = np + 1;
    for (k = km; k <= n; ++k) {
      temp = 0.0;
      iz = k;
      for (i = 1; i <= n; ++i) {
        temp += GS(i) * Z(iz);
        iz += n;
      }
      temp /= dg;
      sum = 0.0;
      iz = k;
      for (i = 1; i <= n; ++i) {
        Z(iz) -= temp * XS(i);
        sum += pow2(Z(iz));
        iz += n;
      }
      if (sum < *zznorm) {
        temp = sqrt(*zznorm/sum);
        iz = k;
        for (i = 1; i <= n; ++i) {
          Z(iz) = temp * Z(iz);
          iz += n;
        }
      }
    }
  }
 L90:
  return;
}


static void
minfun(tolmin_objective fg, void* ctx, const INTEGER n, const INTEGER m,
       const REAL a[], const INTEGER ia, const REAL b[], const REAL xl[],
       const REAL xu[], REAL x[], const REAL acc, INTEGER iact[],
       INTEGER* nact, REAL par[], const INTEGER iprint, int* info,
       REAL g[], REAL z[], REAL u[], REAL xbig[], const REAL relacc,
       REAL* zznorm, const REAL tol, const INTEGER meql, const INTEGER mtot,
       INTEGER* iterc, INTEGER* nfvals, const INTEGER nfmax, REAL reskt[],
       REAL bres[], REAL d[], REAL ztg[], REAL gm[], REAL xs[], REAL gs[])
{
  static REAL f; /* FIXME: static? */
  INTEGER i, k, indxbd, iterk, iterp, msat;
  REAL fprev, relaxf, ssqkt, step, ddotg, stepcb, sum, diff;

  /* FIXME: Fix uninitialized variables. */
  relaxf = 0.0;
  diff = 0.0;

  /* Initialize the minimization calculation. */
  msat = mtot;
  iterk = *iterc;
  if (*nfvals == 0 || *info == 1) {
    f = fg(ctx, x, g);
    ++(*nfvals);
  }
  fprev = abs(f + f + 1.0);
  iterp = -1;
  if (iprint != 0) {
    fprintf(OUTPUT, "\n     NEW VALUE OF TOL =%13.5E\n", (double)tol);
    iterp = (*iterc == 0) ? 0 : *iterc + abs(iprint);
  }

  /* Calculate the next search direction. */
 L10:
  conres(n, m, a, ia, b, xl, xu, x, iact, nact, par, g, z, u, xbig, bres,
         d, ztg, relacc, tol, &stepcb, &ddotg, meql, &msat, mtot, &indxbd,
         gm, reskt, xs, gs);

  /* Calculate the Kuhn Tucker residual vector. */
  ktvec(n, m, a, ia, iact, *nact, par, g, reskt, z, u, bres, &relaxf, meql,
        &ssqkt, xs, gs);

  /* Test for convergence. */
  if (ssqkt <= acc * acc) {
    *info = 1;
    goto L70;
  }
  if (ddotg >= 0.0) {
    *info = 2;
    goto L70;
  }

  /* Test for termination due to no decrease in F. */
  if (f >= fprev) {
    if (tol == relacc || *nact == 0) {
      if (diff > 0.0) {
        goto L20;
      }
    }
    *info = 3;
    goto L70;
  }
 L20:
  diff = fprev - f;
  fprev = f;

  /* Test that more calls of FGCALC are allowed. */
  if (*nfvals == nfmax) {
    *info = 8;
    goto L70;
  }

  /* Test whether to reduce TOL and to provide printing. */
  if (tol > relacc && *iterc > iterk && 0.1*relaxf >= max(diff, -0.5*ddotg)) {
    goto L70;
  }
  if (iterp == *iterc) {
    goto L80;
  }

  /* Calculate the step along the search direction. */
 L40:
  ++(*iterc);
  lsrch(fg, ctx, n, x, g, d, xs, gs, relacc, stepcb, ddotg, &f, &step,
        nfvals, nfmax, bres);
  if (step == 0.0) {
    *info = 3;
    sum = 0.0;
    for (i = 1; i <= n; ++i) {
      sum += abs(D(i) * GS(i));
    }
    if (ddotg + relacc * sum >= 0.0) {
      *info = 2;
    }
    goto L70;
  }

  /* Revise XBIG. */
  for (i = 1; i <= n; ++i) {
    XBIG(i) = max(XBIG(i), abs(X(i)));
  }

  /* Revise the second derivative approximation. */
  zbfgs(n, x, *nact, g, z, ztg, xs, gs, zznorm);

  /* Add a constraint to the active set if it restricts the step. */
  if (step == stepcb) {
    k = IACT(indxbd);
    if (k > m) {
      k -= m;
      if (k <= n) {
        X(k) = XL(k);
      } else {
        X(k-n) = XU(k-n);
      }
    }
    addcon(n, m, a, ia, iact, nact, z, u, relacc, indxbd, xs, gs);
  }
  goto L10;

  /* Printing from the subroutine. */
 L70:
  if (iprint == 0) {
    goto L90;
  }
  iterk = -1;
 L80:
  fprintf(OUTPUT, "\n     ITERS =%4ld     F.VALS =%4ld     F =%15.7E\n",
          (long)(*iterc), (long)(*nfvals), (double)f);
  fputs("  X =", OUTPUT); prtrvec(OUTPUT, NULL, 5, "%14.5E", n, x);
  fputs("  G =", OUTPUT); prtrvec(OUTPUT, NULL, 5, "%14.5E", n, g);
  if (iprint < 0) {
    if (*nact == 0) {
      fprintf(OUTPUT, "     NO ACTIVE CONSTRAINTS\n");
    } else {
      fputs(" IA =", OUTPUT); prtivec(OUTPUT, NULL, 5, "%5ld", *nact, iact);
      fputs(" LP =", OUTPUT); prtrvec(OUTPUT, NULL, 5, "%14.5E", *nact, par);
    }
    if (*nact == n) {
      fprintf(OUTPUT, "     KT RESIDUAL VECTOR IS ZERO\n");
    } else {
      fputs(" KT =", OUTPUT); prtrvec(OUTPUT, NULL, 5, "%14.5E", n, reskt);
    }
  }
  iterp = *iterc + abs(iprint);
  if (iterk >= 0) {
    goto L40;
  }
 L90:
  return;
}


static void
satact(const INTEGER n, const INTEGER m, const REAL a[], const INTEGER ia,
       const REAL b[], const REAL xl[], const REAL xu[], REAL x[],
       INTEGER iact[], INTEGER* nact, int* info, REAL z[], REAL u[],
       REAL xbig[], const REAL relacc, const REAL tol, const INTEGER meql)
{
  INTEGER i, j, k, jx, iz, idrop, n1;
  REAL res, temp, scale, tempa, resbig, resabs, savex;

  /* FIXME: Fix uninitialized variables. */
  savex = 0;
  jx = 0;

  /* Function Body */
  if (*nact == 0) {
    goto L50;
  }
  n1 = *nact;
  for (k = 1; k <= n1; ++k) {

    /* Calculate the next constraint residual. */
    j = IACT(k);
    if (j <= m) {
      res = B(j);
      resabs = abs(B(j));
      resbig = resabs;
      for (i = 1; i <= n; ++i) {
        tempa = A(i,j);
        temp = tempa*X(i);
        res -= temp;
        resabs += abs(temp);
        resbig += abs(tempa)*XBIG(i);
      }
    } else {
      jx = j - m;
      if (jx <= n) {
        res = X(jx) - XL(jx);
        resabs = abs(X(jx)) + abs(XL(jx));
        resbig = XBIG(jx) + abs(XL(jx));
        savex = XL(jx);
      } else {
        jx -= n;
        res = XU(jx) - X(jx);
        resabs = abs(X(jx)) + abs(XU(jx));
        resbig = XBIG(jx) + abs(XU(jx));
        savex = XU(jx);
      }
    }

    /* Shift X if necessary. */
    if (res != 0.0) {
      temp = res/resabs;
      if (k <= meql) {
        temp = -abs(temp);
      }
      if (tol == relacc || temp + relacc < 0.0) {
        *info = 1;
        scale = res*U(k);
        iz = k;
        for (i = 1; i <= n; ++i) {
          X(i) += scale*Z(iz);
          iz += n;
          XBIG(i) = max(XBIG(i), abs(X(i)));
        }
        if (j > m) {
          X(jx) = savex;
        }

        /* Else flag a constraint deletion if necessary. */
      } else if (res/resbig > tol) {
        IACT(k) = -IACT(k);
      }
    }
  }

  /* Delete any flagged constraints and then return. */
  idrop = *nact;
 L40:
  if (IACT(idrop) < 0) {
    IACT(idrop) = -IACT(idrop);
    delcon(n, m, a, ia, iact, nact, z, u, relacc, idrop);
  }
  --idrop;
  if (idrop > meql) {
    goto L40;
  }
 L50:
  return;
}


static void
adjtol(const INTEGER n, const INTEGER m, const REAL a[], const INTEGER ia,
       const REAL b[], const REAL xl[], const REAL xu[], const REAL x[],
       const INTEGER iact[], const INTEGER nact, REAL* xbig,
       const REAL relacc, REAL* tol, const INTEGER meql)
{
  INTEGER i, j, k, kl, jm;
  REAL res, viol, resabs;

  /* Set VIOL to the greatest relative constraint residual of the first
   *   NACT constraints. */
  /* Function Body */
  viol = 0.0;
  if (nact > meql) {
    kl = meql + 1;
    for (k = kl; k <= nact; ++k) {
      j = IACT(k);
      if (j <= m) {
        res = B(j);
        resabs = abs(B(j));
        for (i = 1; i <= n; ++i) {
          res -= A(i,j)*X(i);
          resabs += abs(A(i,j)*XBIG(i));
        }
      } else {
        jm = j - m;
        if (jm <= n) {
          res = X(jm) - XL(jm);
          resabs = XBIG(jm) + abs(XL(jm));
        } else {
          jm -= n;
          res = XU(jm) - X(jm);
          resabs = XBIG(jm) + abs(XU(jm));
        }
      }
      if (res > 0.0) {
        viol = max(viol, res/resabs);
      }
    }
  }

  /* Adjust TOL. */
  *tol = min(*tol,viol)*0.1;
  if (*tol <= relacc + relacc) {
    *tol = relacc;
    for (i = 1; i <= n; ++i) {
      XBIG(i) = abs(X(i));
    }
  }
  return;
}


static void
getfes(const INTEGER n, const INTEGER m, const REAL a[], const INTEGER ia,
       const REAL b[], const REAL xl[], const REAL xu[], REAL x[],
       INTEGER iact[], INTEGER* nact, REAL par[], int* info, REAL g[],
       REAL z[], REAL u[], REAL xbig[], const REAL relacc, REAL* tol,
       const INTEGER meql, INTEGER* msat, const INTEGER mtot, REAL bres[],
       REAL d[], REAL ztg[], REAL gm[], REAL gmnew[], REAL parnew[],
       REAL cgrad[])
{
  INTEGER i, msatk, itest, indxbd;
  REAL stepcb, sumres, sumrsk;

  /* Make the correction to X for the active constraints. */
  /* FIXME: Fix uninitialized variables. */
  itest = 0;

  /* Function Body */
  *info = 0;
 L10:
  satact(n, m, a, ia, b, xl, xu, x, iact, nact, info, z, u, xbig, relacc,
         *tol, meql);
  if (*info > 0) {
    *msat = *nact;
  }
  if (*msat == mtot) {
    goto L60;
  }

  /* Try to correct the infeasibility. */
 L20:
  msatk = *msat;
  sumrsk = 0.0;
 L30:
  conres(n, m, a, ia, b, xl, xu, x, iact, nact, par, g, z, u, xbig, bres, d,
         ztg, relacc, *tol, &stepcb, &sumres, meql, msat, mtot, &indxbd, gm,
         gmnew, parnew, cgrad);

  /* Include the new constraint in the active set. */
  if (stepcb > 0.0) {
    for (i = 1; i <= n; ++i) {
      X(i) += stepcb*D(i);
      XBIG(i) = max(XBIG(i), abs(X(i)));
    }
    addcon(n, m, a, ia, iact, nact, z, u, relacc, indxbd, gmnew, cgrad);
  }

  /* Test whether to continue the search for feasibility. */
  if (*msat < mtot) {
    if (stepcb == 0.0) {
      goto L50;
    }
    if (msatk < *msat) {
      goto L20;
    }
    if (sumrsk == 0.0 || sumres < sumrsk) {
      sumrsk = sumres;
      itest = 0;
    }
    ++itest;
    if (itest <= 2) {
      goto L30;
    }

    /* Reduce TOL if it may be too large to allow feasibility. */
  L50:
    if (*tol > relacc) {
      adjtol(n, m, a, ia, b, xl, xu, x, iact, *nact, xbig, relacc, tol, meql);
      goto L10;
    }
  }
 L60:
  return;
}

static REAL
getrelacc(void)
{
  volatile REAL tempa, tempb, ztpar, relacc;

  ztpar = 100.0;
  relacc = 1.0;
  do {
    relacc *= 0.5;
    tempa = ztpar + relacc*0.5;
    tempb = ztpar + relacc;
  } while (ztpar < tempa && tempa < tempb);
  return relacc;
}

static void
initzu(const INTEGER n, const INTEGER m, const REAL xl[], const REAL xu[],
       REAL x[], INTEGER iact[], INTEGER* meql, int* info, REAL z[],
       REAL u[], REAL xbig[])
{
  INTEGER i, j, nn, iz, jact;

  /* Seek bound inconsistencies and bound equality constraints. */
  *meql = 0;
  for (i = 1; i <= n; ++i) {
    if (XL(i) > XU(i)) {
      goto L50;
    }
    if (XL(i) == XU(i)) {
      ++(*meql);
    }
  }

  /* Initialize U, Z and XBIG. */
  jact = 0;
  nn = n*n;
  for (i = 1; i <= nn; ++i) {
    Z(i) = 0.0;
  }
  iz = 0;
  for (i = 1; i <= n; ++i) {
    if (XL(i) == XU(i)) {
      X(i) = XU(i);
      ++jact;
      U(jact) = 1.0;
      IACT(jact) = i + m + n;
      j = jact;
    } else {
      j = i + *meql - jact;
    }
    Z(iz+j) = 1.0;
    iz += n;
    XBIG(i) = abs(X(i));
  }
  *info = 1;
 L50:
  return;
}

static void
eqcons(const INTEGER n, const INTEGER m, const INTEGER meq, const REAL a[],
       const INTEGER ia, const REAL b[], const REAL xu[], INTEGER iact[],
       INTEGER* meql, int* info, REAL z[], REAL u[], const REAL relacc,
       REAL am[], REAL cgrad[])
{
  INTEGER i, j, k, jm, np, iz, keq;
  REAL rhs, sum, vmult;
  REAL sumabs;


  /* Try to add the next equality constraint to the active set. */
  /* Function Body */
  for (keq = 1; keq <= meq; ++keq) {
    if (*meql < n) {
      np = *meql + 1;
      IACT(np) = keq;
      addcon(n, m, a, ia, iact, meql, z, u, relacc, np, am, cgrad);
      if (*meql == np) {
        goto L50;
      }
    }

    /* If linear dependence occurs then find the multipliers of the
     *   dependence relation and apply them to the right hand sides. */
    sum = B(keq);
    sumabs = abs(B(keq));
    if (*meql > 0) {
      for (i = 1; i <= n; ++i) {
        AM(i) = A(i,keq);
      }
      k = *meql;
    L20:
      vmult = 0.0;
      iz = k;
      for (i = 1; i <= n; ++i) {
        vmult += Z(iz)*AM(i);
        iz += n;
      }
      vmult *= U(k);
      j = IACT(k);
      if (j <= m) {
        for (i = 1; i <= n; ++i) {
          AM(i) -= vmult*A(i,j);
        }
        rhs = B(j);
      } else {
        jm = j - m - n;
        AM(jm) -= vmult;
        rhs = XU(jm);
      }
      sum -= rhs*vmult;
      sumabs += abs(rhs*vmult);
      --k;
      if (k >= 1) {
        goto L20;
      }
    }

    /* Error return if the constraints are inconsistent. */
    if (abs(sum) > relacc*sumabs) {
      *info = 5;
      goto L60;
    }
  L50:
    ;
  }
 L60:
  return;
}


static int
minflc(tolmin_objective fg, void* ctx, const INTEGER n, const INTEGER m,
       const INTEGER meq, const REAL a[], const INTEGER ia, const REAL b[],
       const REAL xl[], const REAL xu[], REAL x[], const REAL acc,
       INTEGER iact[], INTEGER* nact, REAL par[], const INTEGER iprint,
       INTEGER nfmax, REAL g[], REAL z[], REAL u[], REAL xbig[], REAL reskt[],
       REAL bres[], REAL d[], REAL ztg[], REAL gm[], REAL xs[], REAL gs[])
{
  INTEGER i, k, mp, meql, msat, mtot, iterc, nfvals;
  REAL relacc, tol, zznorm;
  int info;

  /* Initialize ZZNORM, ITERC, NFVALS and NFMAX. */
  zznorm = -1.0;
  iterc = 0;
  nfvals = 0;

  /* Check the bounds on N, M and MEQ. */
  info = 4;
  if (max(max(1 - n, -m), meq*(meq - m)) > 0) {
    if (iprint != 0) {
      fprintf(OUTPUT, "\n     ERROR RETURN FROM TOLMIN BECAUSE %s\n",
              "A CONDITION ON N, M OR MEQ IS VIOLATED");
    }
    goto L40;
  }

  /* Initialize RELACC, Z, U and TOL. */
  relacc = getrelacc();
  initzu(n, m, xl, xu, x, iact, &meql, &info, z, u, xbig);
  tol = max(0.01, 10.0*relacc);
  if (info == 4) {
    if (iprint != 0) {
      fprintf(OUTPUT, "\n     ERROR RETURN FROM TOLMIN BECAUSE %s\n",
              "A LOWER BOUND EXCEEDS AN UPPER BOUND");
    }
    goto L40;
  }

  /* Add any equality constraints to the active set. */
  if (meq > 0) {
    eqcons(n, m, meq, a, ia, b, xu, iact, &meql, &info, z, u, relacc, xs, gs);
    if (info == 5) {
      if (iprint != 0) {
        fprintf(OUTPUT, "\n     ERROR RETURN FROM TOLMIN BECAUSE %s\n",
                "THE EQUALITY CONSTRAINTS ARE INCONSISTENT");
      }
      goto L40;
    }
  }
  *nact = meql;
  msat = meql;

  /* Add the bounds to the list of constraints. */
  mtot = *nact;
  for (i = 1; i <= n; ++i) {
    if (XL(i) < XU(i)) {
      mtot += 2;
      IACT(mtot-1) = m + i;
      IACT(mtot) = m + n + i;
    }
  }

  /* Try to satisfy the bound constraints. */
  getfes(n, m, a, ia, b, xl, xu, x, iact, nact, par, &info, g, z, u, xbig,
         relacc, &tol, meql, &msat, mtot, bres, d, ztg, gm, reskt, xs, gs);
  if (msat < mtot) {
    if (iprint != 0) {
      fprintf(OUTPUT, "\n     ERROR RETURN FROM TOLMIN BECAUSE %s\n",
              "THE EQUALITIES AND BOUNDS ARE INCONSISTENT");
    }
    info = 6;
    goto L40;
  }

  /* Add the ordinary inequalities to the list of constraints. */
  if (m > meq) {
    mp = meq + 1;
    for (k = mp; k <= m; ++k) {
      ++mtot;
      IACT(mtot) = k;
    }
  }

  /* Correct any constraint violations. */
 L30:
  getfes(n, m, a, ia, b, xl, xu, x, iact, nact, par, &info, g, z, u, xbig,
         relacc, &tol, meql, &msat, mtot, bres, d, ztg, gm, reskt, xs, gs);
  if (msat < mtot) {
    if (iprint != 0) {
      fprintf(OUTPUT, "\n     ERROR RETURN FROM TOLMIN BECAUSE %s\n",
              "THE CONSTRAINTS ARE INCONSISTENT");
    }
    info = 7;
    goto L40;
  } else if (meql == n) {
    if (iprint != 0) {
      fprintf(OUTPUT, "\n     %s\n",
              "TOLMIN FINDS THAT THE VARIABLES ARE DETERMINED "
              "BY THE EQUALITY CONSTRAINTS");
    }
    goto L40;
  }

  /* Minimize the objective function in the case when constraints are
   *   treated as degenerate if their residuals are less than TOL. */
  minfun(fg, ctx, n, m, a, ia, b, xl, xu, x, acc, iact, nact, par, iprint,
         &info, g, z, u, xbig, relacc, &zznorm, tol, meql, mtot, &iterc,
         &nfvals, nfmax, reskt, bres, d, ztg, gm, xs, gs);

  /* Reduce TOL if necessary. */
  if (tol > relacc && *nact > 0) {
    if (nfvals != nfmax) {
      adjtol(n, m, a, ia, b, xl, xu, x, iact, *nact, xbig, relacc, &tol, meql);
      goto L30;
    } else {
      info = 8;
    }
  }
  if (iprint != 0) {
    if (info == 1) {
      fprintf(OUTPUT, "\n     %s\n",
              "TOLMIN HAS ACHIEVED THE REQUIRED ACCURACY");
    } else if (info == 2) {
      fprintf(OUTPUT, "\n     %s\n",
              "TOLMIN CAN MAKE NO FURTHER PROGRESS BECAUSE "
              "OF ROUNDING ERRORS");
    } else if (info == 3) {
      fprintf(OUTPUT, "\n     %s\n",
              "TOLMIN CAN MAKE NO FURTHER PROGRESS BECAUSE "
              "F WILL NOT DECREASE ANY MORE");
    } else if (info == 8) {
      fprintf(OUTPUT, "\n     %s\n",
              "TOLMIN HAS REACHED THE GIVEN LIMIT ON "
              "THE NUMBER OF CALLS OF FGCALC");
    }
  }
 L40:
  return info;
}

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
 *    available in X(.) on the return from TOLMIN.
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
 *    multipliers PAR(K) K=1(1)NACT on the return from TOLMIN, these
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
 * NFMAX is an integer variable that should be nonnegative.  It should be
 *    zero if the caller does not want to impose an upper bound on the
 *    number of evaluations of the objective function and its gradient;
 *    otherwise it should be set to the value of this bound.
 * W(.) is a working space array of real variables that must be provided
 *    by the user.  Its length must be at least (M+11*N+N**2).  On exit
 *    from the package one can find the final components of GRAD(F) and
 *    RESKT(.) in W(1),...,W(N) and W(N+1),...,W(2*N) respectively.
 *
 * The returned value is an integer, say INFO, with one of the following
 *    integer values to indicate the reason for leaving the optimization
 *    package:
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
 *
 * Note 1.   The variables N, M, MEQ, IA, ACC and IPRINT and the elements
 *    of the arrays A(,.,), B(.), XL(.) and XU(.) are not altered by the
 *    optimization procedure.  Their values, the value of INFO and the
 *    initial components of X(.) must be set on entry to TOLMIN.
 *
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
int
tolmin(tolmin_objective fg, void* ctx, const INTEGER n, const INTEGER m,
       const INTEGER meq, const REAL a[], const INTEGER ia, const REAL b[],
       const REAL xl[], const REAL xu[], REAL x[], const REAL acc,
       INTEGER iact[], INTEGER* nact, REAL par[], const INTEGER iprint,
       INTEGER nfmax, REAL w[])
{
  INTEGER id, ig, iu, iz, igm, igs, ixs, iztg, ixbig, ibres;
  INTEGER ireskt;

  /* Partition the workspace array. */
  ig = 1;
  ireskt = ig + n;
  iz = ireskt + n;
  iu = iz + n*n;
  ixbig = iu + n;
  ibres = ixbig + n;
  id = ibres + m + n + n;
  iztg = id + n;
  igm = iztg + n;
  ixs = igm + n;
  igs = ixs + n;

  /* Call the optimization package. */
  return minflc(fg, ctx, n, m, meq, a, ia, b, xl, xu, x, acc, iact, nact, par,
                iprint, nfmax, &W(ig), &W(iz), &W(iu), &W(ixbig), &W(ireskt),
                &W(ibres), &W(id), &W(iztg), &W(igm), &W(ixs), &W(igs));
}
