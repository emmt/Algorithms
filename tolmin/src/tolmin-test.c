/*
 * tolmin-test.c --
 *
 * Test Mike Powell's TOLMIN algorithm on the pentagon problem.
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
 */

#include <math.h>
#include <stdio.h>

#include "tolmin.h"

#define OUTPUT stdout


/*
 * The following macros are defined to mimic FORTRAN indexing.  All recent
 * compilers will notice the constant offsets and make the necessary
 * adjustments.  So the readability is not to the detriment of efficiency.
 */
#define A(i1,i2)     a[(i2)*ia + (i1) - (ia+1)]
#define B(i1)        b[(i1) - 1]
#define G(i1)        g[(i1) - 1]
#define X(i1)        x[(i1) - 1]
#define XL(i1)       xl[(i1) - 1]
#define XU(i1)       xu[(i1) - 1]


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


static REAL pow2(REAL x)
{
  return x*x;
}


static REAL pow8(REAL x)
{
  REAL x2 = x*x;
  REAL x4 = x2*x2;
  return x4*x4;
}


static REAL pow9(REAL x)
{
  REAL x2 = x*x;
  REAL x4 = x2*x2;
  return x4*x4*x;
}


static REAL
fgcalc(void* ctx, const REAL x[], REAL g[])
{
  REAL wa, wb, wc, f;

  /* Calculate the objective function and its gradient. */
  wa = pow2(X(1) - X(3)) + pow2(X(2) - X(4));
  wb = pow2(X(3) - X(5)) + pow2(X(4) - X(6));
  wc = pow2(X(5) - X(1)) + pow2(X(6) - X(2));
  f = 1.0/pow8(wa) + 1.0/pow8(wb) + 1.0/pow8(wc);
  G(1) = ((X(3) - X(1))/pow9(wa) + (X(5) - X(1))/pow9(wc))*16.0;
  G(2) = ((X(4) - X(2))/pow9(wa) + (X(6) - X(2))/pow9(wc))*16.0;
  G(3) = ((X(5) - X(3))/pow9(wb) + (X(1) - X(3))/pow9(wa))*16.0;
  G(4) = ((X(6) - X(4))/pow9(wb) + (X(2) - X(4))/pow9(wa))*16.0;
  G(5) = ((X(1) - X(5))/pow9(wc) + (X(3) - X(5))/pow9(wb))*16.0;
  G(6) = ((X(2) - X(6))/pow9(wc) + (X(4) - X(6))/pow9(wb))*16.0;
  return f;
}


/* The pentagon problem. */
int main(int argc, char* argv[])
{
  const INTEGER ia = 10;
  const INTEGER n = 6;
  const INTEGER maxm = 15; /* max. numb. of constraints */

  /* Local variables */
  INTEGER i, j, k, m, nact, info, icase, iprint, meq;
  INTEGER iact[27];
  REAL a[ia*maxm], b[maxm], f;
  REAL xl[n], xu[n], acc, cth, par[20];
  REAL sth;
  REAL theta;
  REAL w[1000], x[n];

  /* The two values of ICASE provide two different values of ACC, the latter
   *   accuracy being so demanding that a return with INFO=2 occurs. The final
   *   values of the objective function in the two cases agree well and
   *   constraint violations are negligible, considering the differences
   *   between the final values of the variables. */

  iprint = -1;
  for (icase = 1; icase <= 2; ++icase) {
    acc = 1e-6;
    if (icase == 2) {
      acc = 1e-14;
    }

    /* Set the components of XL, XU and X. */
    for (i = 1; i <= n; ++i) {
      XL(i) = -1e6;
      XU(i) =  1e6;
      X(i) = (REAL)(i - 3)*0.5;
    }
    X(2) =  0.0;
    X(4) = -1.0;
    X(6) =  1.0;

    /* Set the constraints. */
    m = 0;
    meq = 0;
    for (k = 1; k <= 5; ++k) {
      theta = (REAL)(k - 1)*0.4*M_PI;
      cth = cos(theta);
      sth = sin(theta);
      for (j = 2; j <= n; j += 2) {
        ++m;
        for (i = 1; i <= n; ++i) {
          A(i,m) = 0.0;
        }
        A(j-1,m) = cth;
        A(j,m) = sth;
        B(m) = 1.0;
      }
    }

    /* Call the optimization package. */
    info = 0;
    fprintf(OUTPUT, "\n\n     CALL OF TOLMIN WITH  "
            "ACC =%11.4E  AND  IPRINT =%3ld\n",
            (double)acc, (long)iprint);
    tolmin(fgcalc, NULL, n, m, meq, a, ia, b, xl, xu, x, acc, iact,
           &nact, par, iprint, &info, w);
    fprintf(OUTPUT, "\n     RETURN FROM TOLMIN WITH INFO =%2d\n", (int)info);
    f = fgcalc(NULL, x, w);
    fprintf(OUTPUT, "\n     FINAL VALUE OF OBJECTIVE FUNCTION =%20.12E\n",
            (double)f);
    fprintf(OUTPUT, "\n     FINAL COMPONENTS OF X =\n\n");
    prtrvec(OUTPUT, "    ", 3, "%20.12E", n, x);
    for (k = 1; k <= m; ++k) {
      for (i = 1; i <= n; ++i) {
        B(k) -= A(i,k)*X(i);
      }
    }
    fprintf(OUTPUT, "\n     FINAL CONSTRAINT RESIDUALS =\n\n");
    prtrvec(OUTPUT, "   ", 6, "%12.4E", m, b);
  }
  return 0;
}
