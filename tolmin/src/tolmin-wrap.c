/*
 * tolmin-wrap.c --
 *
 * Wrapper for calling the C version of Mike Powell's TOLMIN from FORTRAN code.
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

extern void fgcalc_(const INTEGER *n, const REAL *x, REAL *f, REAL *g);

static REAL
wfg(void* ctx, const REAL x[], REAL g[])
{
  INTEGER n = *(INTEGER*)ctx;
  REAL f;

  fgcalc_(&n, x, &f, g);
  return f;
}

extern void
getmin_(const INTEGER* n, const INTEGER* m, const INTEGER* meq,
        const REAL a[], const INTEGER* ia, const REAL b[],
        const REAL xl[], const REAL xu[], REAL x[], const REAL* acc,
        INTEGER iact[], INTEGER* nact, REAL par[],
        const INTEGER* iprint, INTEGER* info, REAL w[])
{
  INTEGER number = *n;
  int status = tolmin(wfg, &number, *n, *m, *meq, a, *ia, b, xl, xu,
                      x, *acc, iact, nact, par, *iprint, *info, w);
  *info = status;
}
