/*
 * newuoa-test.i -
 *
 * Testing Yorick interface to NEWUOA.
 *
 * Copyright (c) 2015, Éric Thiébaut
 *
 * This file is part of the `Algorithms` (https://github.com/emmt/Algorithms)
 * project.  Read the accompanying `LICENSE` file for details.
 */

func newuoa_test(nil)
/* DOCUMENT newuoa_test;
     Test NEWUOA optimizer.
   SEE ALSO: newuoa_test_objfun, newuoa_create.
 */
{
  iprint = 2;
  maxfun = 5000;
  rhoend = 1e-6;
  for (n = 2; n <= 8; n += 2) {
    npt = 2*n + 1;
    x = indgen(n)/double(n + 1);
    rhobeg = x(1)*0.2;
    write, format="\n\n    Results with N =%2d and NPT =%3d\n", n, npt;

    nevals = 0;
    ctx = newuoa_create(n, npt, rhobeg, rhoend, iprint, maxfun);
    while (ctx.status == NEWUOA_ITERATE) {
      f = newuoa_test_objfun(x);
      newuoa_iterate, ctx, f, x;
      ++nevals;
    }
    //write, format="nevals: %d (%d)\n", nevals, ctx.nevals;
  }
}

func newuoa_test_objfun(x)
/* DOCUMENT f = newuoa_test_objfun(x);
     Return the value of the objective function for the Chebyquad test
     problem (Fletcher, 1965) for N = 2,4,6 and 8, with NPT = 2N+1.
   SEE ALSO: newuoa_test, newuoa_create.
 */
{
  n = numberof(x);
  np = n + 1;
  y = array(double, np, n);
  y(1,) = 1.0;
  y(2,) = x*2.0 - 1.0;
  for (i = 2; i <= n; ++i) {
    y(i+1,) = y(2,)*2.0*y(i,) - y(i-1,);
  }
  f = 0.0;
  iw = 1;
  for (i = 1; i <= np; ++i) {
    a = y(i,avg);
    if (iw > 0) {
      a += 1.0/double(i*i - 2*i);
    }
    iw = -iw;
    f += a*a;
  }
  return f;
}

if (batch()) {
  plug_dir, ".";
  include, "newuoa.i", 1;
  newuoa_test;
  quit;
}

/*---------------------------------------------------------------------------*/
