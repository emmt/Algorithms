/*
 * cobyla.i -
 *
 * Yorick interface to COBYLA.
 *
 * Copyright (c) 2015, Éric Thiébaut
 *
 * This file is part of the `Algorithms` (https://github.com/emmt/Algorithms)
 * project.  Read the accompanying `LICENSE` file for details.
 */

if (is_func(plug_in)) plug_in, "ycobyla";

local cobyla_minimize, cobyla_maximize;
/* DOCUMENT xmin = cobyla_minimize(f, c, x0, rhobeg, rhoend);
         or  obj = cobyla_minimize(f, c, x0, rhobeg, rhoend, all=1);
         or xmax = cobyla_maximize(f, c, x0, rhobeg, rhoend);
         or  obj = cobyla_maximize(f, c, x0, rhobeg, rhoend, all=1);

     Minimize or maximize the multi-variate function F under the inequality
     constraints implemented by C and starting at the initial point X0.  RHOBEG
     and RHOEND are the initial and final values of the trust region radius (0
     < RHOEND <= RHOBEG).

     The objective is to find X which solves one of the problems:

          min f(x)    s.t.   c(x) <= 0   (cobyla_minimize)
          max f(x)    s.t.   c(x) <= 0   (cobyla_maximize)

     Arguments F and C are user defined functions which take a single argument,
     the variables X, and return the function value and the constraints
     respectively.  If there are no contraints, C can be empty.

     Note that the proper scaling of the variables is important for the success
     of the algorithm.  RHOBEG should be set to the typical size of the region
     to explorate and RHOEND should be set to the typical precision.

     Keyword NPT sets the number of of interpolation conditions.  Its default
     value is equal to 2*N+1 (the recommended value).

     Keyword MAXFUN sets the maximum number of function evaluations.  Its
     default value is equal to 30*N.

     If keyword ALL is true, the result is a structured object.  For
     `cobyla_minimize`, the members of the returned object are:

        obj.fmin   = minimal function value found
        obj.cmin   = constraints at the minimum
        obj.xmin   = corresponding parameters
        obj.nevals = number of function evaluations
        obj.status = status of the algorithm upon return
        obj.rho    = final radius of the trust region

     For `cobyla_maximize`, the two first members are:

        obj.fmax   = minimal function value found
        obj.cmax   = constraints at the maximum
        obj.xmax   = corresponding parameters

     Keyword ERR sets the behavior in case of abnormal termination.  If ERR=0,
     anything but a success throws an error (this is also the default
     behavior); if ERR > 0, non-fatal errors are reported by a warning message;
     if ERR < 0, non-fatal errors are silently ignored.

     Keyword VERB set the verbosity level.


   SEE ALSO: cobyla_create, cobyla_error.
 */

func cobyla_minimize(f, c, x0, rhobeg, rhoend, npt=, maxfun=, all=, verb=, err=)
{
  local cx, fx, fmin, cmin, xmin;
  x = double(unref(x0));
  n = numberof(x);
  constrained = (! is_void(c));
  if (constrained) cx = c(x);
  ctx = cobyla_create(n, numberof(cx), rhobeg, rhoend,
                      (is_void(verb) ? 0 : verb),
                      (is_void(maxfun) ? 50*n : maxfun));
  start = 1n;
  do {
    if (constrained && ! start) cx = c(x);
    fx = f(x);
    if (start || fmin > fx) {
      fmin = fx;
      cmin = cx;
      xmin = x;
      start = 0n;
    }
    status = cobyla_iterate(ctx, fx, x, cx);
  } while (status == COBYLA_ITERATE);
  if (status != COBYLA_SUCCESS) cobyla_error, status, err;
  return (all ? save(fmin, cmin, xmin, nevals = ctx.nevals,
                     status, rho = ctx.rho) : x);
}

func cobyla_maximize(f, c, x0, rhobeg, rhoend, npt=, maxfun=, all=, verb=, err=)
{
  local cx, fx, fmax, cmax, xmax;
  x = double(unref(x0));
  n = numberof(x);
  constrained = (! is_void(c));
  if (constrained) cx = c(x);
  ctx = cobyla_create(n, numberof(cx), rhobeg, rhoend,
                      (is_void(verb) ? 0 : verb),
                      (is_void(maxfun) ? 50*n : maxfun));
  start = 1n;
  do {
    if (constrained && ! start) cx = c(x);
    fx = f(x);
    if (start || fmax < fx) {
      fmax = fx;
      cmax = cx;
      xmax = x;
      start = 0n;
    }
    status = cobyla_iterate(ctx, -fx, x, cx);
  } while (status == COBYLA_ITERATE);
  if (status != COBYLA_SUCCESS) cobyla_error, status, err;
  return (all ? save(fmax, cmax, xmax, nevals = ctx.nevals,
                     status, rho = ctx.rho) : x);
}

func cobyla_error(status, errmode)
/* DOCUMENT cobyla_error, status;
         or cobyla_error, status, errmode;

     Report an error in COBYLA according to the value of STATUS.  Nothing is
     done if STATUS is COBYLA_SUCCESS; otherwise, the optional argument ERRMODE
     determines the behavior.  If ERRMODE = 0, the routine throws an error
     (this is also the default behavior); if ERRMODE > 0, non-fatal errors are
     reported by a warning message; if ERRMODE < 0, non-fatal errors are
     silently ignored.

   SEE ALSO: cobyla_reason, error.
 */
{
  if (status != COBYLA_SUCCESS) {
    if (errmode && (status == COBYLA_ROUNDING_ERRORS ||
                    status == COBYLA_TOO_MANY_EVALUATIONS ||
                    status == COBYLA_STEP_FAILED)) {
      if (errmode > 0) {
        write, format="WARNING: Something wrong occured in COBYLA: %s",
          cobyla_reason(status);
      }
    } else {
      error, swrite(format="Something wrong occured in COBYLA: %s",
                    cobyla_reason(status));
    }
  }
}
errs2caller, cobyla_error;

local COBYLA_SUCCESS, COBYLA_ITERATE, COBYLA_ROUNDING_ERRORS;
local COBYLA_TOO_MANY_EVALUATIONS, COBYLA_BAD_ADDRESS;
local COBYLA_CORRUPTED;
extern cobyla_create;
extern cobyla_iterate;
/* DOCUMENT ctx = cobyla_create(n, m, rhobeg, rhoend, iprint, maxfun);
         or status = cobyla_iterate(ctx, f, x, c);

     The function `cobyla_create` makes a new instance for Mike Powell's COBYLA
     algorithm for minimizing a function of a few variables.  The method is
     "derivatives free" (only the function values are needed) and accounts for
     inequality constraints on the variables.  N is the number of variables, M
     is the number of constraints, RHOBEG and RHOEND are the initial and final
     size of the trust region, IPRINT control the verbosity of the method (see
     below) and MAXFUN is the maximum number of function evaluations.

     The objective is to find X which solves the problem:

          min f(x)    s.t.   c(x) <= 0

     Given a penalty parameter SIGMA > 0, COBYLA assumes that a change to X is
     an improvement if it reduces the merit function:

          f(x) + sigma*max(0.0,-min(c(x))),

     where c(x) are the constraint functions that should become nonnegative
     eventually, at least to the precision of RHOEND.  In the printed output
     the displayed term that is multiplied by SIGMA is called MAXCV, which
     stands for 'MAXimum Constraint Violation'.

     IPRINT should be set to 0, 1, 2 or 3, which controls the amount of
     printing during the calculation.  Specifically, there is no output if
     IPRINT=0 and there is output only at the end of the calculation if
     IPRINT=1.  Otherwise each new value of RHO and SIGMA is printed.  Further,
     the vector of variables and some function information are given either
     when RHO is reduced or when each new value of F(X) is computed in the
     cases IPRINT=2 or IPRINT=3 respectively.

     The function `cobyla_iterate` performs an iteration of the algorithm given
     F the function value at X the current variables and C the constraints.

     Typical usage is:

     >   ctx = cobyla_create(n, m, rhobeg, rhoend, iprint, maxfun);
     >   x = ...; // initial solution
     >   while (ctx.status == COBYLA_ITERATE) {
     >     f = ...; // compute function value at X
     >     c = ...; // compute constraints at X, if any
     >     cobyla_iterate, ctx, f, x, c;
     >   }
     >   if (ctx.status != COBYLA_SUCCESS) {
     >     error, swrite(format="Something wrong occured in COBYLA: %s",
     >                   ctx.reason);

     If there are no constraints (M = 0), argument C must be [] (void) or can
     be omitted.

   REFERENCES
     The algorithm is described in:

         M.J.D. Powell, "A direct search optimization method that models the
         objective and constraint functions by linear interpolation," in
         Advances in Optimization and Numerical Analysis Mathematics and Its
         Applications, vol.  275 (eds.  Susana Gomez and Jean-Pierre Hennart),
         Kluwer Academic Publishers, pp. 51-67 (1994).

   SEE ALSO: cobyla_restart, cobyla_reason.
 */

extern cobyla_restart;
/* DOCUMENT cobyla_restart(ctx);
     Restart COBYLA algorithm using the same parameters.  The return value is
     the new status of the algorithm, see `cobyla_get_status` for details.

   SEE ALSO: cobyla_create.
 */

extern cobyla_reason;
/* DOCUMENT cobyla_reason(status);
     Get a textual explanation of the status returned by COBYLA.
   SEE ALSO: cobyla_create.
 */

extern cobyla_init;
/* DOCUMENT cobyla_init;
     Initialize COBYLA interface.  It is automatically called when COBYLA
     plugin is loaded but it can safely be called again to reinitialze
     constants.
   SEE ALSO: cobyla_create.
 */
cobyla_init;

/*---------------------------------------------------------------------------*/
