/*
 * newuoa.i -
 *
 * Yorick interface to NEWUOA.
 *
 * Copyright (c) 2015, Éric Thiébaut
 *
 * This file is part of the `Algorithms` (https://github.com/emmt/Algorithms)
 * project.  Read the accompanying `LICENSE` file for details.
 */

if (is_func(plug_in)) plug_in, "ynewuoa";

local newuoa_maximize;
func newuoa_minimize(f, x0, rhobeg, rhoend, npt=, maxfun=, all=, verb=, err=)
/* DOCUMENT xmin = newuoa_minimize(f, x0, rhobeg, rhoend);
         or  obj = newuoa_minimize(f, x0, rhobeg, rhoend, all=1);
         or xmax = newuoa_maximize(f, x0, rhobeg, rhoend);
         or  obj = newuoa_maximize(f, x0, rhobeg, rhoend, all=1);

     Minimize or maximize the multi-variate function F starting at the initial
     point X0.  RHOBEG and RHOEND are the initial and final values of the trust
     region radius (0 < RHOEND <= RHOBEG).

     Note that the proper scaling of the variables is important for the success
     of the algorithm.  RHOBEG should be set to the typical size of the region
     to explorate and RHOEND should be set to the typical precision.

     Keyword NPT sets the number of of interpolation conditions.  Its default
     value is equal to 2*N+1 (the recommended value).

     Keyword MAXFUN sets the maximum number of function evaluations.  Its
     default value is equal to 30*N.

     If keyword ALL is true, the result is a structured object.  For
     `newuoa_minimize`, the members of the returned object are:

        obj.fmin   = minimal function value found
        obj.xmin   = corresponding parameters
        obj.nevals = number of function evaluations
        obj.status = status of the algorithm upon return
        obj.rho    = final radius of the trust region

     For `newuoa_maximize`, the two first members are:

        obj.fmax   = minimal function value found
        obj.xmax   = corresponding parameters

     Keyword ERR sets the behavior in case of abnormal termination.  If ERR=0,
     anything but a success throws an error (this is also the default
     behavior); if ERR > 0, non-fatal errors are reported by a warning message;
     if ERR < 0, non-fatal errors are silently ignored.

     Keyword VERB set the verbosity level.


   SEE ALSO: newuoa_create, newuoa_error.
 */
{
  x = double(unref(x0));
  n = numberof(x);
  ctx = newuoa_create(n, (is_void(npt) ? 2*n + 1 : npt), rhobeg, rhoend,
                      (is_void(verb) ? 0 : verb),
                      (is_void(maxfun) ? 30*n : maxfun));
  fmin = [];
  do {
    fx = f(x);
    if (is_void(fmin) || fmin > fx) {
      fmin = fx;
      xmin = x;
    }
    status = newuoa_iterate(ctx, fx, x);
  } while (status == NEWUOA_ITERATE);
  if (status != NEWUOA_SUCCESS) newuoa_error, status, err;
  return (all ? save(fmin, xmin, nevals = ctx.nevals,
                     status, rho = ctx.rho) : x);
}

func newuoa_maximize(f, x0, rhobeg, rhoend, npt=, maxfun=, all=, verb=, err=)
{
  x = double(unref(x0));
  n = numberof(x);
  ctx = newuoa_create(n, (is_void(npt) ? 2*n + 1 : npt), rhobeg, rhoend,
                      (is_void(verb) ? 0 : verb),
                      (is_void(maxfun) ? 30*n : maxfun));
  fmax = [];
  do {
    fx = f(x);
    if (is_void(fmax) || fmax < fx) {
      fmax = fx;
      xmax = x;
    }
    status = newuoa_iterate(ctx, -fx, x);
  } while (status == NEWUOA_ITERATE);
  if (status != NEWUOA_SUCCESS) newuoa_error, status, err;
  return (all ? save(fmax, xmax, nevals = ctx.nevals,
                     status, rho = ctx.rho) : x);
}

func newuoa_error(status, errmode)
/* DOCUMENT newuoa_error, status;
         or newuoa_error, status, errmode;

     Report an error in NEWUOA according to the value of STATUS.  Nothing is
     done if STATUS is NEWUOA_SUCCESS; otherwise, the optional argument ERRMODE
     determines the behavior.  If ERRMODE = 0, the routine throws an error
     (this is also the default behavior); if ERRMODE > 0, non-fatal errors are
     reported by a warning message; if ERRMODE < 0, non-fatal errors are
     silently ignored.

   SEE ALSO: newuoa_reason, error.
 */
{
  if (status != NEWUOA_SUCCESS) {
    if (errmode && (status == NEWUOA_ROUNDING_ERRORS ||
                    status == NEWUOA_TOO_MANY_EVALUATIONS ||
                    status == NEWUOA_STEP_FAILED)) {
      if (errmode > 0) {
        write, format="WARNING: Something wrong occured in NEWUOA: %s",
          newuoa_reason(status);
      }
    } else {
      error, swrite(format="Something wrong occured in NEWUOA: %s",
                    newuoa_reason(status));
    }
  }
}
errs2caller, newuoa_error;

local NEWUOA_ITERATE, NEWUOA_SUCCESS, NEWUOA_BAD_NPT, NEWUOA_ROUNDING_ERRORS;
local NEWUOA_TOO_MANY_EVALUATIONS, NEWUOA_STEP_FAILED, NEWUOA_BAD_ADDRESS;
local NEWUOA_CORRUPTED;
extern newuoa_create;
extern newuoa_iterate;
/* DOCUMENT ctx = newuoa_create(n, npt, rhobeg, rhoend, iprint, maxfun);
         or status = newuoa_iterate(ctx, f, x, c);

     The function `newuoa_create` makes a new instance for Mike Powell's NEWUOA
     algorithm for minimizing a function of many variables.  The method is
     "derivatives free" (only the function values are needed).

     N is the number of variables, NPT is the number of interpolation
     conditions. Its value must be in the interval [N+2,(N+1)(N+2)/2].  The
     ecommended number of points for building the quadratic model is NPT=2*N+1.

     RHOBEG and RHOEND are the initial and final values of a trust region
     radius, so both must be positive with RHOEND <= RHOBEG.  Typically RHOBEG
     should be about one tenth of the greatest expected change to a variable,
     and RHOEND should indicate the accuracy that is required in the final
     values of the variables.

     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
     amount of printing. Specifically, there is no output if IPRINT=0 and there
     is output only at the return if IPRINT=1. Otherwise, each new value of RHO
     is printed, with the best vector of variables so far and the corresponding
     value of the objective function.  Further, each new value of F with its
     variables are output if IPRINT=3.

     MAXFUN must be set to an upper bound on the number of objective function
     calls.

     The function `newuoa_iterate` performs an iteration of the algorithm given
     F the function value at the current variables X.

     Typical usage is:

     >   ctx = newuoa_create(n, npt, rhobeg, rhoend, iprint, maxfun);
     >   x = ...; // initial solution
     >   while (ctx.status == NEWUOA_ITERATE) {
     >     f = ...; // compute function value at X
     >     newuoa_iterate, ctx, f, x;
     >   }
     >   if (ctx.status != NEWUOA_SUCCESS) {
     >     error, swrite(format="Something wrong occured in NEWUOA: %s",
     >                   ctx.reason);
     >   }

     The context object CTX returned by the function `newuoa_create` has the
     following members:

         ctx.n       number of variables
         ctx.npt     number of intepolation points
         ctx.rho     radius of the trust region
         ctx.status  current status
         ctx.nevals  number of function evaluations so far
         ctx.reason  textual description of current status


   REFERENCES
     The NEWUOA algorithm is described in:

         M.J.D.  Powell, "The NEWUOA software for unconstrained minimization
         without derivatives", in Large-Scale Nonlinear Optimization, editors
         G. Di Pillo and M. Roma, Springer (2006), pages 255-297.

   SEE ALSO: newuoa_restart, newuoa_reason.
 */

extern newuoa_restart;
/* DOCUMENT newuoa_restart(ctx);
     Restart NEWUOA algorithm using the same parameters.  The return value is
     the new status of the algorithm, see `newuoa_get_status` for details.

   SEE ALSO: newuoa_create.
 */

extern newuoa_reason;
/* DOCUMENT newuoa_reason(status);
     Get a textual explanation of the status returned by NEWUOA.

   SEE ALSO: newuoa_create.
 */

extern newuoa_init;
/* DOCUMENT newuoa_init;
     Initialize NEWUOA interface.  It is automatically called when NEWUOA
     plugin is loaded but it can safely be called again to reinitialze
     constants.

   SEE ALSO: newuoa_create.
 */
newuoa_init;

/*
 * Local Variables:
 * mode: Yorick
 * tab-width: 8
 * indent-tabs-mode: nil
 * c-basic-offset: 2
 * fill-column: 79
 * coding: utf-8
 * ispell-local-dictionary: "american"
 * End:
 */
