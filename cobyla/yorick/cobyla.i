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

local COBYLA_SUCCESS, COBYLA_CALC_FC, COBYLA_ROUNDING_ERRORS;
local COBYLA_TOO_MANY_EVALUATIONS, COBYLA_BAD_ADDRESS;
local COBYLA_CORRUPTED;
extern cobyla_create;
extern cobyla_iterate;
/* DOCUMENT ctx = cobyla_create(n, m, rhobeg, rhoend, iprint, maxfun);
         or status = cobyla_iterate(ctx, f, x, c);

     The function `cobyla_create` makes a new instance for Mike Powell's COBYLA
     algorithm for minimizing a function of a few variables.  The method is
     "derivatives free" (only the function values are needed) and accounts for
     constraints on the variables.  N is the number of variables, M is the
     number of (inequality) constraints, RHOBEG and RHOEND are the initial and
     final size of the trust region, IPRINT control the verbosity of the method
     (see below) and MAXFUN is the maximum number of function evaluations.

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
     >   while (ctx.status == COBYLA_CALC_FC) {
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
