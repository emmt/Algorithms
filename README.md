# Algorithms

This repository is for various useful numerical algorithms which can be
used almost standalone.


## COBYLA

Directory [cobyla](cobyla) contains a C version of Mike Powell's COBYLA
(for **C**onstrained **O**ptimization **BY** **L**inear **A**pproximations)
algorithm for minimizing a function of many variables.  The method is
*derivatives free* (only the function values are needed) and accounts for
constraints on the variables.  The algorithm is described in:

> M.J.D. Powell, "A direct search optimization method that models the
> objective and constraint functions by linear interpolation," in Advances
> in Optimization and Numerical Analysis Mathematics and Its Applications,
> vol. 275 (eds. Susana Gomez and Jean-Pierre Hennart), Kluwer Academic
> Publishers, pp. 51-67 (1994).


## NEWUOA

Directory [newuoa](newuoa) provides a C implementation of Mike Powell's
NEWUOA algorithm for minimizing a function of a few variables.  The method
is *derivatives free* (only the function values are needed) and accounts
for bound constraints on the variables.  The algorithm is described in:

> M.J.D. Powell, "The NEWUOA software for unconstrained minimization
> without derivatives", in Large-Scale Nonlinear Optimization, editors
> G. Di Pillo and M. Roma, Springer (2006), pages 255-297.

NEWUOA builds a quadratic model of the objective function from much less
than `(N+1)(N+2)/2` values of the function (with `N` the number of
variables).  The recommended number of points for building the quadratic
model is `2*N+1`.  For smooth objective functions, NEWUOA is expected to be
more efficient than COBYLA (which exploits a more simple linear model but
implements arbitrary inequality constraints while NEWUOA is unconstrained).
If you have bound constraints, you may consider using BOBYQA instead.


## BOBYQA

Directory [bobyqa](bobyqa) provides a C implementation of Mike Powell's
BOBYQA (for **B**ound constrained **O**ptimization **BY** **Q**uadratic
**A**pproximations) algorithm for minimizing a function of many variables.
The method is *derivatives free* (only the function values are needed) and
accounts for bound constraints on the variables.  The algorithm is
described in:

> M.J.D. Powell, "The BOBYQA Algorithm for Bound Constrained Optimization
> Without Derivatives."  Technical report, Department of Applied
> Mathematics and Theoretical Physics, University of Cambridge (2009).

BOBYQA builds a quadratic model of the objective function from much less
than `(N+1)(N+2)/2` values of the function (with `N` the number of
variables).  The recommended number of points for building the quadratic
model is `2*N+1`.  For smooth objective functions, BOBYQA is expected to be
more efficient than COBYLA (which exploits a more simple linear model but
implements arbitrary inequality constraints).

