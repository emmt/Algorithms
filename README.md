# Algorithms

This repository is for various useful numerical algorithms which can be
used almost standalone.


## COBYLA

Directory [`cobyla`](cobyla) contains a C version of Mike Powell's COBYLA
(for **C**onstrained **O**ptimization **BY** **L**inear **A**pproximations)
algorithm for minimizing a function of many variables.  The method is
*derivatives free* (only the function values are needed) and accounts for
constraints on the variables.  The algorithm is described in:

> M.J.D. Powell, "A direct search optimization method that models the
> objective and constraint functions by linear interpolation," in Advances
> in Optimization and Numerical Analysis Mathematics and Its Applications,
> vol. 275 (eds. Susana Gomez and Jean-Pierre Hennart), Kluwer Academic
> Publishers, pp. 51-67 (1994).

A [Yorick](http://github.com/LLNL/yorick/) interface is also provided in
directory [`cobyla/yorick`](./cobyla/yorick).


## NEWUOA

Directory [`newuoa`](newuoa) provides a C implementation of Mike Powell's
NEWUOA algorithm for minimizing a function of many variables.  The method
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

A [Yorick](http://github.com/LLNL/yorick/) interface is also provided in
directory [`newuoa/yorick`](./newuoa/yorick).


## BOBYQA

Directory [`bobyqa`](bobyqa) provides a C implementation of Mike Powell's
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


## TOLMIN

Directory [`tolmin/src`](tolmin/src) provides a C implementation of Mike
Powell's TOLMIN algorithm for minimizing a function of many variables subject
to linear and bound constraints.  The algorithm is described in:

> M.J.D. Powell, "A tolerant algorithm for linearly constrained optimization
> calculations", Math. Programming B, Vol. 45, pp. 547-566 (1989).

The present code is based on the original FORTRAN version written by Mike
Powell who released his code under the GNU Lesser General Public License.  His
original code is available at
[CCPForge](https://ccpforge.cse.rl.ac.uk/gf/project/powell/).


## LINCOA

Directory [`lincoa/src`](lincoa/src) provides a C implementation of Mike
Powell's LINCOA algorithm for minimizing a function of many variables subject
to linear constraints.  The algorithm is described in:

> M.J.D. Powell, "On fast trust region methods for quadratic models with linear
> constraints", Report of the Department of Applied Mathematics and Theoretical
> Physics, Cambridge University, DAMTP 2014/NA02 (2014).

The present code is based on the original FORTRAN version written by Mike
Powell who released his code under the GNU Lesser General Public License.  His
original code is available at
[CCPForge](https://ccpforge.cse.rl.ac.uk/gf/project/powell/).


## L-BFGS-B

L-BFGS-B (for *Limited memory BFGS method with Bounds*) is an algorithm by
R.H. Byrd, P. Lu, J. Nocedal and C. Zhu (see references below) to minimize
a smooth function of many variables with simple bound constraints.  The
method requires the computation of the function and its gradient.  It exploits
a limited memory approximation of the function Hessian with BFGS updates.
It may be used to solve large scale problems.

Directory [`lbfgsb`](./lbfgsb) contains original code of L-BFGS-B with
simple C and Yorick wrapper code to use L-BFGS-B in these languages.

The L-BFGS-B algorithm is described in:

1. R.H. Byrd, P. Lu, J. Nocedal and C. Zhu, *A limited memory algorithm for
   bound constrained optimization*, SIAM J. Scientific Computing 16 (1995),
   no. 5, pp. 1190--1208.

2. C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, *L-BFGS-B: a limited memory
   FORTRAN code for solving bound constrained optimization problems*,
   Tech. Report, NAM-11, EECS Department, Northwestern University, 1994.

