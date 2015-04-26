# C version of BOBYQA

This provides a C implementation of Mike Powell's BOBYQA algorithm for
minimizing a function of many variables.  The method is *derivatives free*
(only the function values are needed) and accounts for bound constraints on
the variables.  The algorithm is described in:

>  M.J.D. Powell, "The BOBYQA Algorithm for Bound Constrained Optimization
>  Without Derivatives."  Technical report, Department of Applied Mathematics
>  and Theoretical Physics, University of Cambridge (2009).

BOBYQA builds a quadratic model of the objective function from much less than
`(N+1)(N+2)/2` values of the function (with `N` the number of variables).  The
recommended number of points for building the quadratic model is `2*N+1`.  For
smooth objective functions, BOBYQA is expected to be more efficient than COBYLA
(which exploits a more simple linear model but implements arbitrary inequality
constraints).

In addition to being usable from C code, this version of BOBYQA has a few
improvements over the FORTRAN version:

* any objective function can be used (the function is an argument of the
  method) so different objective functions can co-exist in the same code
  (*e.g.* hierarchical optimization with BOBYQA is possible);

* a return code indicates the success of the method or the reason of the
  failure.

Of course, this version produces the same output for the tests and it has a
more C-like interface compared to a version produced by `f2c`.  A FORTRAN
wrapper is available.


## Usage

The code consists in two files [`bobyqa.c`](./bobyqa.c) and [`bobyqa.h`](./bobyqa.h)
and is documented in the header [`bobyqa.h`](./bobyqa.h).

To build the library, edit the `Makefile` to suit your preferences and
do:
```
make
```
Then copy `libbobyqa.a` and `bobyqa.h` to appropriate directories.

A number of macros can defined for the compiler to adapt the type of variables
used in the code (see [`bobyqa.h`](./bobyqa.h) for details).

You may check the code (this requires the original FORTRAN code available
from Mike Powell on request at mjdp@cam.ac.uk):
```
make test
```


## Legal notice

The present code is based on the original FORTRAN version written by Mike
Powell who kindly provides his code on demand (at mjdp@cam.ac.uk) and has
been converted to C by É. Thiébaut.

Copyright (c) 2009, Mike Powell (FORTRAN version).

Copyright (c) 2015, Éric Thiébaut (C version).

Read the accompanying [`LICENSE`](../LICENSE) file for details.
