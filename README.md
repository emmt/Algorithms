# Algorithms

This repository is for various useful numerical algorithms which can be
used almost standalone.


## COBYLA

Directory [cobyla](cobyla) contains a C version of Mike Powell's COBYLA
algorithm for minimizing a function of a few variables.  The method is
*derivatives free* (only the function values are needed) and accounts for
constraints on the variables.  The algorithm is described in:

>   M.J.D. Powell, "A direct search optimization method that models the
>   objective and constraint functions by linear interpolation," in Advances
>   in Optimization and Numerical Analysis Mathematics and Its Applications,
>   vol. 275 (eds. Susana Gomez and Jean-Pierre Hennart), Kluwer Academic
>   Publishers, pp. 51-67 (1994).
