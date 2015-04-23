/*
 * cobyla.h -
 *
 * Definitions for Mike Powell's COBYLA algorithm for minimizing a function of
 * a few variables.  The method is "derivatives free" (only the function values
 * are needed) and accounts for constraints on the variables.  The algorithm is
 * described in:
 *
 *   M.J.D. Powell, "A direct search optimization method that models the
 *   objective and constraint functions by linear interpolation," in Advances
 *   in Optimization and Numerical Analysis Mathematics and Its Applications,
 *   vol. 275 (eds. Susana Gomez and Jean-Pierre Hennart), Kluwer Academic
 *   Publishers, pp. 51-67 (1994).
 *
 * The present code is based on the original FORTRAN version written by Mike
 * Powell who kindly provides his code on demand and has been converted to C by
 * É. Thiébaut.
 *
 * Copyright (c) 1992, Mike Powell (FORTRAN version).
 * Copyright (c) 2015, Éric Thiébaut (C version).
 */

#ifndef _COBYLA_H
#define  _COBYLA_H 1

#ifndef LOGICAL
# define LOGICAL int
#endif

#ifndef INTEGER
# define INTEGER long
#endif

#undef REAL
#ifdef SINGLE_PRECISION
# define REAL      float
#else
# define REAL      double
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Prototype of the objective function assumed by the COBYLA routine.  The
   returned value is the function value at X, N is the number of variables, M
   is the number of constraints, X are the current values of the variables and
   CON is to store the M constraints.  DATA is anything needed by the function
   (unused by COBYLA itself). */
typedef REAL cobyla_calcfc(INTEGER n, INTEGER m, const REAL x[],
			   REAL con[], void* data);

/*
 * This subroutine minimizes an objective function F(X) subject to M inequality
 * constraints on X, where X is a vector of variables that has N components.
 * The algorithm employs linear approximations to the objective and constraint
 * functions, the approximations being formed by linear interpolation at N+1
 * points in the space of the variables.  We regard these interpolation points
 * as vertices of a simplex.  The parameter RHO controls the size of the
 * simplex and it is reduced automatically from RHOBEG to RHOEND.  For each RHO
 * the subroutine tries to achieve a good vector of variables for the current
 * size, and then RHO is reduced until the value RHOEND is reached.  Therefore
 * RHOBEG and RHOEND should be set to reasonable initial changes to and the
 * required accuracy in the variables respectively, but this accuracy should be
 * viewed as a subject for experimentation because it is not guaranteed.  The
 * subroutine has an advantage over many of its competitors, however, which is
 * that it treats each constraint individually when calculating a change to the
 * variables, instead of lumping the constraints together into a single penalty
 * function.  The name of the subroutine is derived from the phrase Constrained
 * Optimization BY Linear Approximations.
 *
 * The user must set the values of N, M, RHOBEG and RHOEND, and must provide an
 * initial vector of variables in X.  Further, the value of IPRINT should be
 * set to 0, 1, 2 or 3, which controls the amount of printing during the
 * calculation.  Specifically, there is no output if IPRINT=0 and there is
 * output only at the end of the calculation if IPRINT=1.  Otherwise each new
 * value of RHO and SIGMA is printed.  Further, the vector of variables and
 * some function information are given either when RHO is reduced or when each
 * new value of F(X) is computed in the cases IPRINT=2 or IPRINT=3
 * respectively.  Here SIGMA is a penalty parameter, it being assumed that a
 * change to X is an improvement if it reduces the merit function:
 *
 *            F(X)+SIGMA*MAX(0.0,-C1(X),-C2(X),...,-CM(X)),
 *
 * where C1,C2,...,CM denote the constraint functions that should become
 * nonnegative eventually, at least to the precision of RHOEND.  In the
 * printed output the displayed term that is multiplied by SIGMA is
 * called MAXCV, which stands for 'MAXimum Constraint Violation'.  The
 * argument MAXFUN is an integer variable that must be set by the user to a
 * limit on the number of calls of CALCFC, the purpose of this routine being
 * given below.  The value of MAXFUN will be altered to the number of calls
 * of CALCFC that are made.  The arguments W and IACT provide real and
 * integer arrays that are used as working space.  Their lengths must be at
 * least N*(3*N+2*M+11)+4*M+6 and M+1 respectively.
 *
 * In order to define the objective and constraint functions, we require a
 * subroutine that has the name and arguments
 *
 *            SUBROUTINE CALCFC (N,M,X,F,CON)
 *            DIMENSION X(*),CON(*)  .
 *
 * The values of N and M are fixed and have been defined already, while X is
 * now the current vector of variables.  The subroutine should return the
 * objective and constraint functions at X in F and CON(1),CON(2), ...,CON(M).
 * Note that we are trying to adjust X so that F(X) is as small as possible
 * subject to the constraint functions being nonnegative.
 */
extern int cobyla(INTEGER n, INTEGER m,
                  cobyla_calcfc* calcfc, void* calcfc_data,
                  REAL x[], REAL rhobeg, REAL rhoend,
                  INTEGER iprint, INTEGER* maxfun, REAL w[], INTEGER iact[]);

/*
 * Subroutine version designed to be callable from FORTRAN code.  Usage is
 * similar to that of the above version except that everything is passed by
 * address and that, in order to define the objective and constraint functions,
 * we require a subroutine that has the name and arguments:
 *
 *            SUBROUTINE CALCFC (N,M,X,F,CON)
 *            DIMENSION X(*),CON(*)
 *
 * The values of N and M are fixed and have been defined already, while X is
 * now the current vector of variables.  The subroutine should return the
 * objective and constraint functions at X in F and CON(1),CON(2), ...,CON(M).
 * Note that we are trying to adjust X so that F(X) is as small as possible
 * subject to the constraint functions being nonnegative.
 */
extern int cobyla_(INTEGER* n, INTEGER* m, REAL x[],
                   REAL* rhobeg, REAL* rhoend, INTEGER* iprint,
                   INTEGER* maxfun, REAL w[], INTEGER iact[]);

/* Possible values returned by COBYLA. */
#define COBYLA_SUCCESS                 (0)
#define COBYLA_ROUNDING_ERRORS        (-1)
#define COBYLA_TOO_MANY_EVALUATIONS   (-2)

/* Get a textual explanation of the status returned by COBYLA. */
extern const char*
cobyla_reason(int status);

extern const char* cobyla_get_error(int code);

#ifdef __cplusplus
}
#endif

#endif /* _COBYLA_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * fill-column: 79
 * coding: utf-8
 * End:
 */
