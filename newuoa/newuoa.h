/*
 * newuoa.h -
 *
 * Definitions for Mike Powell's NEWUOA algorithm for minimizing a function
 * of many variables.  The method is "derivatives free" (only the function
 * values are needed).  The algorithm is described in:
 *
 *   M.J.D. Powell, "The NEWUOA software for unconstrained minimization
 *   without derivatives", in Large-Scale Nonlinear Optimization, editors
 *   G. Di Pillo and M. Roma, Springer (2006), pages 255-297.
 *
 * The present code is based on the original FORTRAN version written by Mike
 * Powell who kindly provides his code on demand (at mjdp@cam.ac.uk) and has
 * been converted to C by É. Thiébaut.
 *
 * Copyright (c) 2004, Mike Powell (FORTRAN version).
 * Copyright (c) 2015, Éric Thiébaut (C version).
 *
 * Read the accompanying `LICENSE` file for details.
 */

#ifndef _NEWUOA_H
#define  _NEWUOA_H 1

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

/* Prototype of the objective function assumed by the NEWUOA routine.  The
   returned value is the function value at X the current variables, N is the
   number of variables and DATA is anything needed by the function (unused by
   NEWUOA itself). */
typedef REAL newuoa_objfun(const INTEGER n, const REAL *x, void* data);


/* This subroutine seeks the least value of a function of many variables, by
   a trust region method that forms quadratic models by interpolation.  There
   can be some freedom in the interpolation conditions, which is taken up by
   minimizing the Frobenius norm of the change to the second derivative of
   the quadratic model, beginning with a zero matrix. The arguments of the
   subroutine are as follows.

   N must be set to the number of variables and must be at least two.  NPT is
   the number of interpolation conditions. Its value must be in the interval
   [N+2,(N+1)(N+2)/2].

   Initial values of the variables must be set in X(1),X(2),...,X(N). They
   will be changed to the values that give the least calculated F.  RHOBEG
   and RHOEND must be set to the initial and final values of a trust region
   radius, so both must be positive with RHOEND<=RHOBEG. Typically RHOBEG
   should be about one tenth of the greatest expected change to a variable,
   and RHOEND should indicate the accuracy that is required in the final
   values of the variables.

   The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
   amount of printing. Specifically, there is no output if IPRINT=0 and there
   is output only at the return if IPRINT=1. Otherwise, each new value of RHO
   is printed, with the best vector of variables so far and the corresponding
   value of the objective function. Further, each new value of F with its
   variables are output if IPRINT=3.

   MAXFUN must be set to an upper bound on the number of calls of CALFUN.
   The array W will be used for working space. Its length must be at least
   (NPT+13)*(NPT+N)+3*N*(N+3)/2.

   SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
   the value of the objective function for the variables X(1),X(2),...,X(N).
  */
extern int newuoa(const INTEGER n, const INTEGER npt,
                  newuoa_objfun* objfun, void* data,
                  REAL* x, const REAL rhobeg, const REAL rhoend,
                  const INTEGER iprint, const INTEGER maxfun,
                  REAL* w);


/* Possible values returned by NEWUOA. */
#define NEWUOA_SUCCESS                 (0) /* algorithm converged */
#define NEWUOA_BAD_NPT                (-1) /* NPT is not in the required
                                              interval */
#define NEWUOA_TOO_CLOSE              (-2) /* insufficient space between the
                                              bounds */
#define NEWUOA_ROUNDING_ERRORS        (-3) /* too much cancellation in a
                                              denominator */
#define NEWUOA_TOO_MANY_EVALUATIONS   (-4) /* maximum number of function
                                              evaluations exceeded */
#define NEWUOA_STEP_FAILED            (-5) /* a trust region step has failed to
                                              reduce Q */

/*---------------------------------------------------------------------------*/
/* FORTRAN WRAPPER */

/* Remark: Depending on your FORTRAN compiler, you may have to change the
           names of the compiled functions (it is assumed below that the
           link name is the name of the FORTRAN subroutine converted to
           lower case letters and with an underscore appended). */

/* This subroutine is a version of NEWUOA that is callable from FORTRAN code.
   The main difference with the C version is that the objective function
   must be provided by the following external subroutine.

   SUBROUTINE CALFUN (N,X,F) has to be provided by the user.  It must set F
   to the value of the objective function for the current values of the
   variables X(1),X(2),...,X(N), which are generated automatically
   by NEWUOA. */
extern int newuoa_(const INTEGER* n, const INTEGER* npt, REAL* x,
                   const REAL* rhobeg, const REAL* rhoend,
                   const INTEGER* iprint, const INTEGER* maxfun,
                   REAL* w);

/* Wrapper function to emulate `newuoa_calfun` function calling
   the user-defined `calfun_` subroutine. */
extern REAL newuoa_calfun_wrapper(const INTEGER n, const REAL* x, void* data);

/* Subroutine that must be defined by the application to use the FORTRAN
   wrapper to NEWUOA. */
extern int calfun_(const INTEGER* n, REAL *x, REAL *f);

/* Test problem for NEWUOA,FIXME: the objective function being the sum of the
   reciprocals of all pairwise distances between the points P_I,
   I=1,2,...,M in two dimensions, where M=N/2 and where the components of
   P_I are X(2*I-1) and X(2*I). Thus each vector X of N variables defines
   the M points P_I. The initial X gives equally spaced points on a
   circle. Four different choices of the pairs (N,NPT) are tried, namely
   (10,16), (10,21), (20,26) and (20,41). Convergence to a local minimum
   that is not global occurs in both the N=10 cases. The details of the
   results are highly sensitive to computer rounding errors. The choice
   IPRINT=2 provides the current X and optimal F so far whenever RHO is
   reduced. The bound constraints of the problem require every component of
   X to be in the interval [-1,1]. */
extern void newuoa_test(void);

#ifdef __cplusplus
}
#endif

#endif /* _NEWUOA_H */

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
