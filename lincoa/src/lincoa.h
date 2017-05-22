#ifndef _LINCOA_H
#define _LINCOA_H 1

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/**
 * Prototype of the objective function assumed by the LINCOA routine.
 *
 * The returned value is the function value at X, the current variables, N is
 * the number of variables and DATA is anything else needed by the objective
 * function (unused by LINCOA itself).  Argument FEASIBLE is set true if the
 * current X satisfies the constraints to working accuracy.
 */
typedef double lincoa_objfun(const int n, const double x[], int feasible,
                             void* data);

/**
 * Get minimum of a multi-variate function subject to linear inequality
 * contraints.
 *
 * LINCOA seeks the least value of a function of many variables, subject to
 * general linear inequality constraints, by a trust region method that forms
 * quadratic models by interpolation. Usually there is much freedom in each new
 * model after satisfying the interpolation conditions, which is taken up by
 * minimizing the Frobenius norm of the change to the second derivative matrix
 * of the model. One new function value is calculated on each iteration,
 * usually at a point where the current model predicts a reduction in the least
 * value so far of the objective function subject to the linear constraints.
 * Alternatively, a new vector of variables may be chosen to replace an
 * interpolation point that may be too far away for reliability, and then the
 * new point does not have to satisfy the linear constraints.  The arguments of
 * the subroutine are as follows.
 *
 * @param n is the number of variables and must be at least two.
 *
 * @param npt is the number of interpolation conditions, which is required to
 *   be in the interval [N+2,(N+1)(N+2)/2]. Typical choices of the author are
 *   NPT=N+6 and NPT=2*N+1. Larger values tend to be highly inefficent when the
 *   number of variables is substantial, due to the amount of work and extra
 *   difficulty of adjusting more points.
 *
 * @param m is the number of linear inequality constraints.
 *
 * @param objfun is the objective function.
 *
 * @param data is anything needed by the objective function in addition to the
 *   variables and the number of variables.
 *
 * @param a is a matrix whose columns are the constraint gradients, which are
 *   required to be nonzero.
 *
 * @param ia is the first dimension of the array A, which must be at least N.*
 *
 * @param b is the vector of right hand sides of the constraints, the J-th
 *   constraint being that the scalar product of A(.,J) with X(.) is at
 *   most B(J). The initial vector X(.) is made feasible by increasing
 *   the value of B(J) if necessary.
 *
 * @param x is the vector of variables. Initial values of X(1),X(2),...,X(N)
 *   must be supplied. If they do not satisfy the constraints, then B
 *   is increased as mentioned above. X contains on return the variables
 *   that have given the least calculated F subject to the constraints.
 *
 * @param rhobeg is the initial trust region radius.
 *
 * @param rhoend is the final trust region radius.  @code rhobeg and @code
 *   rhoend must both be positive with RHOEND<=RHOBEG.  Typically, RHOBEG
 *   should be about one tenth of the greatest expected change to a variable,
 *   and RHOEND should indicate the accuracy that is required in the final
 *   values of the variables.
 *
 * @param the value of IPRINT should be set to 0, 1, 2 or 3, which controls the
 *   amount of printing. Specifically, there is no output if IPRINT=0 and
 *   there is output only at the return if IPRINT=1. Otherwise, the best
 *   feasible vector of variables so far and the corresponding value of
 *   the objective function are printed whenever RHO is reduced, where
 *   RHO is the current lower bound on the trust region radius. Further,
 *   each new value of F with its variables are output if IPRINT=3.
 *
 * @param maxfun must be set to an upper bound on the number of calls of CALFUN,
 *   its value being at least NPT+1.
 *
 * @param ws is an array used for working space. Its size must be at least the
 *   number of bytes given by @link lincoa_storage and properly aligned.  On
 *   return, ((double*)WS)[0] is set to the final value of F, and
 *   ((double*)WS)[1] is set to the total number of function evaluations plus
 *   0.5.
 *
 */
extern void lincoa(const int n, const int npt, const int m,
                   lincoa_objfun* objfun, void* data,
                   const double a[], const int ia,
                   const double b[], double x[],
                   const double rhobeg, const double rhoend,
                   const int iprint, const int maxfun, void* ws);

/**
 * Get the size of the workspace for LINCOA.
 *
 * LINCOA requires a workspace sufficient to store at least M*(2+N) +
 * NPT*(4+N+NPT) + N*(9+3*N) + MAX[ M+3*N, 2*M+N, 2*NPT ] floating-point values
 * plus N integer values.  This function computes the number of required bytes
 * taking care of padding for proper alignment.
 *
 * @param n   the number of variables.
 * @param npt the number of points.
 * @param m   the number of constraints.
 *
 * @return The number of bytes needed by LINCOA workspace (assuming the first
 *         bytes are correctly aligned for integers).
 */
extern size_t
lincoa_storage(int n, int npt, int m);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _LINCOA_H */
