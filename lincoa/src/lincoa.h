#ifndef _LINCOA_H
#define _LINCOA_H 1

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* Prototype of the objective function assumed by the LINCOA routine.  The
   returned value is the function value at X, the current variables, N is the
   number of variables and DATA is anything else needed by the objective
   function (unused by LINCOA itself). */
typedef double lincoa_objfun(const int n, const double x[], void* data);

extern void lincoa(const int n, const int npt, const int m,
                   lincoa_objfun* objfun, void* data,
                   const double a[], const int ia,
                   const double b[], double x[],
                   const double rhobeg, const double rhoend,
                   const int iprint, const int maxfun, void* ws);

/**
 * Get the size of the workspace for LINCOA.
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
