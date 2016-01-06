/*
 * lbfgsb-wrap.c -
 *
 * C wrapper for L-BFGS-B FORTRAN routines.
 *
 *-----------------------------------------------------------------------------
 *
 * This software is licensed under the MIT "Expat" License.
 *
 * Copyright (C) 2002-2005, 2015: Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>
 *
 *-----------------------------------------------------------------------------
 */

#include <stdio.h>
#include <string.h>
#include "lbfgsb.h"

/* Define the following macro in order to remove trailing spaces when
   converting FORTRAN strings into C. */
#undef STR_TRIM

#undef STR_COPY

#ifdef STR_COPY
static char *str_copy(char *dst, const char *src, unsigned int len);
static char *str_copy(char *dst, const char *src, unsigned int len)
{
  unsigned int i;
  for (i = 0; i < len && src[i]; ++i) {
    dst[i] = src[i];
  }
  while (i < len) {
    dst[i++] = ' ';
  }
  dst[len] = 0;
  return dst;
}
#endif /* STR_COPY */

#ifdef STR_TRIM
static char *str_trim(char *s, unsigned int len);
static char *str_trim(char *s, unsigned int len)
{
  unsigned int i, j=len;
  for (i = 0; i < len; ++i) {
    if (s[i] == ' ') {
      if (i < j) {
        j = i;
      }
    } else {
      j = len;
    }
  }
  s[j] = 0;
  return s;
}
#endif /* STR_TRIM */

static char *str_pad(char *s, unsigned int len);
static char *str_pad(char *s, unsigned int len)
{
  unsigned int i;
  s[len] = 0;
  if ((i = strlen(s)) < len) {
    memset(s+i, ' ', len-i);
  }
  return s;
}

/* Prototype of FORTRAN subroutine to call. */
extern void lbsetu_(lbfgsb_integer_t* n, lbfgsb_integer_t* m,
                    double* x, const double* l, const double* u,
                    const lbfgsb_integer_t* bnd,
                    double* f, double* g, double* factr, double* pgtol,
                    double* wa, lbfgsb_integer_t* iwa, char* task,
                    lbfgsb_integer_t* iprint, char* csave, lbfgsb_logical_t* lsave,
                    lbfgsb_integer_t* isave, double* dsave,
                    lbfgsb_strlen_t task_len, lbfgsb_strlen_t csave_len);

int lbfgsb(lbfgsb_integer_t n, lbfgsb_integer_t m, double x[], double* f,
           double g[], const double l[], const double u[],
           const lbfgsb_integer_t bnd[], double factr, double pgtol,
           char task[], char csave[], lbfgsb_integer_t isave[],
           double dsave[], lbfgsb_integer_t iprint)
{
  /*
   * Sizes of arrays in the FORTRAN code:
   *
   *   double  x(n)
   *   double  l(n)
   *   double  u(n)
   *   double  g(n)
   *   integer bnd(n)
   *   integer iwa(3*n)
   *   integer isave(44)
   *   logical lsave(4)
   *   double  dsave(29)
   *   double  wa(2*m*n + 5*n + 11*m*m + 8*m)
   *
   * In the C code, workspace arrays are stored in two arrays:
   *   integer isave(4 + 44 + 3*n) = isave + lsave + iwa
   *   double  dsave(29 + 2*m*n + 5*n + 11*m*m + 8*m) = dsave + wa
   */
  lbfgsb_logical_t* lsave = (lbfgsb_logical_t*)&isave[44];
  lbfgsb_integer_t* iwa   = &isave[48];
  double*           wa    = &dsave[29];

  /* Call FORTRAN subroutine with character arrays fixed. */
  lbsetu_(&n, &m, x, l, u, bnd, f, g, &factr, &pgtol, wa, iwa,
          str_pad(task, LBFGSB_TASK_LENGTH), &iprint,
          str_pad(csave, LBFGSB_TASK_LENGTH),
          lsave, isave, dsave,
          (lbfgsb_strlen_t)LBFGSB_TASK_LENGTH,
          (lbfgsb_strlen_t)LBFGSB_TASK_LENGTH);

  /* Fix character arrays and compute return value (JOB). */
#ifdef STR_TRIM
  str_trim(task,  LBFGSB_TASK_LENGTH);
  str_trim(csave, LBFGSB_TASK_LENGTH);
#else
  task[LBFGSB_TASK_LENGTH] = 0;
  csave[LBFGSB_TASK_LENGTH] = 0;
#endif
  switch (task[0]) {
  case 'F':
    if (task[1] == 'G')
      return LBFGSB_FG;
    break;
  case 'N':
    if (strncmp(task, "NEW_X", 5) == 0)
      return LBFGSB_NEW_X;
    break;
  case 'C':
    if (strncmp(task, "CONVERGENCE", 11) == 0)
      return LBFGSB_CONVERGENCE;
    break;
  case 'W':
    if (strncmp(task, "WARNING", 7) == 0)
      return LBFGSB_WARNING;
    break;
  case 'S':
    if (strncmp(task, "START", 5) == 0)
      return LBFGSB_START;
    break;
  }
  return LBFGSB_ERROR;
}

/*---------------------------------------------------------------------------*/
