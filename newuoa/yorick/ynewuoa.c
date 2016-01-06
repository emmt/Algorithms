/*
 * ynewuoa.c -
 *
 * Implements Yorick interface to NEWUOA.
 */

/* Standard headers. */
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* Yorick headers. */
#include <pstdlib.h>
#include <yapi.h>

/* NEWUOA header. */
#include "newuoa.h"

static void push_string(const char* str);
static void set_global_int(const char* name, int value);

static void ynewuoa_free(void*);
#if 0
static void ynewuoa_print(void*);
static void ynewuoa_eval(void*, int);
#else
#  define ynewuoa_print NULL
#  define ynewuoa_eval NULL
#endif
static void ynewuoa_extract(void*, char*);

typedef struct _ynewuoa_object ynewuoa_object_t;
struct _ynewuoa_object {
  long n, npt;
  newuoa_context_t* ctx;
};

y_userobj_t ynewuoa_type = {
  "NEWUOA workspace", ynewuoa_free, ynewuoa_print,
  ynewuoa_eval, ynewuoa_extract, NULL
};

static void ynewuoa_free(void* addr)
{
  ynewuoa_object_t* obj = (ynewuoa_object_t*)addr;
  if (obj->ctx != NULL) {
    newuoa_delete(obj->ctx);
  }
}

static void ynewuoa_extract(void* addr, char* attr)
{
  ynewuoa_object_t* obj = (ynewuoa_object_t*)addr;
  if (attr != NULL) {
    int c = attr[0];
    if (c == 's') {
      if (strcmp(attr, "status") == 0) {
        ypush_int(newuoa_get_status(obj->ctx));
        return;
      }
    } else if (c == 'n') {
      c = attr[1];
      if (c == '\0') {
        ypush_long(obj->n);
        return;
      }
      if (c == 'p' && strcmp(attr, "npt") == 0) {
        ypush_long(obj->npt);
        return;
      }
      if (c == 'e' && strcmp(attr, "nevals") == 0) {
        ypush_long(newuoa_get_nevals(obj->ctx));
        return;
      }
    } else if (c == 'r') {
      c = attr[1];
      if (c == 'h' && strcmp(attr, "rho") == 0) {
        ypush_double(newuoa_get_rho(obj->ctx));
        return;
      }
      if (c == 'e' && strcmp(attr, "reason") == 0) {
        push_string(newuoa_reason(newuoa_get_status(obj->ctx)));
        return;
      }
    }
  }
  y_error("unknown member");
}

static ynewuoa_object_t*
ynewuoa_get_object(int iarg)
{
  return (ynewuoa_object_t*)yget_obj(iarg, &ynewuoa_type);
}

void Y_newuoa_create(int argc)
{
  long n, npt, iprint, maxfun;
  double rhobeg, rhoend;
  ynewuoa_object_t* obj;

  if (argc != 6) y_error("expecting exactly 5 arguments");
  n = ygets_l(5);
  npt = ygets_l(4);
  rhobeg = ygets_d(3);
  rhoend = ygets_d(2);
  iprint = ygets_l(1);
  maxfun = ygets_l(0);
  if (npt < n + 2 || npt > (n + 2)*(n + 1)/2) {
    y_error("NPT is not in the required interval");
  } else if (rhobeg <= rhoend || rhoend <= 0.0) {
    y_error("invalid RHOBEG and/or RHOEND");
  } else if (maxfun < 1) {
    y_error("invalid MAXFUN");
  }

  obj = (ynewuoa_object_t*)ypush_obj(&ynewuoa_type, sizeof(ynewuoa_object_t));
  obj->n = n;
  obj->npt = npt;
  obj->ctx = newuoa_create(n, npt, rhobeg, rhoend, iprint, maxfun);
  if (obj->ctx == NULL) {
    if (errno == ENOMEM) {
      y_error("insufficient memory");
    } else {
      y_error("unexpected error");
    }
  }
}

void Y_newuoa_iterate(int argc)
{
  double f;
  double* x;
  ynewuoa_object_t* obj;
  long n, index;
  long dims[Y_DIMSIZE];
  int iarg, save, status;

  if (argc != 3) y_error("expecting exactly 3 arguments");
  iarg = argc;

  /* Get object instance. */
  --iarg;
  obj = ynewuoa_get_object(iarg);

  /* Get function value. */
  --iarg;
  f = ygets_d(iarg);

  /* Get variables. */
  --iarg;
  index = yget_ref(iarg);
  if (index < 0) y_error("argument X must be a simple variable");
  x = ygeta_d(iarg, &n, dims);
  if (n != obj->n) y_error("bad number of variables");
  save = (yarg_scratch(iarg) ? iarg : -1);

  /* Call the iterator. */
  status = newuoa_iterate(obj->ctx, f, x);
  if (save >= 0) {
    yput_global(index, save);
  }
  ypush_int(status);
}

void Y_newuoa_restart(int argc)
{
  if (argc != 1) y_error("expecting exactly one argument");
  ypush_int(newuoa_restart(ynewuoa_get_object(0)->ctx));
}

void Y_newuoa_reason(int argc)
{
  if (argc != 1) y_error("expecting exactly one argument");
  push_string(newuoa_reason(ygets_i(0)));
}

void Y_newuoa_init(int argc)
{
#define SET_GLOBAL(x) set_global_int(#x, x)
  SET_GLOBAL(NEWUOA_ITERATE);
  SET_GLOBAL(NEWUOA_SUCCESS);
  SET_GLOBAL(NEWUOA_BAD_NPT);
  SET_GLOBAL(NEWUOA_ROUNDING_ERRORS);
  SET_GLOBAL(NEWUOA_TOO_MANY_EVALUATIONS);
  SET_GLOBAL(NEWUOA_STEP_FAILED);
  SET_GLOBAL(NEWUOA_BAD_ADDRESS);
  SET_GLOBAL(NEWUOA_CORRUPTED);
#undef SET_GLOBAL
  ypush_nil();
}

static void push_string(const char* str)
{
  ypush_q((long *)NULL)[0] = (str ? p_strcpy((char *)str) : NULL);
}

static void set_global_int(const char* name, int value)
{
  ypush_int(value);
  yput_global(yget_global(name, 0), 0);
  yarg_drop(1);
}

/*---------------------------------------------------------------------------*/
