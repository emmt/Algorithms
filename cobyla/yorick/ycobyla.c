/*
 * ycobyla.c -
 *
 * Implements Yorick interface to COBYLA.
 */

/* Standard headers. */
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* Yorick headers. */
#include <pstdlib.h>
#include <yapi.h>

/* COBYLA header. */
#include "cobyla.h"

static void push_string(const char* str);
static void set_global_int(const char* name, int value);

static void ycobyla_free(void*);
#if 0
static void ycobyla_print(void*);
static void ycobyla_eval(void*, int);
#else
#  define ycobyla_print NULL
#  define ycobyla_eval NULL
#endif
static void ycobyla_extract(void*, char*);

typedef struct _ycobyla_object ycobyla_object_t;
struct _ycobyla_object {
  long n, m;
  cobyla_context_t* ctx;
};

y_userobj_t ycobyla_type = {
  "COBYLA workspace", ycobyla_free, ycobyla_print,
  ycobyla_eval, ycobyla_extract, NULL
};

static void ycobyla_free(void* addr)
{
  ycobyla_object_t* obj = (ycobyla_object_t*)addr;
  if (obj->ctx != NULL) {
    cobyla_delete(obj->ctx);
  }
}

static void ycobyla_extract(void* addr, char* attr)
{
  ycobyla_object_t* obj = (ycobyla_object_t*)addr;
  if (attr != NULL) {
    int c = attr[0];
    if (c == 's') {
      if (strcmp(attr, "status") == 0) {
        ypush_int(cobyla_get_status(obj->ctx));
        return;
      }
    } else if (c == 'n') {
      if (attr[1] == '\0') {
        ypush_long(obj->n);
        return;
      }
      if (strcmp(attr, "nevals") == 0) {
        ypush_long(cobyla_get_nevals(obj->ctx));
        return;
      }
    } else if (c == 'm') {
      if (attr[1] == '\0') {
        ypush_long(obj->m);
        return;
      }
    } else if (c == 'r') {
      if (strcmp(attr, "rho") == 0) {
        ypush_double(cobyla_get_rho(obj->ctx));
        return;
      }
      if (strcmp(attr, "reason") == 0) {
        push_string(cobyla_reason(cobyla_get_status(obj->ctx)));
        return;
      }
    }
  }
  y_error("unknown member");
}

static ycobyla_object_t*
ycobyla_get_object(int iarg)
{
  return (ycobyla_object_t*)yget_obj(iarg, &ycobyla_type);
}

void Y_cobyla_create(int argc)
{
  long n, m, iprint, maxfun;
  double rhobeg, rhoend;
  ycobyla_object_t* obj;

  if (argc != 6) y_error("expecting exactly 5 arguments");
  n = ygets_l(5);
  m = ygets_l(4);
  rhobeg = ygets_d(3);
  rhoend = ygets_d(2);
  iprint = ygets_l(1);
  maxfun = ygets_l(0);

  obj = (ycobyla_object_t*)ypush_obj(&ycobyla_type, sizeof(ycobyla_object_t));
  obj->n = n;
  obj->m = m;
  obj->ctx = cobyla_create(n, m, rhobeg, rhoend, iprint, maxfun);
  if (obj->ctx == NULL) {
    if (errno == ENOMEM) {
      y_error("insufficient memory");
    } else {
      y_error("illegal argument");
    }
  }
}

void Y_cobyla_iterate(int argc)
{
  double f;
  double* x;
  double* c;
  ycobyla_object_t* obj;
  long n, m, index;
  long dims[Y_DIMSIZE];
  int iarg, save, status;

  if (argc < 3 || argc > 4) y_error("expecting 3 or 4 arguments");
  iarg = argc;

  /* Get object instance. */
  --iarg;
  obj = ycobyla_get_object(iarg);

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

  /* Get array of constraints. */
  --iarg;
  if (iarg >= 0 && ! yarg_nil(iarg)) {
    c = ygeta_d(iarg, &m, NULL);
  } else {
    c = NULL;
    m = 0;
  }
  if (m != obj->m) y_error("bad number of constraints");

  /* Call the iterator. */
  status = cobyla_iterate(obj->ctx, f, x, c);
  if (save >= 0) {
    yput_global(index, save);
  }
  ypush_int(status);
}

void Y_cobyla_restart(int argc)
{
  if (argc != 1) y_error("expecting exactly one argument");
  ypush_int(cobyla_restart(ycobyla_get_object(0)->ctx));
}

void Y_cobyla_reason(int argc)
{
  if (argc != 1) y_error("expecting exactly one argument");
  push_string(cobyla_reason(ygets_i(0)));
}

void Y_cobyla_init(int argc)
{
#define SET_GLOBAL(x) set_global_int(#x, x)
  SET_GLOBAL(COBYLA_ITERATE);
  SET_GLOBAL(COBYLA_SUCCESS);
  SET_GLOBAL(COBYLA_ROUNDING_ERRORS);
  SET_GLOBAL(COBYLA_TOO_MANY_EVALUATIONS);
  SET_GLOBAL(COBYLA_BAD_ADDRESS);
  SET_GLOBAL(COBYLA_CORRUPTED);
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
