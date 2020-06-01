#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#include "brentq.h"

typedef struct {
  callback_type *f;
  int funcalls;
  void *args;
} callback_wrapper_t;

double callback_wrapper(double x, void *args){
  callback_wrapper_t *c=(callback_wrapper_t *)(args);
  c->funcalls++;
  return (**(c->f))(x,c->args);
}

double brentq(callback_type f, double xa, double xb, double xtol, double rtol, int iter, stats_t *stats, void *args){
  int status;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;
  gsl_function F;
  callback_wrapper_t c;
  c.f=&f;
  c.funcalls=0;
  c.args=args;
  stats->error_num  = CONVERR;
  stats->iterations = 0;
  F.function = callback_wrapper;
  F.params = (void *)(&c);
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_set_error_handler_off();
  status=gsl_root_fsolver_set (s, &F, xa, xb);
  if(status==GSL_EINVAL){
    stats->error_num=SIGNERR;
  }else{
    status=GSL_CONTINUE;
  }
  for(int i=0;i<iter && status==GSL_CONTINUE;i++){
    stats->iterations++;
    gsl_root_fsolver_iterate (s);
    r = gsl_root_fsolver_root (s);
    xa = gsl_root_fsolver_x_lower (s);
    xb = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (xa, xb, xtol, rtol);
    if (status == GSL_SUCCESS){
      stats->error_num=0;
    }
  }
  stats->funcalls=c.funcalls;
  gsl_root_fsolver_free (s);
  return r;
}
