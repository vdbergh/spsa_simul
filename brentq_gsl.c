#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#include "brentq.h"

double brentq(callback_type f, double xa, double xb, double xtol, double rtol, int iter, stats_t *stats, void *args){
  int status;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;
  gsl_function F;
  stats->error_num  = CONVERR;
  stats->funcalls   = 0; /* not used */
  stats->iterations = 0;
  F.function = f;
  F.params = args;
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
  gsl_root_fsolver_free (s);
  return r;
}
