#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#include "brentq.h"

double brentq(callback_type f, double xa, double xb, double xtol, double rtol, int iter, stats_t *stats, void *args){
  int status;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;
  double x_lo = xa;
  double x_hi = xb;
  gsl_function F;
  double f_lo=(*f)(xa,args);
  double f_hi=(*f)(xb,args);
  if (f_lo*f_hi > 0) {
    stats->error_num = SIGNERR;
    return 0;
  }
  if (f_lo == 0) {
    return x_lo;
  }
  if (f_hi == 0) {
    return x_hi;
  }
  stats->error_num  = CONVERR;
  stats->funcalls   = 0; /* not used */
  stats->iterations = 0;
  F.function = f;
  F.params = args;
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  status=GSL_CONTINUE;
  for(int i=0;i<iter && status==GSL_CONTINUE;i++){
    stats->iterations++;
    status = gsl_root_fsolver_iterate (s);
    r = gsl_root_fsolver_root (s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo,x_hi,xtol,rtol);
    if (status == GSL_SUCCESS){
      stats->error_num=0;
      break;
    }
  }
  gsl_root_fsolver_free (s);
  return r;
}
