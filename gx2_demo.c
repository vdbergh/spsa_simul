#include <stdio.h>

#include "gx2.h"

int main(){
  double val     =0.148139684865550; /* MATLAB value */
  double nt      =3;
  double x       =25;
  double coeffs[]={1,5,2};
  int    df[]    ={1,2,3};
  double lambda[]={2,3,7};
  gx2_stats_t stats;
  double p=gx2cdf(nt,x,coeffs,df,lambda,&stats);
  printf("p=%16.16f dev=%g exit=%d chi2_calls=%d err=%g\n",p,p-val,stats.error_num,stats.chi2_calls,stats.truncation_error);
  double x0=gx2ppf(nt,p,coeffs,df,lambda,&stats);
  printf("x=%16.16f dev=%g exit=%d iterations=%d chi2_calls=%d funcalls=%d\n",x0,x0-x,stats.error_num,stats.iterations,stats.chi2_calls, stats.funcalls);
}
