#include <stdio.h>

#include "gx2.h"

int main(){
  double nt      =3;
  double coeffs[]={1,5,2};
  int    df[]    ={1,2,3};
  double lambda[]={2,3,7};
  gx2_stats_t stats;

  double x,x0,p,p0;

  x=25;
  double val=0.148139684865550; /* MATLAB value */

  p=gx2cdf(nt,x,coeffs,df,lambda,&stats);
  printf("p=%16.16g dev=%g exit=%d chi2_calls=%d err=%g\n",p,p-val,stats.error_num,stats.chi2_calls,stats.truncation_error);
  x0=gx2ppf(nt,p,coeffs,df,lambda,&stats);
  printf("x=%16.16g dev=%g exit=%d iterations=%d chi2_calls=%d funcalls=%d\n",x0,x0-x,stats.error_num,stats.iterations,stats.chi2_calls, stats.funcalls);

  p=0.65;
  x=gx2ppf(nt,p,coeffs,df,lambda,&stats);
  printf("x=%16.16g exit=%d iterations=%d chi2_calls=%d funcalls=%d\n",x,stats.error_num,stats.iterations,stats.chi2_calls, stats.funcalls);
  p0=gx2cdf(nt,x,coeffs,df,lambda,&stats);
  printf("p=%16.16g dev=%g exit=%d chi2_calls=%d err=%g\n",p0,p0-p,stats.error_num,stats.chi2_calls,stats.truncation_error);
}
