#include "spsa_sim.h"

double sos_cdf(sos_t *sos, double x){
  params_t coeffs;
  int df[MAX_PARAMS]; /* params_t is double[] */
  params_t lambda;
  double tol=1e-12;
  gx2_stats_t stats;
  for(int i=0;i<sos->num_params;i++){
    coeffs[i]=sos->coeffs[i]*sos->var[i];
    lambda[i]=pow(sos->mu[i],2)/sos->var[i];
    df[i]=1;
  }
  return gx2cdf(sos->num_params, x, coeffs, df, lambda, tol, &stats);
}

double sos_ppf(sos_t *sos, double p){
  params_t coeffs;
  int df[MAX_PARAMS]; /* params_t is double[] */
  params_t lambda;
  double tol=1e-12;
  gx2_stats_t stats;
  for(int i=0;i<sos->num_params;i++){
    coeffs[i]=sos->coeffs[i]*sos->var[i];
    lambda[i]=pow(sos->mu[i],2)/sos->var[i];
    df[i]=1;
  }
  double ret=gx2ppf(sos->num_params, p, coeffs, df, lambda, tol, &stats);
  return ret;
}

void sos_disp(sos_t *sos){
  printf("num_params =%d\n",sos->num_params);
  params_disp("coeffs     =",sos->num_params,&(sos->coeffs));
  params_disp("mu         =",sos->num_params,&(sos->mu));
  params_disp("var        =",sos->num_params,&(sos->var));
}
