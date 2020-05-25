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
  double ret=gx2cdf(sos->num_params, x, coeffs, df, lambda, tol, &stats);
  assert(stats.error_num==GX2_CONVERGED);
  return ret;
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
  assert(stats.error_num==GX2_CONVERGED);
  return ret;
}

void sos_expected(sos_t *sos, double *fixed, double *noise){
  *fixed=0;
  *noise=0;
  for(int j=0;j<sos->num_params;j++){
    *fixed=sos->coeffs[j]*pow(sos->mu[j],2);
    *noise=sos->coeffs[j]*sos->var[j];
  }
}

void sos_disp(sos_t *sos){
  printf("num_params =%d\n",sos->num_params);
  params_disp("coeffs     =",sos->num_params,&(sos->coeffs));
  params_disp("mu         =",sos->num_params,&(sos->mu));
  params_disp("var        =",sos->num_params,&(sos->var));
}

void sos_from_lf_spsa(sos_t *sos, lf_t *lf, spsa_t *s, params_t *p, double t){
  sos->num_params=lf->num_params;
  for(int j=0;j<lf->num_params;j++){
    double ej,decayj,front;
    ej=lf->elos[j]/pow((lf->maxima[j]-lf->minima[j])/2,2);
    sos->coeffs[j]=ej;
    decayj=4*s->r*pow(s->c[j],2)*ej/C;
    sos->mu[j]=exp(-decayj*t)*((*p)[j]-lf->optima[j]);
    front=(s->r)*(1-s->draw_ratio)*C/8;
    sos->var[j]=front*(1-exp(-2*decayj*t))/ej;
  }
}
