#include "spsa_sim.h"

void spsa_disp(spsa_t *s){
  printf("~~~~~~~design~~~~~~~ \n");
  printf("num_params        =%d\n",s->num_params);
  printf("confidence        =%f\n",s->confidence);
  printf("draw_ratio        =%f\n",s->draw_ratio);
  printf("precision         =%f\n",s->precision);
  printf("c_ratio           =%f\n",s->c_ratio);
  printf("lambda_ratio      =%f\n",s->lambda_ratio);
  printf("bounds            =%d\n",s->bounds);
  printf("~~~~~~~computed~~~~~~~ \n");
  printf("r                 =%f\n",s->r);
  params_disp("c                 =",s->num_params,&(s->c));
  printf("num_games         =%d\n",s->num_games);
}

void spsa_init(spsa_t *s){
  s->num_params=0;
  s->confidence=0.95;
  s->draw_ratio=0.61;
  s->precision=0.5;
  s->c_ratio=1.0/6.0;
  s->lambda_ratio=3;
  s->bounds=0; /* do not respect bounds by default */
  s->r=0;
  s->num_games=0;
  for(int i=0;i<MAX_PARAMS;i++){
    s->c[i]=0;
  }
}

void spsa_compute(spsa_t *s, lf_t *est_lf){
  double chi2=gsl_cdf_chisq_Pinv(s->confidence, est_lf->num_params);
  double H_diag,lambda;
  int ng;
  s->r=s->precision/(C*chi2*(1-s->draw_ratio)/8);
  s->num_params=est_lf->num_params;
  s->num_games=0;
  for(int j=0;j<s->num_params;j++){
    s->c[j]=s->c_ratio*(est_lf->maxima[j]-est_lf->minima[j]);
    H_diag=-2*est_lf->elos[j]/pow((est_lf->maxima[j]-est_lf->minima[j])/2,2);
    lambda=-C/(2*s->r*pow(s->c[j],2)*H_diag);
    ng=(int)(s->lambda_ratio*lambda+0.5);
    if(ng>s->num_games){
      s->num_games=ng;
    }
  }
}

void spsa_elo_estimate(spsa_t *s, lf_t *lf, params_t *p0, double t, double *fixed, double *noise, double *asymp){
  sos_t sos, sos_asymp;
  double fixed_, noise_;
  sos_from_lf_spsa(&sos, lf, s, p0, t);
  sos_from_lf_spsa(&sos_asymp, lf, s, p0, INFINITY);
  sos_expected(&sos,&fixed_,&noise_);
  *fixed=-fixed_;
  *noise=-noise_;
  sos_expected(&sos_asymp,&fixed_,&noise_);
  *asymp=-noise_;
}

double spsa_success_estimate(spsa_t *s, lf_t *lf, params_t *p0, double t){
  sos_t sos;
  sos_from_lf_spsa(&sos, lf, s, p0, t);
  return sos_cdf(&sos,s->precision);
}

double spsa_percentile(spsa_t *s, lf_t *lf, params_t *p0, double t, double p){
  sos_t sos;
  sos_from_lf_spsa(&sos, lf, s, p0, t);
  return -sos_ppf(&sos,p);
}

