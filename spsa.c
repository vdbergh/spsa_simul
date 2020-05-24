#include "spsa_sim.h"

#define C 347.43558552260146

void spsa_disp(spsa_t *s){
  printf("num_params        =%d\n",s->num_params);
  printf("confidence        =%f\n",s->confidence);
  printf("draw_ratio        =%f\n",s->draw_ratio);
  printf("elo_per_parameter =%f\n",s->elo_per_parameter);
  printf("precision         =%f\n",s->precision);
  printf("c_ratio           =%f\n",s->c_ratio);
  printf("lambda_ratio      =%f\n",s->lambda_ratio);
  printf("r                 =%f\n",s->r);
  params_disp("c                 =",s->num_params,&(s->c));
  printf("num_games         =%d\n",s->num_games);
}

void spsa_init(spsa_t *s){
  s->num_params=0;
  s->confidence=0.95;
  s->draw_ratio=0.61;
  s->elo_per_parameter=2;
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

void spsa_compute(spsa_t *s, lfd_t *lfd){
  double chi2=gsl_cdf_chisq_Pinv(s->confidence, lfd->num_params);
  double H_diag,lambda;
  int ng;
  s->r=s->precision/(C*chi2*(1-s->draw_ratio)/8);
  s->num_params=lfd->num_params;
  s->num_games=0;
  for(int j=0;j<s->num_params;j++){
    s->c[j]=s->c_ratio*(lfd->maxima[j]-lfd->minima[j]);
    H_diag=-2*s->elo_per_parameter/pow((lfd->maxima[j]-lfd->minima[j])/2,2);
    lambda=-C/(2*s->r*pow(s->c[j],2)*H_diag);
    ng=(int)(s->lambda_ratio*lambda+0.5);
    if(ng>s->num_games){
      s->num_games=ng;
    }
  }
}

double spsa_elo_estimate(spsa_t *s, lfd_t *lfd, params_t *p0, double t){
  double ss=0;
  for(int j=0;j<s->num_params;j++){
    double ej,oj;
    ej=lfd->elos[j]/pow((lfd->maxima[j]-lfd->minima[j])/2,2);
    oj=lfd->optima[j];
    ss+=-ej*exp(-8*s->r*(pow(s->c[j],2))*ej*t/C)*pow(((*p0)[j]-oj),2);
  }
  return ss;
}

double spsa_noise_estimate(spsa_t *s, lfd_t *lfd, params_t *p0, double t){
  double ss=0;
  double front=(s->r)*(1-s->draw_ratio)*C/8;
  for(int j=0;j<s->num_params;j++){
    double ej,decay;
    ej=lfd->elos[j]/pow((lfd->maxima[j]-lfd->minima[j])/2,2);
    decay=8*s->r*(pow(s->c[j],2))*ej/C;
    ss-=front*(1-exp(-decay*t));
  }
  return ss;
}



