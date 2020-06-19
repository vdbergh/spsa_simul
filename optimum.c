#include "spsa_sim.h"
#include "brentq.h"
#include "float.h"

/*
For a given r we find R such that the solution in R to

exp(-R*T)+R*(1-exp(-R*T))=rr

is minimal
*/

/* (reenabled when used)
static double loss(double R, double T){
  return exp(-R*T)+R*(1-exp(-R*T));
}
*/

static double T_from_R(double R, double rr){
  return log((R-1)/(R - rr))/R;
}

static double T_from_R_der(double R, double rr){
  /*
    Calculated by Sage.
    Should be cleaned up.
  */
  return (1/(R - rr) - R/pow(R - rr,2) + 1/pow(R - rr,2))/(R*(R/(R - rr) - 1/(R - rr))) - log(R/(R - rr) - 1/(R - rr))/pow(R,2);
}

static double f(double R,void *args){
  double *rr=(double *)(args);
  return T_from_R_der(R,*rr);
}

static void optimum_R_T(double rr, double *R, double *T){
  stats_t stats={0,0,0};
  double eps=1e-6;
  *R=brentq(f,eps,rr-eps,0,2*DBL_EPSILON,1000,&stats,(void*) &rr);
  assert(stats.error_num==0);
  *T=T_from_R(*R,rr);
}

void optimum_r_t(spsa_t *s, lf_t *lf, double *r, double *t){
  double b=lf->num_params*(1-s->draw_ratio)*C/8;
  double a=INFINITY;
  for(int i=0;i<lf->num_params;i++){
    double ei=lf->elos[i]/pow((lf->maxima[i]-lf->minima[i])/2,2);
    double ai=8*pow(s->c[i],2)*ei/C;
    if(ai< a){
      a=ai;
    }
  }
  double chi2=gsl_cdf_chisq_Pinv(s->confidence, lf->num_params);
  double p=s->precision/chi2*lf->num_params;
  double p0=s->start_elo;
  double rr=p/p0;
  double R,T;
  optimum_R_T(rr,&R,&T);
  *r=R*p0/b;
  *t=(T*b)/(a*p0);
}

