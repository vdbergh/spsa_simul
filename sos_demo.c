#include "spsa_sim.h"

int main(){
  int nt=12;
  sos_t sos={nt,{1,1,1,1,1,1,1,1,1,1,1,1},{0,0,0,0,0,0,0,0,0,0,0,0},{1,1,1,1,1,1,1,1,1,1,1,1}};
  double p_input=0.95;
  double x0= gsl_cdf_chisq_Pinv(p_input, nt);
  double x=sos_ppf(&sos, p_input);
  double chi2=gsl_cdf_chisq_P(x,nt);
  double chi2_=gsl_cdf_chisq_P(x0,nt);
  double p=sos_cdf(&sos, x);
  printf("p_in=%.16g (next line: x=ppf, p=cdf)\n",p_input);
  printf("x(p_in)=%.16g x_gsl(p_in)=%.16g rel_dev=%g, p(x)=%.16g p_gsl(x)=%.16g dev=%g p_gsl(x_gsl)=%.16g\n",x,x0,(x-x0)/x,p,chi2,p-chi2,chi2_);
}
