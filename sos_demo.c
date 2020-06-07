#include "spsa_sim.h"

int main(){
  sos_t sos={8,{1,1,1,1,1,1,1,1},{0,0,0,0,0,0,0,0},{1,1,1,1,1,1,1,1}};
  double x=sos_ppf(&sos, 0.95);
  double chi2=gsl_cdf_chisq_P(x,8);
  double p=sos_cdf(&sos, x);
  printf("p_input=%17.17g\n",0.95);
  printf("x=%17.17g p=%17.17g chi2=%17.17g dev=%g\n",x,p,chi2,p-chi2);
}
