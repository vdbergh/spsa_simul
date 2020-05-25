#include "spsa_sim.h"

int main(){
  sos_t sos={4,{1,1,1,1},{0,0,0,0},{1,1,1,1}};
  double x=sos_ppf(&sos, 0.95);
  printf("%10.10f\n",x);
  printf("%10.10f\n",sos_cdf(&sos, x));
}
