#include "spsa_sim.h"

double L(double x){
  return 1/(1+pow(10,-x/400.0));
}

int pick(uint64_t *prng, double w,double d,double l){
  double s=myrand(prng);
  if(s<=w){
    return WIN;
  }else if(s<=w+d){
    return DRAW;
  }else{
    return LOSS;
  }
}

void wdl(double draw_ratio,double elo,double *wdl_out){
  double d=draw_ratio;
  double s_=L(elo);
  double w=s_-d/2.0;
  double l=1-d-w;
  wdl_out[WIN]=w;
  wdl_out[DRAW]=d;
  wdl_out[LOSS]=l;
}

