#include "spsa_sim.h"

double lf_eval(lf_t *lf, params_t *p){
  double s=0;
  for(int i=0;i<lf->num_params;i++){
    s+=-lf->elos[i]*pow((*p)[i]-lf->optima[i],2)/pow((lf->maxima[i]-lf->minima[i])/2,2);
  }
  return s;
}

void lf_disp(lf_t *lf){
  printf("+---------------+-----------------+-----------------+---------------+\n");
  printf("|      elos     |      minima     |      optima     |     maxima    |\n");
  printf("+---------------+-----------------+-----------------+---------------+\n");
  for(int j=0;j<lf->num_params;j++){
    printf("|  %10.2f   |   %10.2f    |   %10.2f    |   %10.2f  |\n",lf->elos[j],lf->minima[j],lf->optima[j],lf->maxima[j]);
  }
  printf("+---------------+-----------------+-----------------+---------------+\n");
}

int lf_init(lf_t *lf, int num_params, params_t *elos, params_t *optima, params_t *minima, params_t *maxima){
  if(num_params<0 || num_params>MAX_PARAMS){
    return LF_INIT_NUM_PARAMS;
  }
  lf->num_params=num_params;
  if(elos==NULL){
    for(int i=0;i<num_params;i++){
      lf->elos[i]=2;
    }
  }else{
    for(int i=0;i<num_params;i++){
      lf->elos[i]=(*elos)[i];
    }
  }
  for(int i=0;i<num_params;i++){
    if(minima==NULL){
      lf->minima[i]=0;
    }else{
      lf->minima[i]=(*minima)[i];
    }
  }
  if(optima==NULL && maxima==NULL){
    for(int i=0;i<num_params;i++){
      lf->optima[i]=lf->minima[i]+100;
      lf->maxima[i]=lf->minima[i]+200;
    }
  }else if(optima==NULL && maxima!=NULL){
    for(int i=0;i<num_params;i++){
      lf->maxima[i]=(*maxima)[i];
      lf->optima[i]=(lf->minima[i]+lf->maxima[i])/2;
    }
  }else if(optima!=NULL && maxima==NULL){
    for(int i=0;i<num_params;i++){
      lf->optima[i]=(*optima)[i];
      lf->maxima[i]=2*lf->optima[i]-lf->minima[i];
    }
  }else{
    for(int i=0;i<num_params;i++){
      lf->optima[i]=(*optima)[i];
      lf->maxima[i]=(*maxima)[i];
    }
  }
  for(int i=0;i<num_params;i++){
    if(lf->elos[i]<=0){
      return LF_INIT_ELOS;
    }else if(lf->minima[i]>=lf->maxima[i]){
      return LF_INIT_BOUNDS;
    }
  }
  return 0;
}

void lf_start(lf_t *lf, double elo, params_t *p){
  /* create a quadratic function of the form -u*l**2
     with l=0 corresponding to the otimum and l=1
     corresponding to the corner farthest from 
     the optimum
  */
  double u,l;
  params_t corner;
  for(int i=0;i<lf->num_params;i++){
    if(lf->maxima[i]-lf->optima[i]>=lf->optima[i]-lf->minima[i]){
      corner[i]=lf->maxima[i];
    }else{
      corner[i]=lf->minima[i];
    }
  }
  u=-lf_eval(lf,&corner);
  l=sqrt(elo/u);
  for(int i=0;i<lf->num_params;i++){
    (*p)[i]=lf->optima[i]+l*(corner[i]-lf->optima[i]);
  }
}
