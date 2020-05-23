#include "spsa_sim.h"

double loss_function(lfd_t *lfd, params_t *p){
  double s=0;
  for(int i=0;i<lfd->num_params;i++){
    s+=-lfd->elos[i]*pow((*p)[i]-lfd->optima[i],2)/pow((lfd->maxima[i]-lfd->minima[i])/2,2);
  }
  return s;
}

void lfd_disp(lfd_t *lfd){
  printf("num_params       =%d\n",lfd->num_params);
  params_disp("elos             =",lfd->num_params,&(lfd->elos));
  params_disp("optima           =",lfd->num_params,&(lfd->optima));
  params_disp("minima           =",lfd->num_params,&(lfd->minima));
  params_disp("maxima           =",lfd->num_params,&(lfd->maxima));
}

int lfd_init(lfd_t *lfd, int num_params, params_t *elos, params_t *optima, params_t *minima, params_t *maxima){
  if(num_params<0 || num_params>MAX_PARAMS){
    return LFD_INIT_NUM_PARAMS;
  }
  lfd->num_params=num_params;
  if(elos==NULL){
    for(int i=0;i<num_params;i++){
      lfd->elos[i]=2;
    }
  }else{
    for(int i=0;i<num_params;i++){
      lfd->elos[i]=(*elos)[i];
    }
  }
  if(minima==NULL){
    for(int i=0;i<num_params;i++){
      lfd->minima[i]=0;
    }
  }
  if(optima==NULL && maxima==NULL){
    for(int i=0;i<num_params;i++){
      lfd->optima[i]=lfd->minima[i]+100;
      lfd->maxima[i]=lfd->minima[i]+200;
    }
  }else if(optima==NULL && maxima!=NULL){
    for(int i=0;i<num_params;i++){
      lfd->maxima[i]=(*maxima)[i];
      lfd->optima[i]=(lfd->minima[i]+lfd->maxima[i])/2;
    }
  }else if(optima!=NULL && maxima==NULL){
    for(int i=0;i<num_params;i++){
      lfd->optima[i]=(*optima)[i];
      lfd->maxima[i]=2*lfd->optima[i]-lfd->minima[i];
    }
  }else{
    for(int i=0;i<num_params;i++){
      lfd->optima[i]=(*optima)[i];
      lfd->maxima[i]=(*maxima)[i];
    }
  }
  for(int i=0;i<num_params;i++){
    if(lfd->elos[i]<=0){
      return LFD_INIT_ELOS;
    }else if(lfd->minima[i]>=lfd->maxima[i]){
      return LFD_INIT_BOUNDS;
    }
  }
  return 0;
}

void lfd_start(lfd_t *lfd, double elo, params_t *p){
  /* create a quadratic function of the form -u*l**2
     with l=0 corresponding to the otimum and l=1
     corresponding to the corner farthest from 
     the optimum
  */
  double u,l;
  params_t corner;
  for(int i=0;i<lfd->num_params;i++){
    if(lfd->maxima[i]-lfd->optima[i]>=lfd->optima[i]-lfd->minima[i]){
      corner[i]=lfd->maxima[i];
    }else{
      corner[i]=lfd->minima[i];
    }
  }
  u=-loss_function(lfd,&corner);
  l=sqrt(elo/u);
  for(int i=0;i<lfd->num_params;i++){
    (*p)[i]=lfd->optima[i]+l*(corner[i]-lfd->optima[i]);
  }
}
