#include "spsa_sim.h"

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_t threads[MAX_THREADS];

int sim_match(uint64_t *prng,spsa_t *s,double elo0,double elo1){
  double wdl_[3];
  double elo=elo1-elo0;
  wdl(s->draw_ratio,elo,wdl_);
  return pick(prng,wdl_[WIN],wdl_[DRAW],wdl_[LOSS]);
}

void spsa_sim_step(uint64_t *prng, spsa_t *s,lfd_t *lfd,params_t *p){
  params_t flips,p_minus,p_plus;
  double elo0,elo1;
  double r;
  for(int j=0;j<s->num_params;j++){
    r=myrand(prng);
    if(r>=0.5){
      flips[j]=1;
    }else{
      flips[j]=-1;
    }
    p_minus[j]=(*p)[j]-flips[j]*s->c[j];
    p_plus[j]=(*p)[j]+flips[j]*s->c[j];
  }
  elo0=loss_function(lfd,&p_minus);
  elo1=loss_function(lfd,&p_plus);
  r=sim_match(prng,s,elo0,elo1);
  for(int j=0;j<s->num_params;j++){
    if(r==WIN){
      (*p)[j]+=s->r*s->c[j]*flips[j];
    }else if(r==LOSS){
      (*p)[j]-=s->r*s->c[j]*flips[j];
    }
  }
}

void spsa_sim(uint64_t *prng, spsa_t *s,lfd_t *lfd,params_t *p,int quiet){
  for(int i=0;i<s->num_games;i++){
    spsa_sim_step(prng,s,lfd,p);
    if(i%100==0 && !quiet){
      printf("%d",i);
      for(int j=0; j<s->num_params;j++){
	printf("  %f",(*p)[j]);
      }
      printf("  %f\n",loss_function(lfd,p));
    }
  }
}

void* spsa_sims(void *args){
  sim_t *sim;
  sim=(sim_t *)(args);
  uint64_t prng;
  pthread_mutex_lock(&mutex);
  jump(&(sim->prng));
  prng=sim->prng;
  pthread_mutex_unlock(&mutex);
  while(!sim->stop){
    double elo;
    params_t pp;
    for(int j=0;j<(sim->s).num_params;j++){
      pp[j]=(sim->p)[j];
    }
    spsa_sim(&prng,&(sim->s),&(sim->lfd),&pp,1);
    elo=loss_function(&(sim->lfd),&pp);
    pthread_mutex_lock(&mutex);
    sim->count++;
    sim->elo_total+=elo;
    if(fabs(elo)<=(sim->s).precision){
      sim->pass++;
    }
    pthread_mutex_unlock(&mutex);
  }
  return NULL;
}

void mainloop(sim_t *sim){
  while(1){
    double p,ci,elo_avg;
    sleep(2);
    p=sim->pass/(sim->count+0.0);
    ci=3*sqrt(p*(1-p))/sqrt(sim->count);
    elo_avg=sim->elo_total/sim->count;
    printf("sims=%d success=%.4f[%.4f,%.4f] elo_avg=%f\n",sim->count,p,p-ci,p+ci,elo_avg);
    fflush(stdout);
  }
}




int main(int argc, char **argv){
  sim_t sim;
  double v;
  options_t o;
  int ret;
  ret=options_parse(argc,argv,&(sim.s),&(sim.lfd),&o);
  if(ret!=0){
    options_usage();
    return 0;
  }
  if(!o.quiet){
    printf("options\n");
    printf("=======\n");
    options_disp(&o);
    printf("\n");
    printf("loss function\n");
    printf("=============\n");
    lfd_disp(&sim.lfd);
    printf("\n");
  }
  sim.prng=o.seed;
  sim.count=0;
  sim.pass=0;
  sim.elo_total=0;
  sim.stop=0;
  spsa_compute(&(sim.s),&(sim.lfd));
  if(!o.quiet){
    printf("spsa data\n");
    printf("=========\n");
    spsa_disp(&(sim.s));
    printf("\n");
  }
  lfd_start(&(sim.lfd),o.start_elo,&(sim.p));
  if(!o.quiet){
    printf("starting point\n");
    printf("==============\n");
    params_disp("starting point =",sim.lfd.num_params,&(sim.p));
    v=loss_function(&(sim.lfd),&(sim.p));
    printf("loss_function  =%f\n\n",v);
  }
  if(o.quiet){
    printf("num_params        =%d\n",sim.s.num_params);
    printf("start_elo         =%.2f\n",o.start_elo);
    printf("num_games         =%d\n",sim.s.num_games);
    printf("r                 =%f\n",sim.s.r);
    params_disp("c                 =",sim.s.num_params,&(sim.s.c));
    printf("\n");
  }
  printf("sims\n");
  printf("====\n");
  for(int i=0;i<o.num_threads;i++){
    pthread_create(&(threads[i]), NULL, spsa_sims, (void*) (&sim));
  }
  mainloop(&sim);
  for(int i=0;i<o.num_threads;i++){
    pthread_join(threads[i], NULL);
  }
}
