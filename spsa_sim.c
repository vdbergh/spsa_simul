#include "spsa_sim.h"

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_t threads[MAX_THREADS];

void spsa_sim_step(uint64_t *prng, spsa_t *s,lf_t *lf,params_t *p){
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
  elo0=lf_eval(lf,&p_minus);
  elo1=lf_eval(lf,&p_plus);
  r=match(prng,s->draw_ratio,elo0,elo1);
  for(int j=0;j<s->num_params;j++){
    if(r==WIN){
      (*p)[j]+=s->r*s->c[j]*flips[j];
    }else if(r==LOSS){
      (*p)[j]-=s->r*s->c[j]*flips[j];
    }
  }
}

void spsa_sim(uint64_t *prng, spsa_t *s,lf_t *lf,params_t *p,int quiet){
  for(int i=0;i<s->num_games;i++){
    spsa_sim_step(prng,s,lf,p);
    if(i%100==0 && !quiet){
      printf("%d",i);
      for(int j=0; j<s->num_params;j++){
	printf("  %f",(*p)[j]);
      }
      printf("  %f\n",lf_eval(lf,p));
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
    spsa_sim(&prng,&(sim->s),&(sim->true_lf),&pp,1);
    elo=lf_eval(&(sim->true_lf),&pp);
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

void mainloop(sim_t *sim, int truncate){
  while(1){
    double p,ci,elo_avg;
    if(truncate>=0 && sim->count>=truncate){
      sim->stop=1;
      break;
    }
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
  lf_t est_lf;
  ret=options_parse(argc,argv,&(sim.s),&est_lf,&(sim.true_lf),&o);
  if(ret!=0){
    options_usage();
    return 0;
  }
  if(!o.quiet){
    printf("general options\n");
    printf("===============\n");
    options_disp(&o);
    printf("\n");
    printf("guessed loss function\n");
    printf("=====================\n");
    lf_disp(&est_lf);
    printf("\n");
    printf("true loss function\n");
    printf("==================\n");
    lf_disp(&sim.true_lf);
    printf("\n");
  }
  sim.prng=o.seed;
  sim.count=0;
  sim.pass=0;
  sim.elo_total=0;
  sim.stop=0;
  spsa_compute(&(sim.s),&est_lf);
  if(!o.quiet){
    printf("spsa data\n");
    printf("=========\n");
    spsa_disp(&(sim.s));
    printf("\n");
  }
  lf_start(&est_lf,o.start_elo,&(sim.p));
  if(!o.quiet){
    printf("starting point sims\n");
    printf("===================\n");
    params_disp("starting point      =",sim.s.num_params,&(sim.p));
    v=lf_eval(&(sim.true_lf),&(sim.p));
    printf("true loss function  =%f\n\n",v);
  }
  if(!o.quiet){
    double elo_est=spsa_elo_estimate(&(sim.s), &(sim.true_lf), &(sim.p), sim.s.num_games);
    double noise_est=spsa_noise_estimate(&(sim.s), &(sim.true_lf), &(sim.p), sim.s.num_games);
    printf("theoretical expected final Elo values (using the true loss function)\n");
    printf("====================================================================\n");
    printf("fixed part   =%f\n",elo_est);
    printf("noise part   =%f\n",noise_est);
    printf("total        =%f\n",elo_est+noise_est);
    printf("\n");
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
  mainloop(&sim,o.truncate);
  for(int i=0;i<o.num_threads;i++){
    pthread_join(threads[i], NULL);
  }
}
