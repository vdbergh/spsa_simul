#include "spsa_sim.h"

const char *options_usage_s="sprsa_simul [-h] [--num_params NUM_PARAMS] "
  "[--confidence CONFIDENCE] [--draw_ratio DRAW_RATIO] [--seed SEED] "
  "[--truncate TRUNCATE] [--bounds] [--precision PRECISION] [--c_ratio C_RATIO] "
  "[--lambda_ratio LAMBDA_RATIO] [--est_elos EST_ELOS,...] "
  "[--true_elos TRUE_ELOS1,...] [--minima MINIMA1,..] [--optima OPTIMA1,..] [--maxima MAXIMA1,...] "
  "[--start_elo START_ELO] [--quiet] [--threads THREADS]";

void options_usage(){
  printf("%s\n",options_usage_s);
}

void options_disp(options_t *o){
  printf("num_threads       =%d\n",o->num_threads);
  printf("truncate          =%d\n",o->truncate);
  printf("seed              =   %" PRIu64 "\n",o->seed);
  printf("start_elo         =%f\n",o->start_elo);
  printf("quiet             =%d\n",o->quiet);
}

int options_parse(int argc, char **argv, spsa_t *s, lfd_t *est_lfd, lfd_t *true_lfd, options_t *o){
  int ret=0;
  int num_params=1;
  params_t *true_elos_=NULL;
  params_t *est_elos_=NULL;
  params_t *minima_=NULL;
  params_t *optima_=NULL;
  params_t *maxima_=NULL;
  params_t true_elos;
  params_t est_elos;
  params_t minima;
  params_t optima;
  params_t maxima;
  o->seed=(uint64_t) time(0);
  o->truncate=-1;
  o->start_elo=2;
  o->num_threads=nproc();
  o->quiet=0;
  spsa_init(s);
  for(int i=1;i<=argc-1;i++){
    if(strcmp(argv[i],"-h")==0){
      return OPTIONS_PARSE_HELP;
    }else if(strcmp(argv[i],"--num_params")==0){
      if(i<argc-1){
	num_params=atoi(argv[i+1]);
	i++;
      }else{
	return OPTIONS_PARSE_HELP;
      }
    }else if(strcmp(argv[i],"--confidence")==0){
      if(i<argc-1){
	s->confidence=atof(argv[i+1]);
	if(s->confidence<=0|| s->confidence>=1){
	  return OPTIONS_PARSE_CONFIDENCE;
	}
	i++;
      }else{
	return OPTIONS_PARSE_CONFIDENCE;
      }
    }else if(strcmp(argv[i],"--draw_ratio")==0){
      if(i<argc-1){
	s->draw_ratio=atof(argv[i+1]);
	if(s->draw_ratio<=0|| s->draw_ratio>=1){
	  return OPTIONS_PARSE_DRAW_RATIO;
	}
	i++;
      }else{
	return OPTIONS_PARSE_DRAW_RATIO;
      }
    }else if(strcmp(argv[i],"--seed")==0){
      if(i<argc-1){
	o->seed=strtoull(argv[i+1],NULL,0);
	i++;
      }else{
	return OPTIONS_PARSE_SEED;
      }
    }else if(strcmp(argv[i],"--truncate")==0){
      if(i<argc-1){
	o->truncate=atoi(argv[i+1]);
	if(o->truncate<0){
	  return OPTIONS_PARSE_TRUNCATE;
	}
	i++;
      }else{
	return OPTIONS_PARSE_TRUNCATE;
      }
    }else if(strcmp(argv[i],"--bounds")==0){
      s->bounds=1;
    }else if(strcmp(argv[i],"--precision")==0){
      if(i<argc-1){
	s->precision=atof(argv[i+1]);
	if(s->precision<=0){
	  return OPTIONS_PARSE_PRECISION;
	}
	i++;
      }else{
	return OPTIONS_PARSE_PRECISION;
      }
    }else if(strcmp(argv[i],"--c_ratio")==0){
      if(i<argc-1){
	s->c_ratio=atof(argv[i+1]);
	if(s->c_ratio<=0 || s->c_ratio>=0){
	  return OPTIONS_PARSE_C_RATIO;
	}
	i++;
      }else{
	return OPTIONS_PARSE_C_RATIO;
      }
    }else if(strcmp(argv[i],"--lambda_ratio")==0){
      if(i<argc-1){
	s->lambda_ratio=atof(argv[i+1]);
	if(s->lambda_ratio<=0){
	  return OPTIONS_PARSE_LAMBDA_RATIO;
	}
	i++;
      }else{
	return OPTIONS_PARSE_LAMBDA_RATIO;
      }
    }else if(strcmp(argv[i],"--est_elos")==0){
      if(i<argc-1){
	params_from_string(argv[i+1],&est_elos);
	est_elos_=&est_elos;
	i++;
      }else{
	return OPTIONS_PARSE_EST_ELOS;
      }
    }else if(strcmp(argv[i],"--true_elos")==0){
      if(i<argc-1){
	params_from_string(argv[i+1],&true_elos);
	true_elos_=&true_elos;
	i++;
      }else{
	return OPTIONS_PARSE_TRUE_ELOS;
      }
    }else if(strcmp(argv[i],"--minima")==0){
      if(i<argc-1){
	params_from_string(argv[i+1],&minima);
	minima_=&minima;
	i++;
      }else{
	return OPTIONS_PARSE_MINIMA;
      }
    }else if(strcmp(argv[i],"--optima")==0){
      if(i<argc-1){
	params_from_string(argv[i+1],&optima);
	optima_=&optima;
	i++;
      }else{
	return OPTIONS_PARSE_OPTIMA;
      }
    }else if(strcmp(argv[i],"--maxima")==0){
      if(i<argc-1){
	params_from_string(argv[i+1],&maxima);
	maxima_=&maxima;
	i++;
      }else{
	return OPTIONS_PARSE_MAXIMA;
      }
    }else if(strcmp(argv[i],"--start_elo")==0){
      if(i<argc-1){
	o->start_elo=atof(argv[i+1]);
	if(o->start_elo<0){
	  return OPTIONS_PARSE_START_ELO;
	}
	i++;
      }else{
	return OPTIONS_PARSE_START_ELO;
      }
    }else if(strcmp(argv[i],"--quiet")==0){
      o->quiet=1;
    }else if(strcmp(argv[i],"--threads")==0){
      if(i<argc-1){
	o->num_threads=atoi(argv[i+1]);
	if(o->num_threads<1 || o->num_threads>MAX_THREADS){
	  return OPTIONS_PARSE_THREADS;
	}
	i++;
      }else{
	return OPTIONS_PARSE_THREADS;
      }
    }else{
      return OPTIONS_PARSE_UNKNOWN;
    }
  }
  ret=lfd_init(est_lfd,num_params,true_elos_,optima_,minima_,maxima_);
  if(ret!=0){
    return ret;
  }
  ret=lfd_init(true_lfd,num_params,est_elos_,optima_,minima_,maxima_);
  if(ret!=0){
    return ret;
  }
  return 0;
}
