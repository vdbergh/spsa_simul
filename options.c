#include "spsa_sim.h"

const char * options_messages[] = {
  NULL,
  "Illegal number of parameters",
  "Illegal elo list",
  "Illegal bounds",
  NULL,
  "Illegal confidence level",
  "Illegal draw ratio",
  "Illegal seed",
  "Illegal truncation",
  "Illegal precision",
  "Illegal c_ratio",
  "Illegal lambda_ratio",
  "Illegal estimated elos",
  "Illegal true elos",
  "Illegal minima",
  "Illegal true optima",
  "Illegal maxima",
  "Illegal estimated start elo",
  "Illegal true start elo",
  "Illegal number of threads",
  "Illegal heuristic",
  "Illegal number of parameters",
  "Unknown option"
};

const char *options_usage_s="spsa_simul [-h] [--num_params NUM_PARAMS] "
  "[--confidence CONFIDENCE] [--draw_ratio DRAW_RATIO] [--seed SEED] "
  "[--truncate TRUNCATE] [--bounds] [--precision PRECISION] [--c_ratio C_RATIO] "
  "[--lambda_ratio LAMBDA_RATIO] [--est_elos EST_ELOS1,...] "
  "[--true_elos TRUE_ELOS1,...] [--minima MINIMA1,..] [--true_optima OPTIMA1,..] [--maxima MAXIMA1,...] "
  "[--est_start_elo EST_START_ELO] [--true_start_elo TRUE_START_ELO] [--heuristic HEURISTIC] "
  "[--quiet] [--threads THREADS]";

void options_usage(){
  printf("%s\n",options_usage_s);
}

void options_disp(options_t *o){
  printf("num_threads       =%d\n",o->num_threads);
  printf("truncate          =%d\n",o->truncate);
  printf("true_start_elo    =%f\n",o->true_start_elo);
  printf("heuristic         =%d\n",o->heuristic);
  printf("seed              =   %" PRIu64 "\n",o->seed);
  printf("quiet             =%d\n",o->quiet);
}

int options_parse(int argc, char **argv, spsa_t *s, lf_t *est_lf, lf_t *true_lf, options_t *o, const char **option){
  int ret=0;
  int num_params=1;
  params_t *true_elos_=NULL;
  params_t *est_elos_=NULL;
  params_t *minima_=NULL;
  params_t *true_optima_=NULL;
  params_t *maxima_=NULL;
  params_t true_elos;
  params_t est_elos;
  params_t minima;
  params_t true_optima;
  params_t maxima;
  prng_init(&(o->seed));
  o->truncate=-1;
  o->true_start_elo=2;
  o->num_threads=nproc();
  o->quiet=0;
  o->heuristic=OPTIONS_HEURISTIC_LAMBDA_RATIO;
  spsa_init(s);
  for(int i=1;i<=argc-1;i++){
    const char *option_=argv[i];
    *option=option_;
    if(strcmp(option_,"-h")==0){
      return OPTIONS_PARSE_HELP;
    }else if(strcmp(option_,"--num_params")==0){
      if(i<argc-1){
	num_params=atoi(argv[i+1]);
	if(num_params<=0){
	  return OPTIONS_PARSE_NUM_PARAMS;
	}
	i++;
      }else{
	return OPTIONS_PARSE_NUM_PARAMS;
      }
    }else if(strcmp(option_,"--confidence")==0){
      if(i<argc-1){
	s->confidence=atof(argv[i+1]);
	if(s->confidence<=0|| s->confidence>=1){
	  return OPTIONS_PARSE_CONFIDENCE;
	}
	i++;
      }else{
	return OPTIONS_PARSE_CONFIDENCE;
      }
    }else if(strcmp(option_,"--draw_ratio")==0){
      if(i<argc-1){
	s->draw_ratio=atof(argv[i+1]);
	if(s->draw_ratio<=0|| s->draw_ratio>=1){
	  return OPTIONS_PARSE_DRAW_RATIO;
	}
	i++;
      }else{
	return OPTIONS_PARSE_DRAW_RATIO;
      }
    }else if(strcmp(option_,"--seed")==0){
      if(i<argc-1){
	o->seed=strtoull(argv[i+1],NULL,0);
	i++;
      }else{
	return OPTIONS_PARSE_SEED;
      }
    }else if(strcmp(option_,"--truncate")==0){
      if(i<argc-1){
	o->truncate=atoi(argv[i+1]);
	if(o->truncate<0){
	  return OPTIONS_PARSE_TRUNCATE;
	}
	i++;
      }else{
	return OPTIONS_PARSE_TRUNCATE;
      }
    }else if(strcmp(option_,"--bounds")==0){
      s->bounds=1;
    }else if(strcmp(option_,"--precision")==0){
      if(i<argc-1){
	s->precision=atof(argv[i+1]);
	if(s->precision<=0){
	  return OPTIONS_PARSE_PRECISION;
	}
	i++;
      }else{
	return OPTIONS_PARSE_PRECISION;
      }
    }else if(strcmp(option_,"--heuristic")==0){
      if(i<argc-1){
	o->heuristic=atoi(argv[i+1]);
	if(o->heuristic<=0 || o->heuristic>=3){
	  return OPTIONS_PARSE_HEURISTIC;
	}
	i++;
      }else{
	return OPTIONS_PARSE_HEURISTIC;
      }
    }else if(strcmp(option_,"--c_ratio")==0){
      if(i<argc-1){
	s->c_ratio=atof(argv[i+1]);
	if(s->c_ratio<=0 || s->c_ratio>=0.5){
	  return OPTIONS_PARSE_C_RATIO;
	}
	i++;
      }else{
	return OPTIONS_PARSE_C_RATIO;
      }
    }else if(strcmp(option_,"--lambda_ratio")==0){
      if(i<argc-1){
	s->lambda_ratio=atof(argv[i+1]);
	if(s->lambda_ratio<=0){
	  return OPTIONS_PARSE_LAMBDA_RATIO;
	}
	i++;
      }else{
	return OPTIONS_PARSE_LAMBDA_RATIO;
      }
    }else if(strcmp(option_,"--est_elos")==0){
      if(i<argc-1){
	params_from_string(argv[i+1],&est_elos);
	est_elos_=&est_elos;
	i++;
      }else{
	return OPTIONS_PARSE_EST_ELOS;
      }
    }else if(strcmp(option_,"--true_elos")==0){
      if(i<argc-1){
	params_from_string(argv[i+1],&true_elos);
	true_elos_=&true_elos;
	i++;
      }else{
	return OPTIONS_PARSE_TRUE_ELOS;
      }
    }else if(strcmp(option_,"--minima")==0){
      if(i<argc-1){
	params_from_string(argv[i+1],&minima);
	minima_=&minima;
	i++;
      }else{
	return OPTIONS_PARSE_MINIMA;
      }
    }else if(strcmp(option_,"--true_optima")==0){
      if(i<argc-1){
	params_from_string(argv[i+1],&true_optima);
	true_optima_=&true_optima;
	i++;
      }else{
	return OPTIONS_PARSE_TRUE_OPTIMA;
      }
    }else if(strcmp(option_,"--maxima")==0){
      if(i<argc-1){
	params_from_string(argv[i+1],&maxima);
	maxima_=&maxima;
	i++;
      }else{
	return OPTIONS_PARSE_MAXIMA;
      }
    }else if(strcmp(option_,"--est_start_elo")==0){
      if(i<argc-1){
	s->start_elo=atof(argv[i+1]);
	if(s->start_elo<0){
	  return OPTIONS_PARSE_EST_START_ELO;
	}
	i++;
      }else{
	return OPTIONS_PARSE_EST_START_ELO;
      }
    }else if(strcmp(option_,"--true_start_elo")==0){
      if(i<argc-1){
	o->true_start_elo=atof(argv[i+1]);
	if(o->true_start_elo<0){
	  return OPTIONS_PARSE_TRUE_START_ELO;
	}
	i++;
      }else{
	return OPTIONS_PARSE_TRUE_START_ELO;
      }
    }else if(strcmp(option_,"--quiet")==0){
      o->quiet=1;
    }else if(strcmp(option_,"--threads")==0){
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
  ret=lf_init(true_lf,num_params,true_elos_,true_optima_,minima_,maxima_);
  if(ret!=0){
    return ret;
  }
  ret=lf_init(est_lf,num_params,est_elos_,NULL,minima_,maxima_);
  if(ret!=0){
    return ret;
  }
  return 0;
}
