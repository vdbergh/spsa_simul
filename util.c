#include "spsa_sim.h"

void params_disp(const char *prompt, int num_params, params_t *p){
  printf("%s",prompt);
  for(int i=0;i<num_params;i++){
    printf("%f  ",(*p)[i]);
  }
  printf("\n");
}

void params_from_string(const char *str_in,  params_t *p){
  char *str_copy=strdup(str_in);
  char *token;
  double pp;
  token=strtok(str_copy,":,");
  pp=atof(token);
  for(int i=0;i<MAX_PARAMS;i++){
    (*p)[i]=pp;
    if(token!=NULL){
      token=strtok(NULL,":,");
      if(token!=NULL){
	pp=atof(token);
      }
    }
  }
  free(str_copy);
}

int nproc(void){
  char *line=NULL;
  size_t len=0;
  ssize_t read;
  int nproc_;
  FILE *f=popen(NPROC_COMMAND,"r");
  read=getline(&line,&len,f);
  pclose(f);
  nproc_=read>0?atoi(line):1;
  free(line);
  return nproc_;
}
