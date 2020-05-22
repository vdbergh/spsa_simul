#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_cdf.h>
#include <time.h>
#include <inttypes.h>
#include <pthread.h>
#include <unistd.h>

#define MAX_THREADS 128

#define NPROC_COMMAND "/usr/bin/nproc"

int nproc(void);

/*
  https://nuclear.llnl.gov/CNP/rng/rngman/node4.html
*/
double myrand(uint64_t *prng);

/* do 2^48 steps */
void jump(uint64_t *prng);

#define WIN 2
#define DRAW 1
#define LOSS 0

double L(double x);
int pick(uint64_t *prng, double w,double d,double l);
void wdl(double draw_ratio,double elo,double *wdl_out);

#define MAX_PARAMS 20
typedef double params_t[MAX_PARAMS];

/* loss function data */
typedef struct {
  int num_params;
  params_t elos;
  params_t optima;
  params_t minima;
  params_t maxima;
} lfd_t;

double loss_function(lfd_t *lfd, params_t *p);


/* TODO: make this enum */
#define LFD_INIT_NUM_PARAMS              1
#define LFD_INIT_ELOS                    2
#define LFD_INIT_BOUNDS                  3
#define OPTIONS_PARSE_HELP               4
#define OPTIONS_PARSE_CONFIDENCE         5
#define OPTIONS_PARSE_DRAW_RATIO         6
#define OPTIONS_PARSE_SEED               7
#define OPTIONS_PARSE_TRUNCATE           8
#define OPTIONS_PARSE_PRECISION          9
#define OPTIONS_PARSE_C_RATIO           10
#define OPTIONS_PARSE_LAMBDA_RATIO      11
#define OPTIONS_PARSE_ELO_PER_PARAMETER 12
#define OPTIONS_PARSE_ELOS              13
#define OPTIONS_PARSE_MINIMA            14
#define OPTIONS_PARSE_OPTIMA            15
#define OPTIONS_PARSE_MAXIMA            16
#define OPTIONS_PARSE_START_ELO         17
#define OPTIONS_PARSE_THREADS           18
#define OPTIONS_PARSE_UNKNOWN           19


int lfd_init(lfd_t *lfd, int num_params, params_t *elos, params_t *optima, params_t *minima, params_t *maxima);
void lfd_disp(lfd_t *lfd);
void lfd_start(lfd_t *lfd, double elo, params_t *p);

typedef struct {
  int num_params;
  double confidence;
  double draw_ratio;
  double elo_per_parameter;
  double precision;
  double c_ratio;
  double lambda_ratio;
  double bounds;
  /* computed */
  double r;
  params_t c;
  int num_games;
} spsa_t;

void spsa_init(spsa_t *s);
void spsa_disp(spsa_t *s);
void spsa_compute(spsa_t *s, lfd_t *lfd);

void params_disp(const char *prompt, int num_params, params_t *p);
void params_from_string(const char* str_in, params_t *p);

typedef struct {
  /* identical for every thread */
  spsa_t s;
  lfd_t  lfd;
  params_t p;
  /* in/out data */
  /* read and written with mutex */
           uint64_t prng;
  /* flag, read and written without mutex */
  /* volatile for safety */
  volatile int      stop;
  /* feedback variables, read without mutex */
  /* volatile for safety */
  volatile int      count;
  volatile int      pass;
  volatile double   elo_total;
} sim_t;

typedef struct {
  int num_threads;
  int truncate;
  uint64_t seed;
  double start_elo;
  int quiet;
} options_t;

int options_parse(int argc, char **argv, spsa_t *s, lfd_t *lfd, options_t *o);
void options_usage();
void options_disp(options_t *o);