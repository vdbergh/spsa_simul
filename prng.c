#include <time.h>
#include <assert.h>

#include "prng.h"

static const double   pow2_64 = 0x1.0000000000000p+64;
static const uint64_t pow10_9 = UINT64_C(1000000000);
static const uint64_t a       = UINT64_C(2862933555777941757);
static const uint64_t b       = UINT64_C(3037000493);
static const uint64_t a48     = UINT64_C(3311271626024157185);
static const uint64_t b48     = UINT64_C(8774982398954700800);

void prng_init(prng_t *prng){
  struct timespec t;
  int ret=clock_gettime(CLOCK_REALTIME,&t);
  assert(ret==0);
  *prng=pow10_9*((uint64_t) t.tv_sec)+((uint64_t) t.tv_nsec);
}

double prng_get(prng_t *prng){
/*
  https://nuclear.llnl.gov/CNP/rng/rngman/node4.html
*/
  uint64_t current=*prng;
  *prng=a*(*prng)+b;
  return current/pow2_64;
}

void prng_split(prng_t *master, prng_t *out){
  /* This is a poor man's splittable prng.
     It can only do a linear tree.
  */
  *master=a48*(*master)+b48;   /* do 2^48 steps */
  *out=*master;
}
