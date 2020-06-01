#include <time.h>
#include <assert.h>

#include "prng.h"

static const double   pow2_64 = 0x1.0000000000000p+64;
static const uint64_t pow10_9 = UINT64_C(1000000000);
static const uint64_t a       = UINT64_C(2862933555777941757);
static const uint64_t b       = UINT64_C(3037000493);
static const uint64_t a48     = UINT64_C(3311271626024157185);
static const uint64_t b48     = UINT64_C(8774982398954700800);

prng_t prng_seed(void){
  struct timespec t;
  int ret=clock_gettime(CLOCK_REALTIME,&t);
  assert(ret==0);
  return pow10_9*((uint64_t) t.tv_sec)+((uint64_t) t.tv_nsec);
}

double prng_get(prng_t *prng){
/*
  https://nuclear.llnl.gov/CNP/rng/rngman/node4.html
*/
  uint64_t current=*prng;
  *prng=a*(*prng)+b;
  return current/pow2_64;
}

void prng_jump(prng_t *prng){
  /* do 2^48 steps */
  *prng=a48*(*prng)+b48;
}
