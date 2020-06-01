#include "spsa_sim.h"

static const double pow2_64=0x1.0000000000000p+64;
static const uint64_t a=UINT64_C(2862933555777941757);
static const uint64_t b=UINT64_C(3037000493);
static const uint64_t a48=UINT64_C(3311271626024157185);
static const uint64_t b48=UINT64_C(8774982398954700800);

double myrand(uint64_t *prng){
  uint64_t current=*prng;
  *prng=a*(*prng)+b;
  return current/pow2_64;
}

void jump(uint64_t *prng){
  /* do 2^48 steps */
  *prng=a48*(*prng)+b48;
}
