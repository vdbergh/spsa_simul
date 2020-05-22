#include "spsa_sim.h"

double myrand(uint64_t *prng){
  uint64_t a=UINT64_C(2862933555777941757);
  uint64_t b=UINT64_C(3037000493);
  uint64_t current=*prng;
  *prng=a*(*prng)+b;
  return current/pow(2,64);
}

void jump(uint64_t *prng){
  /* do 2^48 steps */
  uint64_t a=UINT64_C(3311271626024157185);
  uint64_t b=UINT64_C(8774982398954700800);
  *prng=a*(*prng)+b;
}
