#ifndef PRNG_H
#define PRNG_H

#include <inttypes.h>

typedef uint64_t prng_t;

void prng_init(prng_t *prng);
double prng_get(prng_t *prng);
void prng_split(prng_t *master, prng_t *out);

#endif /* PRNG_H */
