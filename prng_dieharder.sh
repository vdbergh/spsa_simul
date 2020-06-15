#!/bin/sh
gcc -Wall prng_dieharder.c prng.c prng.h -o prng_dieharder
./prng_dieharder | dieharder -g 200 -m 10 -a
