/* This file is a slightly modified version of one provided by SciPy */

/* Written by Charles Harris charles.harris@sdl.usu.edu */

/* Modified to not depend on Python everywhere by Travis Oliphant.
 */

#ifndef BRENTQ_H
#define BRENTQ_H

typedef struct {
    int funcalls;
    int iterations;
    int error_num;
} stats_t;

#define SIGNERR -1
#define CONVERR -2

typedef double (*callback_type)(double,void*);
extern double brentq(callback_type f, double xa, double xb, double xtol, double rtol, int iter, stats_t *stats, void *args);

#endif
