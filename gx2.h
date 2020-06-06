#ifndef GX2_H
#define GX2_H

#define GX2_CONVERGED              0
#define GX2_NEGATIVE_TERMS         1
#define GX2_NEGATIVE_COEFFS        2
#define GX2_NEGATIVE_DF            3
#define GX2_NEGATIVE_LAMBDA        4
#define GX2CDF_UNDERFLOW           5
#define GX2CDF_NOT_CONVERGED       6
#define GX2PPF_SIGN_ERROR          7
#define GX2PPF_NOT_CONVERGED       8
#define GX2PPF_NOT_A_PROBABILITY   9

typedef struct {
  int error_num;
  int iterations;
  int funcalls;
  int chi2_calls;
  double truncation_error;
} gx2_stats_t;

void gx2_stats_disp(gx2_stats_t *stats);

double gx2cdf(int nt, double x, double *coeffs, int *df, double *lambda, gx2_stats_t *stats);
double gx2ppf(int nt, double p, double *coeffs, int *df, double *lambda, gx2_stats_t *stats);

#endif /* GX2_H */
