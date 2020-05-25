#include <math.h>
#include <stdio.h>
#include <float.h>
#include <gsl/gsl_cdf.h>
#include <assert.h>

#include "brentq.h"

#include "gx2.h"


static const int K=1000;

double gx2cdf(int nt, double x, double *coeffs, int *df, double *lambda, double tol, gx2_stats_t *stats){
  /*
    Returns the CDF of a generalized chi-squared (a weighted sum of
    non-central chi-squares with positive weights), using Ruben's
    [1962] algorithm. See (3.1)(3.5)(3.6) in:

    Harold Ruben, Probability content of regions under spherical normal distributions IV:
    The distribution of homogeneous and non-homogeneous quadratic functions of normal
    variables.

    Example:
    double nt      =3;
    double x       =25;
    double coeffs[]={1,5,2};
    int    df[]    ={1,2,3};
    double lambda[] ={2,3,7};
    double c       =0;
    int    K       =100;
    double p;
    cx2cdf stats_t stats;

    double p=gx2cdf(nt,x,coeffs,df,lambda,tol,&stats);

    Inputs:
    nt        number of iterations in the weighted sum
    x         point at which to evaluate the CDF
    coeffs    row vector of coefficients of the non-central chi-squares
    df        row vector of degrees of freedom of the non-central chi-squares
    lambda     row vector of non-centrality parameters (sum of squares of
              means) of the non-central chi-squares
    tol       requested tolerance
    Outputs:
    p         computed CDF
    stats     statistcs (convergence, chi2_calls, truncation error)

    Author:
    Abhranil Das <abhranil.das@utexas.edu>
    Center for Perceptual Systems, University of Texas at Austin
    If you use this code, you may cite:
    A new method to compute classification error
    jov.arvojournals.org/article.aspx?articleid=2750251

    Hand translated to C by Michel Van den Bergh.
  */

  stats->error_num=GX2CDF_NOT_CONVERGED;
  stats->iterations=0; /* not used */
  stats->chi2_calls=0; 
  stats->truncation_error=0;
  stats->funcalls=0; /* not used */

  if(nt<=0){
    stats->error_num=GX2CDF_NEGATIVE_ITERATIONS;
    return 0;
  }
  
  
  for(int i=0;i<nt;i++){
    if(coeffs[i]<=0){
       stats->error_num=GX2CDF_NEGATIVE_COEFFS;
       return 0;
    }
    if(df[i]<=0){
       stats->error_num=GX2CDF_NEGATIVE_DF;
       return 0;
    }
    if(lambda[i]<0){
       stats->error_num=GX2CDF_NEGATIVE_LAMBDA;
       return 0;
    }
  }

  double coeffs_min=-1;
  for(int i=0;i<nt;i++){
    if(coeffs_min<0 || coeffs[i]<coeffs_min){
      coeffs_min=coeffs[i];
    }
  }

  double beta=0.90625*coeffs_min; /* Can be arbitrary in loc. cit. */

  double M=0;
  double D=0;
  for(int i=0;i<nt;i++){
    M+=df[i];
    D+=lambda[i];
  }

  double a[K];
  double prod=exp(-D)*pow(beta,M);
  for(int i=0;i<nt;i++){
    prod*=pow(coeffs[i],-df[i]);
  }
  prod=sqrt(prod);
  
  a[0]=prod; 
  if(a[0]<FLT_MIN){
     stats->error_num=GX2CDF_UNDERFLOW;
     return 0;
  }

  double s=a[0];
  double p=a[0]*gsl_cdf_chisq_P(x/beta, M); 
  stats->chi2_calls++; 

  double g[K-1];
  for(int k=0;k<=K-2;k++){
    double chi2=gsl_cdf_chisq_P(x/beta, M+2*(k+1));
    stats->chi2_calls++;
    double err_=(1-s)*chi2;
    stats->truncation_error=err_;
    if(fabs(err_)<tol){
      stats->error_num=GX2_CONVERGED;
      return p;
    }
    g[k]=0;
    for(int i=0;i<nt;i++){
      g[k]+=df[i]*pow((1-beta/coeffs[i]),k+1);
      g[k]+=beta*(k+1)*pow(1-beta/coeffs[i],k)*(lambda[i]/coeffs[i]);
    }
    a[k+1]=0;
    for(int u=0;u<=k;u++){
      a[k+1]+=g[u]*a[k-u];
    }
    a[k+1]/=2*(k+1);
    s+=a[k+1];
    p+=a[k+1]*chi2;
  }
  return 0;
}

typedef struct {
  int nt;
  double *coeffs;
  int *df;
  double *lambda;
  int *error_num;
  int *chi2_calls;
  double p;
} p_t;

static double f(double x, void *args){
  gx2_stats_t stats;
  double tol=1e-17;
  p_t *ps;
  ps=(p_t *)(args);
  double ret=gx2cdf(ps->nt, x, ps->coeffs, ps->df, ps->lambda, tol, &stats)-ps->p;
  *(ps->chi2_calls)+=stats.chi2_calls;
  if(*(ps->error_num)==GX2_CONVERGED){
    *(ps->error_num)=stats.error_num;
  }
  return ret;
}

double gx2ppf(int nt, double p, double *coeffs, int *df, double *lambda, double tol, gx2_stats_t *stats){
  p_t args;
  int error_num;
  int chi2_calls;

  if(p<=0 || p>=1){
    stats->error_num=GX2PPF_NOT_A_PROBABILITY;
  }
  
  stats->truncation_error=0;  /* not used */
  
  args.nt=nt;
  args.coeffs=coeffs;
  args.df=df;
  args.lambda=lambda;
  args.error_num=&error_num;
  error_num=GX2_CONVERGED;
  args.chi2_calls=&chi2_calls;
  chi2_calls=0;
  args.p=p;
  stats_t stats_={0,0,0};
  double rb=10;
  double x0;
  for(double rb_=rb;rb_<1500;rb_*=2){
    stats_.error_num=0;
    x0=brentq(f,0,rb_,tol,tol,K,&stats_,&args);
    if(stats_.error_num!=SIGNERR){
      break;
    }
  }
  if(*(args.error_num)!= GX2_CONVERGED){
    stats->error_num=*(args.error_num);
    return 0;
  }else{
    if(stats_.error_num==SIGNERR){
      stats->error_num=GX2PPF_SIGN_ERROR;
      return 0;
    }else if(stats_.error_num==CONVERR){
      stats->error_num=GX2PPF_NOT_CONVERGED;
      return 0;
    }else{
      stats->error_num=GX2_CONVERGED;
    }
  }
  stats->chi2_calls+=*(args.chi2_calls);
  stats->iterations=stats_.iterations;
  stats->funcalls=stats_.funcalls;
  return x0;
}

void gx2_stats_disp(gx2_stats_t *stats){
  printf("error_num        =%d\n",stats->error_num);
  printf("iterations       =%d\n",stats->iterations);
  printf("funcalls         =%d\n",stats->funcalls);
  printf("chi2_calls       =%d\n",stats->chi2_calls);
  printf("truncation_error =%f\n",stats->truncation_error);
}
