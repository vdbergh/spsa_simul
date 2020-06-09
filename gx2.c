#include <math.h>
#include <stdio.h>
#include <float.h>
#include <gsl/gsl_cdf.h>
#include <assert.h>

#include "gx2.h"
#include "brentq.h"

static const int K=1000;

static void kahan_sum(double *sum, double *c, double term){
  double y=term-*c;
  double t=*sum+y;
  *c=(t-*sum)-y;
  *sum=t;
}

void gx2_validate(int nt, double *coeffs, int *df, double *lambda, gx2_stats_t *stats){
  if(nt<=0){
    stats->error_num=GX2_NEGATIVE_TERMS;
    return;
  }
  for(int i=0;i<nt;i++){
    if(coeffs[i]<=0){
       stats->error_num=GX2_NEGATIVE_COEFFS;
       return;
    }
    if(df[i]<=0){
       stats->error_num=GX2_NEGATIVE_DF;
       return;
    }
    if(lambda[i]<0){
       stats->error_num=GX2_NEGATIVE_LAMBDA;
       return;
    }
  }
}

static void gx2stats(int nt, double *coeffs, int *df, double *lambda, double *mu, double *sigma2){
  *mu=0;
  *sigma2=0;
  for(int i=0;i<nt;i++){
    *mu+=coeffs[i]*(df[i]+lambda[i]);
    *sigma2+=coeffs[i]*coeffs[i]*(2*df[i]+4*lambda[i]);
  }
}

double gx2cdf(int nt, double x, double *coeffs, int *df, double *lambda, gx2_stats_t *stats){
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
    cx2cdf stats_t stats;

    double p=gx2cdf(nt,x,coeffs,df,lambda,&stats);

    Inputs:
    nt        number of terms in the weighted sum
    x         point at which to evaluate the CDF
    coeffs    array of coefficients of the non-central chi-squares
    df        array of degrees of freedom of the non-central chi-squares
    lambda    array of non-centrality parameters (sum of squares of
              means) of the non-central chi-squares
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

  stats->iterations=0; /* not used */
  stats->chi2_calls=0;
  stats->truncation_error=0;
  stats->funcalls=0; /* not used */
  stats->error_num=0;

  gx2_validate(nt,coeffs,df,lambda,stats);
  if(stats->error_num!=0){
    return 0;
  }

  if(x<=0){
    return 0;
  }

  if(x==INFINITY){
    return 1;
  }

  stats->error_num=GX2CDF_NOT_CONVERGED;

  double coeffs_min=-1;
  for(int i=0;i<nt;i++){
    if(coeffs_min<0 || coeffs[i]<coeffs_min){
      coeffs_min=coeffs[i];
    }
  }

  double beta=coeffs_min;

  double M=0;
  double D=0;
  for(int i=0;i<nt;i++){
    M+=df[i];
    D+=lambda[i];
  }

  double prod=exp(-D)*pow(beta,M);
  for(int i=0;i<nt;i++){
    prod*=pow(coeffs[i],-df[i]);
  }
  prod=sqrt(prod);
  
  if(prod<DBL_MIN){
     stats->error_num=GX2CDF_UNDERFLOW;
     return 0;
  }

  double a[K];
  a[0]=prod; 
  double p=a[0]*gsl_cdf_chisq_P(x/beta, M); 
  stats->chi2_calls++; 

  double sum_a=a[0];
  double g[K-1];
  double  c_p=0;
  for(int k=0;k<=K-2;k++){
    g[k]=0;
    for(int i=0;i<nt;i++){
      g[k]+=df[i]*pow((1-beta/coeffs[i]),k+1);
      g[k]+=beta*(k+1)*pow(1-beta/coeffs[i],k)*(lambda[i]/coeffs[i]);
    }
    a[k+1]=0;
    double c_a=0;
    for(int u=0;u<=k;u++){
      kahan_sum(&(a[k+1]),&c_a,g[u]*a[k-u]);
    }
    a[k+1]/=2*(k+1);
    sum_a+=a[k+1];
    double chi2=gsl_cdf_chisq_P(x/beta, M+2*(k+1));
    stats->chi2_calls++;
    stats->truncation_error=(1-sum_a)*chi2;
    double p_old=p;
    kahan_sum(&p,&c_p,a[k+1]*chi2);
    if(p==p_old){
      stats->error_num=GX2_CONVERGED;
      break;
    }
  }
  return p;
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
  p_t *ps;
  ps=(p_t *)(args);
  double ret=gx2cdf(ps->nt, x, ps->coeffs, ps->df, ps->lambda, &stats)-ps->p;
  *(ps->chi2_calls)+=stats.chi2_calls;
  if(*(ps->error_num)==GX2_CONVERGED){
    *(ps->error_num)=stats.error_num;
  }
  return ret;
}

double gx2ppf(int nt, double p, double *coeffs, int *df, double *lambda, gx2_stats_t *stats){
  /*
    Returns the PPF of a generalized chi-squared (a weighted sum of
    non-central chi-squares with positive weights) by numerically inverting the CDF.

    Example:
    double nt      =3;
    double p       =0.1481396848655497;
    double coeffs[]={1,5,2};
    int    df[]    ={1,2,3};
    double lambda[] ={2,3,7};
    cx2cdf stats_t stats;

    double x=gx2ppf(nt,p,coeffs,df,lambda,&stats);

    Inputs:
    nt        number of terms in the weighted sum
    p         point at which to evaluate the PPF
    coeffs    array of coefficients of the non-central chi-squares
    df        array of degrees of freedom of the non-central chi-squares
    lambda    array of non-centrality parameters (sum of squares of
              means) of the non-central chi-squares
    Outputs:
    x         computed PPF
    stats     statistcs (convergence, chi2_calls, funcalls, iterations)
  */

  p_t args;
  int error_num;
  int chi2_calls;

  stats->chi2_calls=0;
  stats->iterations=0;
  stats->funcalls=0;
  stats->truncation_error=0;  /* not used */
  stats->error_num=0;

  gx2_validate(nt,coeffs,df,lambda,stats);
  if(stats->error_num!=0){
    return 0;
  }
  
  if(p<0 || p>1){
    stats->error_num=GX2PPF_NOT_A_PROBABILITY;
    return 0;
  }

  if(p==0){
    return 0;
  }else if(p==1){
    return INFINITY;
  }

  args.nt=nt;
  args.coeffs=coeffs;
  args.df=df;
  args.lambda=lambda;
  args.error_num=&error_num;
  error_num=GX2_CONVERGED;
  args.chi2_calls=&chi2_calls;
  chi2_calls=0;
  args.p=p;

  stats_t brentq_stats={0,0,0};
  double mu;
  double sigma2;
  double sigma;
  gx2stats(nt,coeffs,df,lambda,&mu,&sigma2);
  sigma=sqrt(sigma2);
  double x0;
  double rb_=mu;
  while(1){
    x0=brentq(f,0,rb_,0,2*DBL_EPSILON,K,&brentq_stats,&args);
    if((*(args.error_num)!= GX2_CONVERGED) || (brentq_stats.error_num!=SIGNERR)){
      break;
    }
    rb_+=sigma;
  }
  if(*(args.error_num)!= GX2_CONVERGED){
    stats->error_num=*(args.error_num);
  }else if(brentq_stats.error_num==SIGNERR){
      stats->error_num=GX2PPF_SIGN_ERROR;
  }else if(brentq_stats.error_num==CONVERR){
    stats->error_num=GX2PPF_NOT_CONVERGED;
  }else{
    stats->error_num=GX2_CONVERGED;
  }
  stats->chi2_calls+=*(args.chi2_calls);
  stats->iterations=brentq_stats.iterations;
  stats->funcalls=brentq_stats.funcalls;
  return x0;
}

void gx2_stats_disp(gx2_stats_t *stats){
  printf("error_num        =%d\n",stats->error_num);
  printf("iterations       =%d\n",stats->iterations);
  printf("funcalls         =%d\n",stats->funcalls);
  printf("chi2_calls       =%d\n",stats->chi2_calls);
  printf("truncation_error =%f\n",stats->truncation_error);
}
