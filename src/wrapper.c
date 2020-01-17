#include <R.h>
#include <Rmath.h>

/******************/
/* random numbers */
/******************/
 
/* simulation */
void F77_SUB(rndstart)(void) 
{ GetRNGstate(); }

void F77_SUB(rndend)(void) 
{ PutRNGstate(); }

void F77_SUB(ggrmultinom)(int n, double* prob, int K, int* rN)
{ rmultinom(n, prob, K, rN); }

double F77_SUB(ggrnorm)(double *mu, double *sigma) 
{ return rnorm(*mu, *sigma); }

double F77_SUB(ggrexp)(double *scale)
{return rexp(*scale);}

double F77_SUB(ggrgam)(double *a, double *scale)
{return rgamma(*a, *scale);}

double F77_SUB(ggrbet)(double *a, double *b)
{return rbeta(*a, *b);}

double F77_SUB(ggrunif)(double *a, double *b)
{return runif(*a, *b);}

double F77_SUB(ggrbinom)(double *n, double *p)
{return rbinom(*n, *p);}

double F77_SUB(ggrpois)(double *lambda)
{return rpois(*lambda);}


/* density */
/* normal df */
double F77_SUB(ggdnorm)(double *x, double *mu, double *sigma, int *give_log)
{return dnorm( *x, *mu, *sigma, *give_log) ; }

/* gamma df */
double F77_SUB(ggdgamma)(double *x, double *shape, double *scale, int *give_log)
{return dgamma( *x, *shape, *scale, *give_log) ; }



/*********************/
/* special functions */
/*********************/

/* gamma function  */
double F77_SUB(gggamfn)(double *x) 
{ return gammafn(*x); }

/* log of gamma function */
double F77_SUB(gglgamfn)(double *x) 
{ return lgammafn(*x); }

/* Bessel K  */
double F77_SUB(ggbesselk)(double *x, double *nu, double *expo)
{return bessel_k(*x,*nu,*expo);}
 

/* normal cdf */
double F77_SUB(ggpnorm)(double *x, double *mu, double *sigma,
			int *lower_tail, int *give_log)
{return pnorm( *x, *mu, *sigma, *lower_tail, *give_log); }


/* log(1+x) for |x| << 1  */
double F77_SUB(gglog1p)(double *x) 
{ return log1p(*x); }


/*********************/
/* Writing to files  */
/*********************/
/* void F77_SUB(ggwint1)(char *filename, int *output)  */
/* { */
/*   FILE *pFile; */
/*   pFile = fopen(filename,"a");  */
/*   fprintf(pFile, "%d \n", *output); */
/*   fclose(pFile); */
/* } */

/* void F77_SUB(ggwdble1)(char *filename, double *output)  */
/* { */
/*   FILE *pFile; */
/*   pFile = fopen(filename,"a");  */
/*   fprintf(pFile, "%f \n", *output); */
/*   fclose(pFile); */
/* } */
