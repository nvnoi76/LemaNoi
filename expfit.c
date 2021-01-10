/*
MunozGroup1-paper
Exanple 7 - Page 11
By Nguyen Van Noi
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 
#define HAVE_LAPACK

#include "levmar.h"

#define DBL_RAND_MAX (double)(RAND_MAX)


#define M_PI   3.14159265358979323846


/* Gaussian noise with mean m and variance s */
double gNoise(double m, double s)
{
	double r1, r2, val;
	
	r1=((double)rand())/DBL_RAND_MAX;
	
	r2=((double)rand())/DBL_RAND_MAX;

	val=sqrt(-2.0*log(r1))*cos(2.0*M_PI*r2);
	
	val=s*val+m;
	
	return val;
}

/* model to be fitted to measurements: x_i = p[0]*exp(-p[1]*i) + p[2], i=0...n-1 */
/*void expfunc(double *p, double *x, int m, int n, void *data)
{
	register int i;
	
	for(i=0; i<n; ++i)
	{
		x[i]=p[0]*exp(-p[1]*i) + p[2];
	}
}*/

/*void expfunc(double *p, double *x, int m, int n, void *data)
{
	register int i;
	
	for(i=0; i<n; ++i)
	{
		x[i]=p[0]*exp(p[1]*(i+1));
	}
}
*/

void expfunc(double *p, double *x, int m, int n, void *data)
{
	register int i;
	
	for(i=0; i<n; ++i)
	{
		x[i]=p[0]* sin(p[1]*(i+1)+p[2]) + p[3]; //p[0]*exp(p[1]*(i+1));
	}
}


/* Jacobian of expfunc() */
/*void jacexpfunc(double *p, double *jac, int m, int n, void *data)
{   
	register int i, j;	

	for(i=j=0; i<n; ++i){
		jac[j++]=exp(p[1]*(i+1));
		jac[j++]=p[0]*(i+1)*exp(p[1]*(i+1));
	}
}*/

void jacexpfunc(double *p, double *jac, int m, int n, void *data)
{   
	register int i, j;
	
	/* fill Jacobian row by row */
	for(i=j=0; i<n; ++i){
		jac[j++]= sin(p[1]*(i+1)+p[2]);    
		jac[j++]= (i+1)*p[0]* cos(p[1]*(i+1)+p[2]);  
		jac[j++]= p[0]* cos(p[1]*(i+1)+p[2]);
		jac[j++]= 1;
	}
}




#define n 12
#define m 4

int main2()
{
	//const int n=40, m=3; // 40 measurements, 3 parameters
	double p[m], x[n], opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	double population[]={61, 65, 72, 78, 85, 90, 92, 92, 88, 81, 72, 63};
	register int i;
	int ret;
	
	/* generate some measurement using the exponential model with
	* parameters (5.0, 0.1, 1.0), corrupted with zero-mean
	* Gaussian noise of s=0.1
	*/
	//  INIT_RANDOM(0);
	
	srand(time(NULL));
	
	for(i=0; i<n; ++i)
		x[i]= population[i];// (5.0*exp(-0.1*i) + 1.0) + gNoise(0.0, 0.1);
	
	/* initial parameters estimate: (1.0, 0.0, 0.0) */
	//p[0]=1.0; p[1]=0.0; p[2]=0.0;
	
	//p[0]=6.0; p[1]=0.3;

	p[0]=17.0; p[1]=0.5, 	p[2]=10.5; p[3]=77.0;

	/* optimization control parameters; passing to levmar NULL instead of opts reverts to defaults */
	opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
	opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used 
	
	/* invoke the optimization function */
	ret=dlevmar_der(expfunc, jacexpfunc, p, x, m, n, 1000, opts, info, NULL, NULL, NULL); // with analytic Jacobian
	//ret=dlevmar_dif(expfunc, p, x, m, n, 1000, opts, info, NULL, NULL, NULL); // without Jacobian
	printf("Levenberg-Marquardt returned in %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]);
	printf("Best fit parameters: %.7g %.7g %.7g %.7g \n", p[0], p[1], p[2], p[3]);
	
	exit(0);
}
