/*
LAMANOI by Nguyen Van Noi
Date : 11/01/2015
Modified : 13/07/2017
Version: 1.0.1
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 
#include "levmar.h"

#define HAVE_LAPACK
#define DBL_RAND_MAX (double)(RAND_MAX)
#define M_PI   3.14159265358979323846

//So tham so hoi quy
#define NUMPARAM  3

#define	_X0	p[0]
#define	_Y0	p[1]
#define	_R	p[2]


//So tham so bang du lieu
#define NUMCOL	  2

// So dong bang du lieu 
#define NUMLINE   82
#define NUMROW	  82

double gdata[NUMLINE][NUMCOL];

#define	_X		gdata[i][0]
#define	_Y		gdata[i][1]

double sqr(double x);
int ReadData(int n, char * Filename);
//Ham f(x)
void fxFunc(double *p, double *x, int m, int n, void *data)
{
	register int i;	
	for(i=0; i<n; ++i)
		x[i]=  sqr(_X-_X0)+sqr(_Y-_Y0) -sqr(_R);
}
void JacfxFunc(double *p, double *jac, int m, int n, void *data)
{   
	register int i, j;	
	for(i=j=0; i<n; ++i)
	{				
		jac[j++]=2*(_X0-_X);
		jac[j++]=2*(_Y0-_Y);
		jac[j++]=(-2)*_R;							
	}
}
void fxFunc1(double *p, double *x, int m, int n, void *data)
{
	register int i;	
	for(i=0; i<n; ++i)
		x[i]=  sqrt(sqr(_X-_X0)+sqr(_Y-_Y0)) -_R;
}
// Ma tran Jacobi cua f(x)
void JacfxFunc1(double *p, double *jac, int m, int n, void *data)
{   
	register int i, j;	
	for(i=j=0; i<n; ++i)
	{				
		jac[j++]=(_X0-_X)/sqrt(sqr(_X-_X0)+sqr(_Y-_Y0));
		jac[j++]=(_Y0-_Y)/sqrt(sqr(_X-_X0)+sqr(_Y-_Y0));
		jac[j++]=-1;							
	}
}
void main()
{
	const int n= NUMROW, m= NUMPARAM;
	double  x[NUMROW], opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	int ret;
	int i;		
	//double p[NUMPARAM]={5.0, 6.5, 21};
	double p[NUMPARAM]={13, 15.5, 35};
	char * filein = "datavd2.txt";	
	if (!ReadData(NUMLINE,filein))
		printf("Loi doc file %s!\n", filein), exit(-1);
	for(i=0;i<NUMROW;i++)
		x[i] = 0;
	opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; 
	opts[3]=1E-20;	opts[4]=LM_DIFF_DELTA; 
	ret=dlevmar_der(fxFunc, JacfxFunc, p, x, m, n, 10000, opts, info, NULL, NULL, NULL); 				
	printf("Finished %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]); 
	printf("X=[%.14g, %.14g, %.14g]T\n",_X0,_Y0,_R);
}
int ReadData(int n, char * Filename)
{
	register int i;	
	FILE * file = fopen(Filename,"rt");
	if (!file) 
		return 0;	
	for(i=0; i<n; ++i)
	{
		fscanf(file,"%lf",&_X);
		fscanf(file,"%lf",&_Y);
	}	
	fclose(file);
	
	file = fopen(".\\lemaNOIout1.csv","wt");
	if (!file)
		return 0;
	
	for(i=0; i<n; ++i)
	{
		fprintf(file," %.14g,",_X);
		fprintf(file," %.14g\n",_Y);
	}
	fclose(file);
	return 1;	
}
double sqr(double x)
{
	return x*x;
}