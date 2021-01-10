/*
LAMANOI by Nguyen Van Noi
Date : 11/01/2015
Modified : 13/07/2017
Version: 1.0.1
*/
#include "lemaNoi.h"
int ReadData(int n, char * Filename);
//Ham f(x)
void fxFunc(double *p, double *x, int m, int n, void *data)
{
	register int i;	
	for(i=0; i<n; ++i)
		x[i]=  _X1* sin(_X2* _T + _X3) + _X4; //x1*sin(x2*t + x3) + x4
}
// Ma tran Jacobi cua f(x)
void JacfxFunc(double *p, double *jac, int m, int n, void *data)
{   
	register int i, j;	
	for(i=j=0; i<n; ++i)
	{		
		
		jac[j++]=sin(_X2*_T + _X3); 		//df/dx1=sin(x2*t + x3)
		jac[j++]=_T*_X1*cos(_X2*_T + _X3);  //df/dx2=t*x1*cos(x2*t + x3)
		jac[j++]=_X1*cos(_X2*_T + _X3);		//df/dx3=x1*cos(x2*t + x3)
		jac[j++]=1;							//df/dx4=1
	}
}
void main()
{
	const int n= NUMROW, m= NUMPARAM;
	double  x[NUMROW], opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	int ret;
	int i;		
	double p[NUMPARAM]={15.5, 0, 12, 78.25};
	char * filein = "nhietdo1.txt";	
	if (!ReadData(NUMLINE,filein))
		printf("Loi doc file %s!\n", filein), exit(-1);
	for(i=0;i<NUMROW;i++)
		x[i] = _Y;
	opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; 
	opts[3]=1E-20;	opts[4]=LM_DIFF_DELTA; 
	ret=dlevmar_der(fxFunc, JacfxFunc, p, x, m, n, 10000, opts, info, NULL, NULL, NULL); 				
	printf("Finished %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]); 
	printf("X=[%.14g, %.14g, %.14g, %.14g]T\n",_X1,_X2,_X3,_X4);
}
int ReadData(int n, char * Filename)
{
	register int i;	
	FILE * file = fopen(Filename,"rt");
	if (!file) 
		return 0;	
	for(i=0; i<n; ++i)
	{
		fscanf(file,"%lf",&_T);
		fscanf(file,"%lf",&_Y);
	}	
	fclose(file);
	
	file = fopen(".\\lemaNOIout.csv","wt");
	if (!file)
		return 0;
	
	for(i=0; i<n; ++i)
	{
		fprintf(file," %.14g,",_T);
		fprintf(file," %.14g\n",_Y);
	}
	fclose(file);
	return 1;	
}