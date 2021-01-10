
#ifndef _LEMANOI_H
#define _LEMANOI_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 
#include "levmar.h"

#define HAVE_LAPACK
#define DBL_RAND_MAX (double)(RAND_MAX)
#define M_PI   3.14159265358979323846

//So tham so hoi quy
#define NUMPARAM  4

#define	_X1	p[0]
#define	_X2	p[1]
#define	_X3	p[2]
#define	_X4	p[3]

//So tham so bang du lieu
#define NUMCOL	  12

// So dong bang du lieu 
#define NUMLINE   12
#define NUMROW	  12

double gdata[NUMLINE][NUMCOL];

#define	_T		gdata[i][0]
#define	_Y		gdata[i][1]

double gNoise(double, double);	

#endif _LEMANOI_H
