#include "mex.h"

#define _Wrap( x, y )                          \
{                                                \
if ( x > 0.5*y ) x -= y;                         \
else if  ( x < -0.5*y ) x += y;                  \
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

		
	if ( nlhs > 1 || nrhs != 5 ){
		mexErrMsgTxt("Input error for calcdr2");
	}

	
	double *r = mxGetPr(prhs[0]);
	double *r0 = mxGetPr(prhs[1]);
	int *cross = (int *)mxGetPr(prhs[2]);
	double *lbox = mxGetPr(prhs[3]);
	unsigned int npart = (unsigned int)mxGetScalar(prhs[4]);

	double maxdr2 = 0.0f;
	for ( unsigned n=0; n<npart; n++ ){
		
		double dr2 = 0.0f; 
		for ( unsigned k=0; k<3; k++ ){
			unsigned idx = k*npart + n;  
			double xtrue = r[idx] + lbox[k]*cross[idx];
			double dr = xtrue - r0[idx];
			_Wrap( dr, lbox[k] );
			dr2 += dr*dr;
		}
		if ( dr2 > maxdr2 ) maxdr2 = dr2;
	}	

	plhs[0] = mxCreateDoubleScalar(maxdr2);
}


