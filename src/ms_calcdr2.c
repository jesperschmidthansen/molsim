#include "mex.h"
#include "ms_misc.h"

#define HELPTXT "Usage: maxdrsq = ms_calcdr2(pos, pos0, lbox)"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

		
	if ( nlhs > 1 || nrhs != 3 ) mexErrMsgTxt(HELPTXT);

	double *r = mxGetPr(prhs[0]);
	double *r0 = mxGetPr(prhs[1]);
	double *lbox = mxGetPr(prhs[2]);

	unsigned int npart = mxGetM(prhs[0]);

	double maxdr2 = 0.0f;
	for ( unsigned n=0; n<npart; n++ ){
		
		double dr2 = 0.0f; 
		for ( unsigned k=0; k<3; k++ ){
			unsigned idx = k*npart + n;  
			double dr = r[idx]-r0[idx];	
			_Wrap( dr, lbox[k] );
			dr2 += dr*dr;
		}
		if ( dr2 > maxdr2 ) maxdr2 = dr2;
	}	

	plhs[0] = mxCreateDoubleScalar(maxdr2);
}


