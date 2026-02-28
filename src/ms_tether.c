
#include "mex.h"
#include <math.h>

#include "ms_misc.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	if ( nlhs > 0 || nrhs != 8 ){
		mexErrMsgTxt("Input error for tether");
		plhs[0] = mxCreateDoubleScalar(0.0);
	}

	
	double *f = mxGetPr(prhs[0]);
	double *r = mxGetPr(prhs[1]);
	double *rl = mxGetPr(prhs[2]);
	char *ptype = (char *)mxGetData(prhs[3]);
	char *types = (char*)mxGetData(prhs[4]);
	double kspring = mxGetScalar(prhs[5]);
	double *lbox = mxGetPr(prhs[6]);
	unsigned int npart = (unsigned int)mxGetScalar(prhs[7]);

	double dr[3];	
	for ( unsigned n=0; n<npart; n++ ){
		if ( types[n] == ptype[0] ){
      		for ( unsigned k=0; k<3; k++ ){
				unsigned idx = k*npart + n;
				dr[k] = rl[idx] - r[idx];
				_Wrap( dr[k], lbox[k] );

				f[idx] += kspring * dr[k];
      		}
    	}
	}


}
