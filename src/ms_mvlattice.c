
#include "mex.h"
#include <math.h>

#include "ms_misc.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	if ( nlhs > 0 || nrhs != 6 ){
		mexErrMsgTxt("Input error for mvlattice");
		plhs[0] = mxCreateDoubleScalar(0.0);
	}

	double *rl = mxGetPr(prhs[0]);
	char *ptype = (char *)mxGetData(prhs[1]);
	double *dr = mxGetPr(prhs[2]);
	char *types = (char*)mxGetData(prhs[3]);
	double *lbox = mxGetPr(prhs[4]);
	unsigned int npart = (unsigned int)mxGetScalar(prhs[5]);

	for ( unsigned n=0; n<npart; n++ ){
		if ( types[n] == ptype[0] ){
      		for ( unsigned k=0; k<3; k++ ){
				unsigned idx = k*npart + n;
				rl[idx] = rl[idx] + dr[k];
				_Periodic0( rl[idx], lbox[k] );
      		}
    	}
	}


}
