#include "mex.h"
#include <omp.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	if ( nlhs > 0 || nrhs != 1 ){
		mexErrMsgTxt("Input error for setomp");
		plhs[0] = mxCreateDoubleScalar(0.0);
	}

	unsigned int nthreads = (unsigned int)mxGetScalar(prhs[0]);

	omp_set_num_threads(nthreads);
}
	

