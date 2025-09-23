#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	
	if ( nlhs > 1 || nrhs != 5 ){
		mexErrMsgTxt("Input error for nose_hoover_all");
	}

	double *f = mxGetPr(prhs[0]);
	double xi = mxGetScalar(prhs[1]);
	double *v = mxGetPr(prhs[2]);
	unsigned npart = (unsigned)mxGetScalar(prhs[3]); 
	double *mass = mxGetPr(prhs[4]);


	double ekin = 0.0f;	
	for ( unsigned n=0; n<npart; n++ ){
		for ( int k=0; k<3; k++ ){
			f[k*npart + n] += -xi*mass[n]* v[k*npart + n];	
			ekin += mass[n]*v[k*npart + n]*v[k*npart + n];	
		}

	}
	
	plhs[0] = mxCreateDoubleScalar(0.5*ekin);

}

