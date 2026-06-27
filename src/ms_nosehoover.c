#include "mex.h"
#include <math.h>

#define HELPTXT "Usage: ekin = ms_nosehoover(force, xi, vel, particle types, all types, all mass)"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	
	if ( nlhs > 1 || nrhs != 6 ) mexErrMsgTxt(HELPTXT);

	double *f = mxGetPr(prhs[0]);
	double xi = mxGetScalar(prhs[1]);
	double *v = mxGetPr(prhs[2]);
	char ptype = (char)mxGetScalar(prhs[3]);
	char *types = (char*)mxGetPr(prhs[4]);
	double *mass = mxGetPr(prhs[5]);

	unsigned npart = (unsigned)mxGetM(prhs[0]); 

	double ekin = 0.0f;	
	for ( unsigned n=0; n<npart; n++ ){
		
		if ( types[n]==ptype ){
			for ( int k=0; k<3; k++ ){
				f[k*npart + n] += -xi*mass[n]* v[k*npart + n];	
				ekin += mass[n]*v[k*npart + n]*v[k*npart + n];	
			}
		}

	}
	
	plhs[0] = mxCreateDoubleScalar(0.5*ekin);

}

