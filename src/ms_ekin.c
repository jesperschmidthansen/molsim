#include "mex.h"
#include <math.h>

#define HELPTXT "Usage: ekin = ms_ekin(vel, particle types, all types, all mass)"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	
	if ( nlhs > 1 || nrhs != 4 ) mexErrMsgTxt("HELPTXT");

	double *v = mxGetPr(prhs[0]);
	char ptype = (char)mxGetScalar(prhs[1]);
	char *types = (char*)mxGetPr(prhs[2]);
	double *mass = mxGetPr(prhs[3]);

	unsigned npart = (unsigned)mxGetM(prhs[0]); 

	double ekin = 0.0f;	
	for ( unsigned n=0; n<npart; n++ ){
		if ( types[n]==ptype ){
			for ( int k=0; k<3; k++ )	
				ekin += mass[n]*v[k*npart + n]*v[k*npart + n];	
		}
	}

	
	plhs[0] = mxCreateDoubleScalar(0.5*ekin);

}

