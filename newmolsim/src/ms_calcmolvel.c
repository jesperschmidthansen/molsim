#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

		
	if ( nlhs > 1 || nrhs != 5 ){
		mexErrMsgTxt("Input error for calcmolvel");
	}


	double *v = mxGetPr(prhs[0]);
	double *amass = mxGetPr(prhs[1]);
	unsigned int npart = (unsigned int)mxGetScalar(prhs[2]);
	double *atom_idx = mxGetPr(prhs[3]);
	unsigned int nuau = (unsigned)mxGetScalar(prhs[4]);

	unsigned nmols = npart/nuau;

	plhs[0] = mxCreateDoubleMatrix(nmols, 3, mxREAL);
	double *ptr = (double *)mxGetPr(plhs[0]);

	for ( unsigned n=0; n<nmols; n++){
        
		for ( int k=0; k<3; k++ ) ptr[k*nmols + n] = 0.0;
        double mass = 0.0;

		for ( unsigned i=0; i<nuau; i++ ){
            int aidx = (int)atom_idx[i*nmols+n]-1;
		   	mass += amass[aidx];	
            for ( int k=0; k<3; k++ ){
				ptr[k*nmols + n] += amass[aidx]*v[aidx];
			}
		}

		for ( int k=0; k<3; k++ ) ptr[k*nmols + n] /= mass;
	}

}
