#include "mex.h"
#include <stdio.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

		
	if ( nlhs > 2 || nrhs != 7 ){
		printf("%d %d\n", nlhs, nrhs);
		mexErrMsgTxt("Input error for calcmolpos");
	}


	double *r = mxGetPr(prhs[0]);
	double *amass = mxGetPr(prhs[1]);
	unsigned int npart = (unsigned int)mxGetScalar(prhs[2]);
	double *atom_idx = mxGetPr(prhs[3]);
	unsigned int nuau = (unsigned)mxGetScalar(prhs[4]);
	int *cross = (int *)mxGetPr(prhs[5]);
	double *lbox = mxGetPr(prhs[6]);

	unsigned nmols = npart/nuau;

	plhs[0] = mxCreateDoubleMatrix(nmols, 3, mxREAL);
	double *ptr = (double *)mxGetPr(plhs[0]);

	plhs[1] = mxCreateDoubleMatrix(nmols, 1, mxREAL);
	double *ptr_m = (double *)mxGetPr(plhs[1]);

	for ( unsigned n=0; n<nmols; n++){
        
		for ( int k=0; k<3; k++ ) ptr[k*nmols + n] = 0.0;
        double mass = 0.0;

		for ( unsigned i=0; i<nuau; i++ ){
            unsigned int aidx = (unsigned int)atom_idx[i*nmols+n]-1;
			if ( aidx > npart-1 ){
				mexErrMsgTxt("Boom");
			}
		   	mass += amass[aidx];	
            for ( int k=0; k<3; k++ ){
				double rtrue =  r[k*npart+aidx] + cross[k*npart+aidx]*lbox[k]; 
				ptr[k*nmols + n] += amass[aidx]*rtrue;
			}
		}

		for ( int k=0; k<3; k++ ) ptr[k*nmols + n] /= mass;

		ptr_m[n] = mass;
	}

}
