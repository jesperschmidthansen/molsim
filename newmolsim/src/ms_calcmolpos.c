#include "mex.h"
#include <stdio.h>
#include <math.h>
#include "ms_misc.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

		
	if ( nlhs > 2 || nrhs != 7 ){
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

	double pos[100][3]; 

	for ( unsigned n=0; n<nmols; n++){
        
		unsigned int aidx = (unsigned int)atom_idx[n]-1;

		for ( int k=0; k<3; k++ ) pos[0][k] = r[k*npart+aidx]; 
		double mass = amass[aidx];

		for ( unsigned i=1; i<nuau; i++ ){
			if ( aidx > npart-1 ){
				mexErrMsgTxt("Atomic index exceeds the number of atoms");
			}
		

            for ( int k=0; k<3; k++ ){
				double dr =  r[k*npart+aidx] - r[k*npart+aidx-1]; 
				_Wrap(dr, lbox[k]);				
				pos[i][k] = pos[i-1][k] + dr;
			}		
			mass += amass[aidx];
		}

		for ( int k=0; k<3; k++ ){
			int idx = k*nmols + n; 
			ptr[idx] = 0.0;
		}
		
		for ( unsigned i=0; i<nuau; i++ ){
			aidx = (unsigned int)atom_idx[i*nmols+n]-1;
			for ( int k=0; k<3; k++ ){
				int idx = k*nmols + n; 
				ptr[idx] += amass[aidx]*pos[i][k];
			}
		}
					
		for ( int k=0; k<3; k++ ){
			int idx = k*nmols + n;
			ptr[idx] = ptr[idx]/mass;
			_Periodic0( ptr[idx], lbox[k] );
		}

		ptr_m[n] = mass;

	}

}
