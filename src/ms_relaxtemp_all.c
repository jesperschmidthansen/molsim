#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <unistd.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	
	if ( nlhs > 0 || nrhs != 5 ){
		mexErrMsgTxt("Input error for relaxtemp");
		plhs[0] = mxCreateDoubleScalar(0.0);
	}

	double *v = mxGetPr(prhs[0]);
	unsigned npart = (unsigned)mxGetScalar(prhs[1]); 
	double *mass = mxGetPr(prhs[2]);
	double tau = mxGetScalar(prhs[3]);
	double T0 = mxGetScalar(prhs[4]);

	const double rational = 1.0/3.0;
	double total_momentum[3]={0.0f};
	double total_mass=0.0f;

	for ( unsigned n=0; n<npart; n++ ){
		double Tkin_particle = 0.0;
		for ( int k=0; k<3; k++ ){
			unsigned idx = k*npart + n;
			Tkin_particle += v[idx]*v[idx];
		}
		Tkin_particle *= rational*mass[n];
	
		double fac = sqrt( 1.0 + tau*(T0/Tkin_particle - 1.0) );
		printf("%f\n", Tkin_particle);
		for ( int k=0; k<3; k++ ){
			unsigned idx = k*npart + n;
			v[idx] *= fac;
			total_momentum[k] += v[idx]*mass[n];
		}
		total_mass += mass[n];
	}


	// Abusing total_momentum variable
	for ( int k=0; k<3; k++ ) 
		total_momentum[k] = total_momentum[k]/total_mass;

	for ( unsigned n=0; n<npart; n++ ){
		for ( int k=0; k<3; k++ ){
			unsigned idx = k*npart + n;
			v[idx] = v[idx] - total_momentum[k];	
		}
	}

}


