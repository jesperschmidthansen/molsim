#include "mex.h"
#include <math.h>
#include <omp.h>

#include "ms_misc.h"

void _verlet_dpd(double *ekin, double *r, double *v, double *f, double *pacc, double *pv, 
				const double *mass,  const unsigned npart, const double *lbox, int *cross, 
				const double dt, const unsigned stepnow, const double lambda);

	
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	if ( nlhs > 1 || nrhs != 12 )
		mexErrMsgTxt("Input error for verlet_dpd");

	double ekin = 0.0f;

	double *r = mxGetPr(prhs[0]);
	double *v = mxGetPr(prhs[1]);
	double *f = mxGetPr(prhs[2]);
	double *pa = mxGetPr(prhs[3]);
	double *pv = mxGetPr(prhs[4]);
	double *masses = mxGetPr(prhs[5]);
	int *cross = (int *)mxGetPr(prhs[6]);
	double *lbox = mxGetPr(prhs[7]);
	unsigned int npart = (unsigned int)mxGetScalar(prhs[8]);
	double dt = mxGetScalar(prhs[9]);
	unsigned int istep = (unsigned int)mxGetScalar(prhs[10]);
	double lambda = mxGetScalar(prhs[11]);

	_verlet_dpd(&ekin, r, v, f, pa, pv, masses, npart, lbox, cross, dt, istep, lambda);

	plhs[0] = mxCreateDoubleScalar(ekin);

}

void _verlet_dpd(double *ekin, double *r, double *v, double *f, double *pacc, double *pv, 
				const double *mass,  const unsigned npart, const double *lbox, int *cross, 
				const double dt, const unsigned stepnow, const double lambda){

	int cross_idx;
	for ( unsigned n=0; n<npart; n++ ){
		double invmass = 1.0/mass[n];
		for (int k=0; k<3; k++){
  			unsigned idx = k*npart + n;  
			double acc = f[idx]*invmass;

  			// Previous time step
  			if ( stepnow > 0 )
				v[idx] += 0.5*dt*(acc + pacc[idx]);

  			r[idx] += v[idx]*dt + 0.5*dt*dt*acc;
 			pv[idx] = v[idx] + lambda*dt*acc;
	    	pacc[idx] = acc;

			_Periodic(cross_idx, r[idx], lbox[k]); 
			cross[idx] += cross_idx;
  
  			*ekin = *ekin + 0.5*v[idx]*v[idx]*mass[n];
		} 
	}

}

