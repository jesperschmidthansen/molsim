#include "mex.h"

#define _Periodic( cross, x, y )                  \
 {                                                \
 	cross = 0;                                    \
	if ( x > y ) { x -= y; cross = 1; }           \
 	else if  ( x < 0.0f ) { x += y; cross = -1;}  \
 }


void _leapfrog(double *ekin, double *v, double *r, double *f, int *cross, const double lbox[3], const unsigned npart, const double dt);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	
	double *v = mxGetPr(prhs[0]);
	double *r = mxGetPr(prhs[1]);
	double *f = mxGetPr(prhs[2]);
	int *cross = (int *)mxGetPr(prhs[3]);
	double *lbox = mxGetPr(prhs[4]);
	unsigned int npart = (unsigned int)mxGetScalar(prhs[5]);
	double dt = mxGetScalar(prhs[6]);

	double ekin = 0.0f;

	_leapfrog(&ekin, v, r, f, cross, lbox, npart, dt);

	plhs[0] = mxCreateDoubleScalar(ekin);
}


void _leapfrog(double *ekin, double *v, double *r, double *f, int *cross, const double lbox[3], const unsigned npart, const double dt){
	
	// NOTE_ THIS VRESION SUPPORTS MASS UNITY ONLY
	int cross_idx[3];
	for ( unsigned n=0; n<npart; n++ ){
		for ( int k=0; k<3; k++ ){
    		unsigned idx = k*npart + n;  
      		v[idx] += f[idx]*dt;
      		r[idx] += v[idx]*dt;
     		
		   	_Periodic(cross_idx[k], r[idx], lbox[k]); 
			cross[idx] += cross_idx[k];

     	 	double vhalf = v[idx] - 0.5*f[idx]*dt;
      		*ekin = *ekin + 0.5*vhalf*vhalf;
    	}
	}

}


