#include "mex.h"

#define _Periodic( x, y )                 \
 {                                           \
 if ( x > y ) x -= y;                        \
 else if  ( x < 0.0f ) x += y;                  \
 }


void _leapfrog(double *ekin, double *v, double *r, double *f, const double lbox[3], const unsigned npart, const double dt);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	
	double *v = mxGetPr(prhs[0]);
	double *r = mxGetPr(prhs[1]);
	double *f = mxGetPr(prhs[2]);
	double *lbox = mxGetPr(prhs[3]);
	unsigned int npart = (unsigned int)mxGetScalar(prhs[4]);
	double dt = mxGetScalar(prhs[5]);

	double ekin = 0.0f;

	_leapfrog(&ekin, v, r, f, lbox, npart, dt);

	plhs[0] = mxCreateDoubleScalar(ekin);
}


void _leapfrog(double *ekin, double *v, double *r, double *f, const double lbox[3], const unsigned npart, const double dt){
	
	// NOTE MASS UNITY ONLY
	for ( unsigned n=0; n<npart; n++ ){
    	for ( int k=0; k<3; k++ ){
    		unsigned idx = k*npart + n;  
      		v[idx] += f[idx]*dt;
      		r[idx] += v[idx]*dt;
     		
		   	_Periodic(r[idx], lbox[k]); 

     	 	double vhalf = v[idx] - 0.5*f[idx]*dt;
      		*ekin = *ekin + 0.5*vhalf*vhalf;
    	}
	}

}


