#include "mex.h"

#define _Periodic( cross, x, y )                  \
 {                                                \
 	cross = 0;                                    \
	if ( x > y ) { x -= y; cross = 1; }           \
 	else if  ( x < 0.0f ) { x += y; cross = -1;}  \
 }


void _leapfrog(double *ekin, double *Pkin, double *v, double *r, double *f, double *mass, 
		int *cross, const double lbox[3], const unsigned npart, const double dt);
	

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	
	if ( nlhs > 2 || nrhs != 8 )
		mexErrMsgTxt("Input error for leapfrog");

	double *v = mxGetPr(prhs[0]);
	double *r = mxGetPr(prhs[1]);
	double *f = mxGetPr(prhs[2]);
	double *m = mxGetPr(prhs[3]);
	int *cross = (int *)mxGetPr(prhs[4]);
	double *lbox = mxGetPr(prhs[5]);
	unsigned int npart = (unsigned int)mxGetScalar(prhs[6]);
	double dt = mxGetScalar(prhs[7]);

	double ekin = 0.0f;
	plhs[1] = mxCreateDoubleMatrix(3,3, mxREAL);
	double *ptr = (double *)mxGetPr(plhs[1]);
	for ( int k=0; k<3; k++)
			for ( int kk=0; kk<3; kk++ ) ptr[k*3+kk] = 0.0; 

	_leapfrog(&ekin, ptr, v, r, f, m, cross, lbox, npart, dt);

	const double ivol = 1.0/(lbox[0]*lbox[1]*lbox[2]);
	for ( int k=0; k<3; k++)
		for ( int kk=0; kk<3; kk++ ) ptr[k*3+kk] = ptr[k*3+kk]*ivol; 

	plhs[0] = mxCreateDoubleScalar(ekin);
	
}


void _leapfrog(double *ekin, double *Pkin, double *v, double *r, double *f, double *mass, 
		int *cross, const double lbox[3], const unsigned npart, const double dt){
	
	int cross_idx[3]; double vhalf[3];
	for ( unsigned n=0; n<npart; n++ ){
		double invmass = 1.0/mass[n];
		for ( int k=0; k<3; k++ ){
    		unsigned idx = k*npart + n;  
      		v[idx] += f[idx]*invmass*dt;
      		r[idx] += v[idx]*dt;
     		
		   	_Periodic(cross_idx[k], r[idx], lbox[k]); 
			cross[idx] += cross_idx[k];

     	 	vhalf[k] = v[idx] - 0.5*f[idx]*dt;
      		*ekin = *ekin + 0.5*vhalf[k]*vhalf[k]*mass[n];
    	}
		
		for ( int k=0; k<3; k++ )
			for ( int kk=0; kk<3; kk++ ) 
				Pkin[3*k+kk] += vhalf[k]*vhalf[kk]*mass[n];

	}

}


