#include "mex.h"
#include "ms_misc.h"

#define HELPTXT "Usage [ekin Pkin] = ms_leapfrog(vel, pos, force, mass, crossings, lbox, dt)"

void _leapfrog(double *ekin, double *Pkin, double *v, double *r, double *f, double *mass, 
		int *cross, const double lbox[3], const unsigned npart, const double dt);
	

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	
	if ( nlhs > 2 || nrhs != 7) mexErrMsgTxt(HELPTXT);

	double *v = mxGetPr(prhs[0]);
	double *r = mxGetPr(prhs[1]);
	double *f = mxGetPr(prhs[2]);
	double *m = mxGetPr(prhs[3]);
	int *cross = (int *)mxGetPr(prhs[4]);
	double *lbox = mxGetPr(prhs[5]);
	double dt = mxGetScalar(prhs[6]);

	unsigned int natoms = mxGetM(prhs[0]);

	// Kinetic energy
	double ekin = 0.0f;

	// Kinetic contribution to pressure
	plhs[1] = mxCreateDoubleMatrix(3,3, mxREAL);
	double *Pkin = (double *)mxGetPr(plhs[1]);
	for ( int k=0; k<3; k++)
		for ( int kk=0; kk<3; kk++ ) Pkin[k*3+kk] = 0.0; 

	// Integrate
	_leapfrog(&ekin, Pkin, v, r, f, m, cross, lbox, natoms, dt);

	// Normalize wrt system volume
	const double ivol = 1.0/(lbox[0]*lbox[1]*lbox[2]);
	for ( int k=0; k<3; k++)
		for ( int kk=0; kk<3; kk++ ) Pkin[k*3+kk] = Pkin[k*3+kk]*ivol; 

	plhs[0] = mxCreateDoubleScalar(ekin);
	
}


void _leapfrog(double *ekin, double *Pkin, double *v, double *r, double *f, double *mass, 
		int *cross, const double lbox[3], const unsigned npart, const double dt){
	
	int cross_idx; double vhalf[3];
	for ( unsigned n=0; n<npart; n++ ){
		double invmass = 1.0/mass[n];
		for ( int k=0; k<3; k++ ){
    		unsigned idx = k*npart + n;  
      		v[idx] += f[idx]*invmass*dt;
      		r[idx] += v[idx]*dt;
     		
		   	_Periodic(cross_idx, r[idx], lbox[k]); 
			cross[idx] += cross_idx;

     	 	vhalf[k] = v[idx] - 0.5*f[idx]*dt;
      		*ekin = *ekin + 0.5*vhalf[k]*vhalf[k]*mass[n];
    	}
		
		for ( int k=0; k<3; k++ )
			for ( int kk=0; kk<3; kk++ ) 
				Pkin[3*k+kk] += vhalf[k]*vhalf[kk]*mass[n];

	}

}


