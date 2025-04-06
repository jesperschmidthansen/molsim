#include "mex.h"
#include <math.h>
#include <omp.h>
#include <float.h>

#include "ms_misc.h"


void _sf(double *epot, double *pconf, double *force, double *pos, double *charges, 
		const int *neighb_list, double *lbox, unsigned int npart, double cf);
	
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	if ( nlhs > 2 || nrhs != 7 )
		mexErrMsgTxt("Input error for sf");

	double *f = mxGetPr(prhs[0]);
	double *r = mxGetPr(prhs[1]);
	double *z = mxGetPr(prhs[2]);
	int *neighb_list = (int*)mxGetData(prhs[3]);
	double *lbox = mxGetPr(prhs[4]);
	unsigned int npart = (unsigned int)mxGetScalar(prhs[5]);
	double cf = mxGetScalar(prhs[6]);

	plhs[1] = mxCreateDoubleMatrix(3,3, mxREAL);
	double *ptr = (double *)mxGetPr(plhs[1]);
	for ( int k=0; k<3; k++)
			for ( int kk=0; kk<3; kk++ ) ptr[k*3+kk]=0.0; 
	
	double epot;	
	_sf(&epot, ptr, f, r, z, neighb_list, lbox, npart, cf);

	const double ivol = 1.0/(lbox[0]*lbox[1]*lbox[2]);
		for ( int k=0; k<3; k++)
			for ( int kk=0; kk<3; kk++ ) ptr[k*3+kk] = ptr[k*3+kk]*ivol; 

	plhs[0] = mxCreateDoubleScalar(epot);

}

void _sf(double *epot, double *pconf, double *force, double *pos, double *charges, 
		const int *neighb_list, double *lbox, unsigned int npart, double cf){
	int i1, i2, n, k, kk;
  	double r[3], r2, ft, f[3], ecoul, rij, zizj;
  	const double cf2 = cf*cf, icf2=1.0/cf2, icf=1.0/cf;
  	size_t lvec = npart*3;

  	ecoul = 0.0;
  
#pragma omp parallel for schedule(dynamic) private(i1, n, i2, k, kk, r, r2, ft, f)			\
  reduction(+:ecoul, force[:lvec], pconf[:9]) 
	for ( i1=0; i1<npart; i1++ ){
    
		if ( fabs(charges[i1]) < DBL_EPSILON ) continue;
		
		n = 0;
		while (1){
			i2 = neighb_list[n*npart + i1];
		  	
			if ( i2 == -1 ) break; 
		  
		  	for ( k=0; k<3; k++ ){
				r[k] = pos[k*npart + i1] - pos[k*npart + i2];
				_Wrap( r[k], lbox[k] );
		  	}
		  
		  	r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
		  
		  	if ( r2 < cf2 ){	
				zizj = charges[i2]*charges[i1]; 
				rij = sqrt(r2);
				ft = zizj*(1.0/r2 - icf2)/rij; 
	
				for ( k=0; k<3; k++ ){
	  				f[k] = ft*r[k];
	  				force[i1 + k*npart] += f[k];
	  				force[i2 + k*npart] += -f[k];
				}		

				ecoul +=  zizj*(1.0/rij + (rij-cf)*icf2 - icf);
				for ( k=0; k<3; k++ )
	  				for ( kk=0; kk<3; kk++ )
	    				pconf[k*3+kk] += f[k]*r[kk];
      		}		
      		n++;
    	}
  	}

  	*epot = ecoul;

}	
