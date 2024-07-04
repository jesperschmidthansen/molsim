
#include "mex.h"
#include <math.h>
#include <omp.h>

#define _Wrap( x, y )                          \
{                                                \
if ( x > 0.5*y ) x -= y;                         \
else if  ( x < -0.5*y ) x += y;                  \
}



void _lj(double *epot, double *f, double *r, const double cf, const double lbox[3], const unsigned npart);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	
	double *f = mxGetPr(prhs[0]);
	double *r = mxGetPr(prhs[1]);
	double *lbox = mxGetPr(prhs[2]);
	unsigned int npart = (unsigned int)mxGetScalar(prhs[3]);

	double epot = 0.0f;

	_lj(&epot, f, r, 2.5, lbox, npart);

	plhs[0] = mxCreateDoubleScalar(epot);
}


void _lj(double *epot, double *f, double *r, const double cf, const double lbox[3], const unsigned npart){
	double dr[3], fij[3];

	const double shift = 4.0*(pow(1.0/cf, 12.) - pow(1.0/cf,6.));
	const double cf2 = cf*cf;

	for ( unsigned n=0; n<npart-1; n++ ){
		for ( unsigned m=n+1; m<npart; m++ ){
			
			double dr2 = 0.0;
			for ( unsigned k=0; k<3; k++ ){
				dr[k] = r[k*npart + n] - r[k*npart + m];
				_Wrap( dr[k], lbox[k] );
				dr2 += dr[k]*dr[k];
			}

			if ( dr2 < cf2 ){ 

				double rri = 1.0/dr2; double rri3 = rri*rri*rri;
				double ft = 48.0*rri3*(rri3 - 0.5)*rri;
				
				for ( unsigned k=0; k<3; k++ ){
					fij[k] = ft*dr[k];
					f[k*npart + n] += fij[k];
					f[k*npart + m] -= fij[k];
				}
									 
				*epot = *epot + 4.0*rri3*(rri3 - 1.0) - shift; 
			}

		}
	}

}	

