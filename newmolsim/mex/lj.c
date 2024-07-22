
#include "mex.h"
#include <math.h>
#include <omp.h>

#define _Wrap( x, y )                          \
{                                                \
if ( x > 0.5*y ) x -= y;                         \
else if  ( x < -0.5*y ) x += y;                  \
}



void _lj_brute(double *epot, double *f, double *r, const double cf, const double lbox[3], const unsigned npart);
void  _lj_neighb(double *epot, double *force, double *pconf, const double *pos, const char *ptypes, const double *param, 
					const double *lbox, const char *types, const int *neighb_list, const unsigned npart); 


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	if ( nlhs > 1 || (nrhs != 8 && nrhs != 4) )
		mexErrMsgTxt("Input error for lj");

	double epot = 0.0f;
	
	if ( nrhs == 8 ){
		double *f = mxGetPr(prhs[0]);
		char *ptypes = (char *)mxGetData(prhs[1]);
		double *params = mxGetPr(prhs[2]);
		double *r = mxGetPr(prhs[3]);
		char *types = (char*)mxGetData(prhs[4]);
		int *neighb_list = (int*)mxGetData(prhs[5]);
		double *lbox = mxGetPr(prhs[6]);
		unsigned int npart = (unsigned int)mxGetScalar(prhs[7]);

		double pressure[9] = {0.0f};
		
		_lj_neighb(&epot, f, pressure, r, ptypes, params, lbox, types, neighb_list, npart);
	}
	else if ( nrhs == 4 ){
		double *f = mxGetPr(prhs[0]);
		double *r = mxGetPr(prhs[1]);
		double *lbox = mxGetPr(prhs[2]);
		unsigned int npart = (unsigned int)mxGetScalar(prhs[3]);

		_lj_brute(&epot, f, r, 2.5, lbox, npart);
	}

	plhs[0] = mxCreateDoubleScalar(epot);
}


void _lj_brute(double *epot, double *f, double *r, const double cf, const double lbox[3], const unsigned npart){
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

// ptypes length 2 - interaction
// types length npart - particles' type
void  _lj_neighb(double *epot, double *force, double *pconf, const double *pos, const char *ptypes, const double *param, 
					const double *lbox, const char *types, const int *neighb_list, const unsigned npart) {
	int i1, i2, n, k, kk;
	double r2, ft, f[3], r[3], rri, rri3;
	const double cf = param[0], eps=param[1], sigma=param[2], aw=param[3];
	const double shift = 4.0*eps*(pow(sigma/cf, 12.) - aw*pow(sigma/cf,6.));
	const double cf2 = cf*cf;
	const double eps48 = 48.0*eps;
	const double eps4 = 4.0*eps;
	const double awh = 0.5*aw;
	const double sigmasqr = sigma*sigma;
	const unsigned lvec = 3*npart;


#pragma omp parallel for schedule(dynamic)			\
		private(i1, n, i2, k, kk, r, r2, ft, f, rri, rri3)	\
		reduction(+:epot[:1], force[:lvec], pconf[:9]) 
	for (i1=0; i1<npart; i1++){

		if ( types[i1] != ptypes[0] && types[i1] != ptypes[1] ) continue;

		n = 0;
		while ( 1 ) {

			i2 = neighb_list[n*npart + i1];

			if ( i2 == -1 ) break; 
		
			if ( (types[i1] == ptypes[0] && types[i2] == ptypes[1]) || 
					(types[i1] == ptypes[1] && types[i2] == ptypes[0]) ){
				
				r2 = 0.0;
				for ( k=0; k<3; k++ ){
					r[k] = pos[k*npart + i1] - pos[k*npart + i2];
					_Wrap( r[k], lbox[k] );
					r2 += r[k]*r[k];
				}

				if ( r2 < cf2 ){
					rri = sigmasqr/r2; rri3 = rri*rri*rri;
					ft = eps48*rri3*(rri3 - awh)*rri;

					for ( k=0; k<3; k++ ){
						f[k] = ft*r[k];

						force[i1 + k*npart] += f[k];
						force[i2 + k*npart] += -f[k];
					}

					*epot = *epot + eps4*rri3*(rri3 - aw) - shift; 

					for (k=0; k<3; k++)
						for ( kk=0; kk<3; kk++ ) pconf[k*3+kk] += f[k]*r[kk];
				}
			}
			n++;
		}
	}

}

