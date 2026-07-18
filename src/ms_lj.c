
#include "mex.h"
#include <math.h>
#include <omp.h>

#include "ms_misc.h"

#define HELPTXT "Usage [epot Pconf]= ms_lj(force, itypes, params, pos, atypes, neighbl, lbox)"


void _lj_neighb(double *epot, double *force, double *pconf,
	   			const double *pos, const char *ptypes, const double *param, 
				const double *lbox, const char *types, const int *neighb_list, const unsigned npart); 


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	if ( nlhs > 2 || nrhs != 7 ) mexErrMsgTxt(HELPTXT);

	double *f = mxGetPr(prhs[0]); 
	char *itypes = (char *)mxGetData(prhs[1]);
	double *params = mxGetPr(prhs[2]);
	double *r = mxGetPr(prhs[3]);
	char *atypes = (char*)mxGetData(prhs[4]);
	int *neighbl = (int*)mxGetData(prhs[5]);
	double *lbox = mxGetPr(prhs[6]);

	unsigned int natoms = mxGetM(prhs[0]);

	// Potential energy
	double epot = 0.0f;
	
	// Configurational contribution to presse	
	plhs[1] = mxCreateDoubleMatrix(3,3, mxREAL);
	double *Ppot = (double *) mxGetPr(plhs[1]);
	for ( int k=0; k<3; k++ )
		for ( int kk=0; kk<3; kk++ ) Ppot[k*3+kk] = 0.0; 

	// Calculate forces 
	_lj_neighb(&epot, f, Ppot, r, itypes, params, lbox, atypes, neighbl, natoms);
	
	// Normalizing wrt volume
	const double ivol = 1.0/(lbox[0]*lbox[1]*lbox[2]);
	for ( int k=0; k<3; k++)
		for ( int kk=0; kk<3; kk++ ) Ppot[k*3+kk] = Ppot[k*3+kk]*ivol; 
	
	
	plhs[0] = mxCreateDoubleScalar(epot);

}

void  _lj_neighb(double *epot, double *force, double *pconf, const double *pos, const char *itypes, const double *param, 
				const double *lbox, const char *atypes, const int *neighb_list, const unsigned npart) {
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

#pragma omp parallel for schedule(static)			        \
		private(i1, n, i2, k, kk, r, r2, ft, f, rri, rri3)	\
		reduction(+:epot[:1], force[:lvec], pconf[:9]) 
	for (i1=0; i1<npart; i1++){

		if ( atypes[i1] != itypes[0] && atypes[i1] != itypes[1] ) continue;

		n = 0;
		while ( 1 ) {

			i2 = neighb_list[n*npart + i1];

			if ( i2 == -1 ) break; 
		
			if ( atypes[i2] == itypes[0] || atypes[i2] == itypes[1] ){	
				
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

					for ( k=0; k<3; k++ )
						for ( kk=0; kk<3; kk++ ) pconf[k*3+kk] += f[k]*r[kk];
					
				}

			}
			n++;
		}	
	}
}

/*
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
*/

