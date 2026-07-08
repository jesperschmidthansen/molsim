
#include "mex.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>

#include "ms_misc.h"

#define HELPTXT "Usage: [epot Pconf] = ms_ljexcl(force, exclude type, params, pos, all types, neighblist, lbox)"

void  _lj_neighb(double *epot, double *force, double *pconf, const double *pos, const char *ptypes, const double *param, 
					const double *lbox, const char *types, const int *neighb_list, const unsigned npart); 


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	if ( nlhs > 2 || nrhs != 7 ) mexErrMsgTxt(HELPTXT);

	double epot = 0.0f;
	
	double *f = mxGetPr(prhs[0]);
	char *excludetype = (char *)mxGetData(prhs[1]);
	double *params = mxGetPr(prhs[2]);
	double *r = mxGetPr(prhs[3]);
	char *types = (char*)mxGetData(prhs[4]);
	int *neighb_list = (int*)mxGetData(prhs[5]);
	double *lbox = mxGetPr(prhs[6]);

	unsigned int npart = mxGetM(prhs[0]);
	
	plhs[1] = mxCreateDoubleMatrix(3,3, mxREAL);
	double *ptr = (double *)mxGetPr(plhs[1]);
	for ( int k=0; k<3; k++)
		for ( int kk=0; kk<3; kk++ ) ptr[k*3+kk]=0.0; 

	_lj_neighb(&epot, f, ptr, r, excludetype, params, lbox, types, neighb_list, npart);

	const double ivol = 1.0/(lbox[0]*lbox[1]*lbox[2]);
	for ( int k=0; k<3; k++)
		for ( int kk=0; kk<3; kk++ ) ptr[k*3+kk] = ptr[k*3+kk]*ivol; 

	plhs[0] = mxCreateDoubleScalar(epot);

}


// ptype length 1 - exclude
// types length npart - particles' type
void  _lj_neighb(double *epot, double *force, double *pconf, const double *pos, const char *excludetype, const double *param, 
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

		if ( types[i1] == *excludetype ) continue;
		
		n = 0;
		while ( 1 ) {

			i2 = neighb_list[n*npart + i1];
	
			if ( i2 == -1 ) break; 
			if ( types[i2] == *excludetype ){ n++; continue; }

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
				
			n++;
		}
	}


}

