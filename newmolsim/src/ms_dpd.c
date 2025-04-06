#include "mex.h"
#include <math.h>
#include <omp.h>

#include "ms_misc.h"

void _dpd_neighb(double *epot, double *force, const double *pos, const double *vel, const char *ptypes, const double *param, 
				const double temperature, const double *lbox, const char *types, const int *neighb_list, const unsigned npart);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	if ( nlhs > 1 || nrhs != 10  )
		mexErrMsgTxt("Input error for lj");

	double epot = 0.0f;
	
	double *f = mxGetPr(prhs[0]);
	double *r = mxGetPr(prhs[1]);
	double *v = mxGetPr(prhs[2]);
	char *ptypes = (char *)mxGetData(prhs[3]); // The two types interacting
	double *params = mxGetPr(prhs[4]); //cf, aij, sigma, dt
	double temperature = mxGetScalar(prhs[5]);
	double *lbox = mxGetPr(prhs[6]);
	char *types = (char*)mxGetData(prhs[7]);
	int *neighb_list = (int*)mxGetData(prhs[8]);
	unsigned int npart = (unsigned int)mxGetScalar(prhs[9]);
		
	_dpd_neighb(&epot, f, r, v, ptypes, params, temperature, lbox, types, neighb_list, npart);
	
	plhs[0] = mxCreateDoubleScalar(epot);

}

void _dpd_neighb(double *epot, double *force, const double *pos, const double *vel, const char *ptypes, const double *param, 
				const double temperature, const double *lbox, const char *types, const int *neighb_list, const unsigned npart){ 
		
	int i1, i2, n, k;
	double r2,  r[3], rhat[3], vij[3], fC[3], fD[3], fR[3], dij, one_dij, randnum, dotrv;
	
	const double cf = param[0];
	const double aij = param[1];
	const double sigma = param[2];
	const double dt = param[3];

	const double cf2 = cf*cf, isqrtdt = 1.0/sqrt(dt);
	const double facchk = 2.0*sqrt(3.0);
	const double gamma = sigma*sigma/(2.0*temperature);
	const size_t lvec = npart*3;
	
#pragma omp parallel for schedule(dynamic)			\
	private(i1, n, i2, k, r, r2, dij, one_dij, dotrv, rhat, vij, randnum, fC, fD, fR) \
	reduction(+:epot[:1], force[:lvec]) 
	for ( i1=0; i1<npart; i1++ ){

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

					dij = sqrt(r2);	one_dij = (1.0-dij);
					
					dotrv = 0.0;
					for ( k=0; k<3; k++ ){
						rhat[k] = r[k]/dij;
						vij[k] = vel[k*npart + i1] - vel[k*npart + i2]; 
						dotrv += rhat[k]*vij[k];
					}

					randnum = (_Rand()-0.5)*facchk;

					for ( k=0; k<3; k++ ){
						// Conservative force
						fC[k] = aij*one_dij*rhat[k];
						// Dissipative force
						fD[k] =-gamma*one_dij*one_dij*dotrv*rhat[k];
						// Random force
						fR[k] = sigma*one_dij*rhat[k]*isqrtdt*randnum; 
						// Summing up
						force[k*npart + i1] += fC[k] + fD[k] + fR[k];
						force[k*npart + i2] -= fC[k] + fD[k] + fR[k];
					}

					// potential energy (conservative force)
					*epot = *epot + 0.5*aij*one_dij*one_dij;
				}  

			}    //  if atom type is considered
			n++;
		}
	}        //  neighbor list loops

}

