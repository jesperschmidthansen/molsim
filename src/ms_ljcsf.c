
#include "mex.h"
#include <math.h>
#include <omp.h>
#include <float.h>

#include "ms_misc.h"

void  _ljsf(double *epot, double *force, double *pconf, const double *pos, const double *z, const char *ptypes, const double *param, 
					const double *lbox, const char *types, const int *neighb_list, const unsigned npart); 


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	if ( nlhs > 2 || nrhs != 9 )
		mexErrMsgTxt("Input error for ljcsf");

	double epot = 0.0f;
	
	double *f = mxGetPr(prhs[0]);
	char *ptypes = (char *)mxGetData(prhs[1]);
	double *params = mxGetPr(prhs[2]);
	double *r = mxGetPr(prhs[3]);
	double *z = mxGetPr(prhs[4]);
	char *types = (char*)mxGetData(prhs[5]);
	int *neighb_list = (int*)mxGetData(prhs[6]);
	double *lbox = mxGetPr(prhs[7]);
	unsigned int npart = (unsigned int)mxGetScalar(prhs[8]);
	
	plhs[1] = mxCreateDoubleMatrix(3,3, mxREAL);
	double *ptr = (double *)mxGetPr(plhs[1]);
	for ( int k=0; k<3; k++)
		for ( int kk=0; kk<3; kk++ ) ptr[k*3+kk]=0.0; 

	_ljsf(&epot, f, ptr, r, z, ptypes, params, lbox, types, neighb_list, npart);

	const double ivol = 1.0/(lbox[0]*lbox[1]*lbox[2]);
	for ( int k=0; k<3; k++)
		for ( int kk=0; kk<3; kk++ ) ptr[k*3+kk] = ptr[k*3+kk]*ivol; 

	plhs[0] = mxCreateDoubleScalar(epot);
}

// ptypes length 2 - interaction
// types length npart - particles' type
void  _ljsf(double *epot, double *force, double *pconf, const double *pos, const double *charges, 
		const char *ptypes, const double *param, const double *lbox, const char *types, const int *neighb_list, const unsigned npart) {
	int i1, i2, n, k, kk;
	double r2, ft, f[3], r[3], rri, rri3;
	const double cf = param[0], eps=param[1], sigma=param[2], aw=param[3];
	const double cf_sf = param[4];
	const double shift = 4.0*eps*(pow(sigma/cf, 12.) - aw*pow(sigma/cf,6.));
	const double cf2 = cf*cf; 
	const double cf2_sf = cf_sf*cf_sf; const double icf2_sf = 1.0/cf2_sf;
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

				if ( fabs(charges[i1]) > DBL_EPSILON && r2 < cf2_sf ){
					double zizj = charges[i2]*charges[i1]; 
					double rij = sqrt(r2);
					ft = zizj*(1.0/r2 - icf2_sf)/rij; 
	
					for ( k=0; k<3; k++ ){
	  					f[k] = ft*r[k];
	  					force[i1 + k*npart] += f[k];
	  					force[i2 + k*npart] += -f[k];
					}		

					*epot +=  zizj*(1.0/rij + (rij-cf_sf)*icf2_sf - 1.0/cf_sf);
					for ( k=0; k<3; k++ )
	  					for ( kk=0; kk<3; kk++ )
	    					pconf[k*3+kk] += f[k]*r[kk];

				}

			}
			n++;
		}
	}

}

