#include "mex.h"
#include "ms_misc.h"
#include <math.h>


double _harmonic(double *f, double *r, unsigned int nbonds, int ttype, double *btype, double *ks, 
					double *lbond, double *pidx, double *lbox, unsigned npart);
	

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	if ( nlhs > 1 || nrhs != 10 )
		mexErrMsgTxt("Input error for harmonic");

	double *f = mxGetPr(prhs[0]);
	double *r = mxGetPr(prhs[1]);
	unsigned int nbonds = (unsigned int)mxGetScalar(prhs[2]);
	int this_type = (int)mxGetScalar(prhs[3]);
	double *btypes = mxGetData(prhs[4]);
	double *ksprings = mxGetPr(prhs[5]);
	double *lbonds = mxGetPr(prhs[6]);
	double *pidx = mxGetPr(prhs[7]);
	double *lbox = mxGetPr(prhs[8]);
	unsigned int npart = (unsigned int)mxGetScalar(prhs[9]);

	 
	double epot = _harmonic(f, r, nbonds, this_type, btypes, ksprings, lbonds, pidx, lbox, npart);  

	plhs[0] = mxCreateDoubleScalar(epot);
}


double _harmonic(double *f, double *r, unsigned int nbonds, int ttype, double *btype, double *ks, 
											double *lbond, double *pidx, double *lbox, unsigned npart){
	double dr[3];

	double epot = 0.0f;
	for ( unsigned n=0; n<nbonds; n++ ){

		if ( (int)btype[n] == ttype ){
			unsigned a = (unsigned)pidx[n]-1;
			unsigned b = (unsigned)pidx[nbonds + n]-1;
			
			double dr2 = 0.0;
			for ( unsigned k=0; k<3; k++ ){
				dr[k] = r[k*npart + a] - r[k*npart + b];
				_Wrap( dr[k], lbox[k] );
				dr2 += dr[k]*dr[k];
			}

			double dist = sqrt(dr2);
			double ft = -ks[n]*(dist - lbond[n])/dist;

			for ( int k=0; k<3; k++ ){
				f[k*npart+a] += ft*dr[k];
				f[k*npart+b] -= ft*dr[k];
			}

			epot += 0.5*ks[n]*(dist - lbond[n])*(dist - lbond[n]);
		}
	}

	return epot;
}

