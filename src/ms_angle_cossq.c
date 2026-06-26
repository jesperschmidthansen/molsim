#include "mex.h"
#include "ms_misc.h"
#include <math.h>

#define HELPTXT "Usage: epot = ms_angle_cossq(force, pos, angle type, all angle types, spring, zero force angle, part. idxs, lbox)"
	
double _cossq(double *f, double *r, unsigned nangles, int ttype, double *atypes, double *pidx, double *lbox, 
		double *angles0, double *ks, unsigned npart);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	if ( nlhs > 1 || nrhs != 8 ) mexErrMsgTxt(HELPTXT);

	double *f = mxGetPr(prhs[0]);
	double *r = mxGetPr(prhs[1]);
	int this_type = (int)mxGetScalar(prhs[2]);
	double *atypes = mxGetData(prhs[3]);
	double *ksprings = mxGetPr(prhs[4]);
	double *angles = mxGetPr(prhs[5]);
	double *pidx = mxGetPr(prhs[6]);
	double *lbox = mxGetPr(prhs[7]);

	unsigned int npart = mxGetM(prhs[0]);
	unsigned int nangles = mxGetM(prhs[6]);
	 
	double epot = _cossq(f, r, nangles, this_type, atypes, pidx, lbox, angles, ksprings, npart);  

	plhs[0] = mxCreateDoubleScalar(epot);
}


double _cossq(double *force, double *r, unsigned nangles, int ttype, double *atypes, double *pidx, double *lbox, 
		double *angles0, double *ks, unsigned npart){
	double dr1[3], dr2[3];
	double c11, c12, c22;

	double epot = 0.0f;
	for ( unsigned n=0; n<nangles; n++ ){
		if ( (int)atypes[n] == ttype ){
			unsigned a = (unsigned)pidx[n] - 1;
			unsigned b = (unsigned)pidx[nangles + n] - 1;
			unsigned c = (unsigned)pidx[2*nangles + n] - 1;

			for ( unsigned k=0; k<3; k++ ){
				dr1[k] = r[k*npart + b] - r[k*npart + a]; _Wrap( dr1[k], lbox[k] );
				dr2[k] = r[k*npart + c] - r[k*npart + b]; _Wrap( dr2[k], lbox[k] );
			}

			_Dot3(c11, dr1, dr1); _Dot3(c12, dr1, dr2);	_Dot3(c22, dr2, dr2);

			double cD = sqrt(c11*c22);	double angle = M_PI - acos(c12/cD); 

			double dangle = angle-angles0[n];
			double f = -ks[n]*dangle;

			for ( unsigned k=0; k<3; k++ ){
				double f1 = f*((c12/c11)*dr1[k] - dr2[k])/cD;
				double f2 = f*(dr1[k] - (c12/c22)*dr2[k])/cD;

				force[k*npart + a] += f1;
				force[k*npart + b] += (-f1-f2);
				force[k*npart + c] += f2;
			}

			epot += 0.5*ks[n]*dangle*dangle;
		}
	}

	return epot;
}
