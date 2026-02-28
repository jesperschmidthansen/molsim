#include "mex.h"
#include "ms_misc.h"
#include <math.h>
 
double _ryckbell(double *force, double *r, unsigned ndihedrals, int ttype, double *dtypes, double *pidx, 
		double *lbox, double *coeffs, unsigned npart);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	if ( nlhs > 1 || nrhs != 9 )
		mexErrMsgTxt("Input error for Ryckaert-Bellemann");

	double *f = mxGetPr(prhs[0]);
	double *r = mxGetPr(prhs[1]);
	unsigned int ndihedrals = (unsigned int)mxGetScalar(prhs[2]);
	int this_type = (int)mxGetScalar(prhs[3]);
	double *dtypes = mxGetData(prhs[4]);
	double *coeffs = mxGetPr(prhs[5]);
	double *pidx = mxGetPr(prhs[6]);
	double *lbox = mxGetPr(prhs[7]);
	unsigned int npart = (unsigned int)mxGetScalar(prhs[8]);

	
	double epot = _ryckbell(f, r, ndihedrals, this_type, dtypes, pidx, lbox, coeffs, npart);

	plhs[0] = mxCreateDoubleScalar(epot);
}


double _ryckbell(double *force, double *r, unsigned ndihedrals, int ttype, double *dtypes, double *pidx, 
		double *lbox, double *coeffs, unsigned npart){
	double dr1[3], dr2[3], dr3[3];
	double c11, c12, c13, c22, c23, c33, g[6];

 
	double epot = 0.0f;
	for ( unsigned n=0; n<ndihedrals; n++ ){

		if ( (int)dtypes[n] == ttype ){

			for ( int k=0; k<6; k++ ) g[k] = coeffs[n + k*ndihedrals];
			
			unsigned a = (unsigned)pidx[n] - 1;
			unsigned b = (unsigned)pidx[ndihedrals + n] - 1;
			unsigned c = (unsigned)pidx[2*ndihedrals + n] - 1;
			unsigned d = (unsigned)pidx[3*ndihedrals + n] - 1;

			for ( int k=0; k<3; k++ ){
				dr1[k] = r[k*npart + b] - r[k*npart + a]; _Wrap( dr1[k], lbox[k] );
				dr2[k] = r[k*npart + c] - r[k*npart + b]; _Wrap( dr2[k], lbox[k] );
				dr3[k] = r[k*npart + d] - r[k*npart + c]; _Wrap( dr3[k], lbox[k] );
			}

			_Dot3(c11, dr1, dr1); _Dot3(c12, dr1, dr2); _Dot3(c13, dr1, dr3);
			_Dot3(c22, dr2, dr2); _Dot3(c23, dr2, dr3); _Dot3(c33, dr3, dr3);
			

			double cA = c13*c22 - c12*c23;
			double cB1 = c11*c22 - c12*c12;
			double cB2 = c22*c33 - c23*c23;

			// To avoid the special problem of aligned bonds 
			// (can happen for som initial conditions)
			if ( fabs(cB1) > DIHEDRAL_THRESHOLD && fabs(cB2) > DIHEDRAL_THRESHOLD ){
				
				double cD = sqrt(cB1*cB2); 	double cc = cA/cD;

				double f = -(g[1]+(2.*g[2]+(3.*g[3]+(4.*g[4]+5.*g[5]*cc)*cc)*cc)*cc);
				double t1 = cA; 
				double t2 = c11*c23 - c12*c13;
				double t3 = -cB1; 
				double t4 = cB2;
				double t5 = c13*c23 - c12*c33; 
				double t6 = -cA;
				double cR1 = c12/c22; 
				double cR2 = c23/c22;

				for ( int k=0; k<3; k++ ){
					double f1 = f*c22*(t1*dr1[k] + t2*dr2[k] + t3*dr3[k])/(cD*cB1);
					double f2 = f*c22*(t4*dr1[k] + t5*dr2[k] + t6*dr3[k])/(cD*cB2);

					force[k*npart + a] += f1;
					force[k*npart + b] += - (1.0 + cR1)*f1 + cR2*f2;
					force[k*npart + c] += cR1*f1 - (1.0 + cR2)*f2;
					force[k*npart + d] += f2;
				}

				epot += g[0]+(g[1]+(g[2]+(g[3]+(g[4]+g[5]*cc)*cc)*cc)*cc)*cc;
			}
		}

	}

	return epot;
}


