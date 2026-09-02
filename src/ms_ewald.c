/*
 * ewald.c
 *
 * Reference Ewald implementation for point charges in an orthorhombic
 * periodic simulation box.
 *
 * Features:
 *   - Separate real-space and reciprocal-space functions
 *   - User-controlled real-space cutoff
 *   - User-controlled number of reciprocal vectors through kmax
 *   - Energy and forces
 *   - Optional exclusion of intramolecular charge interactions
 *
 * Assumptions:
 *   - 3D periodic boundary conditions
 *   - Orthorhombic box
 *   - Electrically neutral simulation cell
 *   - Tin-foil (conducting) boundary conditions
 *
 * Compile:
 *     gcc -O3 -Wall ewald.c -lm
 */

#include "mex.h"
#include <math.h>
#include "ms_misc.h"
#include <float.h>

#define HELPTXT "Usage: epot = ms_ewald(force, positions, charges, neighbour list, lbox, [alpha, cutoff, col. fac], exclusion list)"

double _ewald_short_range(double *force, const double *pos,  const double *charges,  
						const double *lbox,	double alpha,   double cutoff,  double coulomb,  
		 				const int *neighb_list, int npart); 
double _ewald_long_range(double *force, const double *pos, const double *q, const double *lbox,
    					double alpha,  int kmax, double coulomb, int npart);
double _ewald_remove_exclusion(double *force, double *pos, const double *charges, 
								const int *exclusion_list, const int maxexcl, const double coulomb, const double alpha,
								const double *lbox, const int npart);




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	if ( nrhs > 7 || nlhs > 1 ) mexErrMsgTxt(HELPTXT);

	double epot = 0.0f;
	
	double *f = mxGetPr(prhs[0]);
	double *r = mxGetPr(prhs[1]);
	double *z = mxGetPr(prhs[2]);
	int *neighb_list = (int*)mxGetData(prhs[3]);
	double *lbox = mxGetPr(prhs[4]);
	double *params = mxGetPr(prhs[5]);
	int *exclusion_list = (int *)mxGetPr(prhs[6]);

	int npart = mxGetM(prhs[0]);
	int max_exclusion = mxGetN(prhs[6]);

	double alpha = params[0]; 
	double cutoff = params[1]; 
	double coulomb = params[2];

	epot = _ewald_short_range(f, r,  z, lbox, alpha,  cutoff,  coulomb,  neighb_list, npart);
	epot += _ewald_long_range(f, r, z, lbox, alpha, 10 , coulomb, npart);
	epot += _ewald_remove_exclusion(f, r, z,  exclusion_list, max_exclusion, coulomb, alpha, lbox, npart);

	plhs[0] = mxCreateDoubleScalar(epot);

}

/***********************************************************************
 * REAL-SPACE EWALD CONTRIBUTION
 *
 * Potential:
 *
 *     U_real = sum_{i<j} q_i q_j erfc(alpha*r_ij) / r_ij
 *
 * for r_ij < cutoff.
 *
 * Excluded intramolecular pairs are simply omitted here.
 *
 * Parameters
 * ----------
 * n          : number of particles
 * r          : particle positions
 * q          : particle charges
 * molecule   : molecule IDs
 * box        : box dimensions
 * alpha      : Ewald splitting parameter
 * cutoff     : real-space cutoff
 * coulomb    : Coulomb prefactor, e.g. 1 or 138.935456
 * exclude_same_molecule : nonzero => exclude intramolecular interactions
 * force      : forces are ADDED to this array
 *
 * Returns
 * -------
 * Real-space electrostatic energy.
 *
 * Notice that exclusion is made via exclusion list and the neighbour list
 ***********************************************************************/
double _ewald_short_range(double *force, const double *pos,  const double *charges,  
						const double *lbox,	double alpha,   double cutoff,  double coulomb,  
		 				const int *neighb_list, int npart){

    double energy = 0.0;
	int i, j, k, n;
	double r2, r[3], f[3];
    
	const double rc2 = cutoff * cutoff;
    const double alpha2 = alpha * alpha;
    const double two_alpha_over_sqrt_pi = 2.0 * alpha / sqrt(M_PI);

	for ( i=0; i<npart; i++ ){
    
		if ( fabs(charges[i]) < DBL_EPSILON ) continue;
		
		n = 0;
		while (1){
			j = neighb_list[n*npart + i];
		  	
			if ( j == -1 ) break; 

		  	for ( k=0; k<3; k++ ){
				r[k] = pos[k*npart + i] - pos[k*npart + j];
				_Wrap( r[k], lbox[k] );
		  	}
		  
		  	r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];

            if ( r2 <= rc2 ) {

				double rij = sqrt(r2);
				double inv_r = 1.0 / rij;
				double inv_r2 = 1.0 / r2;
				double inv_r3 = inv_r * inv_r2;

				double ar = alpha * rij;

				double erfc_ar = erfc(ar);
				double exp_term = exp(-alpha2 * r2);

				double qq = coulomb * charges[i] * charges[j];

				energy += qq * erfc_ar * inv_r;

				double ft =  qq * (erfc_ar * inv_r3 +two_alpha_over_sqrt_pi *exp_term * inv_r2);

				for ( k=0; k<3; k++ ){
					f[k] = ft*r[k];
					force[i + k*npart] += f[k];
					force[j + k*npart] += -f[k];
				}		
	        }
    		n++;	
		} // end neighbourloop
	} // end particle loop

    return energy;
}


/***********************************************************************
 * RECIPROCAL-SPACE EWALD CONTRIBUTION
 *
 * Reciprocal vectors:
 *
 *   k = 2*pi ( nx/Lx, ny/Ly, nz/Lz )
 *
 * with
 *
 *   -kmax <= nx,ny,nz <= kmax
 *
 * excluding (0,0,0).
 *
 * Thus the number of reciprocal vectors is:
 *
 *   (2*kmax + 1)^3 - 1
 *
 * before any symmetry optimizations.
 *
 *
 * Reciprocal energy:
 *
 * U_k =
 *
 *   (2*pi/V) sum_{k != 0}
 *
 *       exp[-k^2/(4 alpha^2)]
 *       -------------------- |S(k)|^2
 *               k^2
 *
 * where
 *
 *   S(k) = sum_j q_j exp(i k.r_j)
 *
 *
 * The self energy
 *
 *   -alpha/sqrt(pi) sum_i q_i^2
 *
 * is included.
 *
 *
 * INTRAMOLECULAR EXCLUSIONS
 * ------------------------
 *
 * Simply removing atoms from S(k) does NOT correctly implement pair
 * exclusions.
 *
 * If a pair is excluded from the real-space calculation, the
 * corresponding direct-space part still contained in the reciprocal
 * sum is
 *
 *       q_i q_j erf(alpha*r)/r .
 *
 * We therefore subtract this term for every excluded pair. This is the
 * usual Ewald exclusion correction.
 *
 *
 * Returns the reciprocal + self + exclusion correction energy.
 ***********************************************************************/
double _ewald_long_range(double *force, const double *pos, const double *q, const double *lbox,
    					double alpha,  int kmax, double coulomb, int npart){
    const double volume = lbox[0]*lbox[1]*lbox[2];
    const double twopi = 2.0 * M_PI;
    const double fourpi_over_V = 4.0 * M_PI / volume;
    const double alpha2 = alpha * alpha;

    double retval = 0.0;
	double *epot = &retval;
	int i;
	int nx, ny, nz;
	double kx, ky, kz, k2, kr, damping, C, D, s, c, prefactor;

#pragma omp parallel for schedule(static) collapse(3)			        \
		private(nx, ny, nz, kx, ky, kz, k2, kr, damping, C, D,s,c, prefactor)	\
		reduction(+:epot[:1], force[:3*npart]) 
    for ( nx = -kmax; nx <= kmax; ++nx ) {

        kx = twopi * nx / lbox[0];

        for ( ny = -kmax; ny <= kmax; ++ny ) {

            ky = twopi * ny / lbox[1];

            for ( nz = -kmax; nz <= kmax; ++nz ) {

                // k = 0 is excluded.
                if ( nx == 0 && ny == 0 && nz == 0 )  continue;

                kz = twopi * nz / lbox[2];
				k2 = kx * kx + ky * ky + kz * kz;

                //Ewald reciprocal-space damping.
                damping = exp(-k2 / (4.0 * alpha2)) / k2;

                // Calculate the structure factor
                C = 0.0;  D = 0.0;
                for ( i = 0; i < npart; ++i ) {
                    kr = kx*pos[i] + ky*pos[npart + i] + kz*pos[2*npart + i];
                    C += q[i] * cos(kr); D += q[i] * sin(kr);
                }

                // Because both +k and -k are explicitly included, the energy coefficient is 2*pi/V.
                *epot +=  coulomb * (2.0 * M_PI / volume)*damping*(C * C + D * D);

                for ( i = 0; i < npart; ++i ) {
					
					kr = kx*pos[i] + ky*pos[npart + i] + kz*pos[2*npart + i];

                    s = sin(kr); c = cos(kr);

                    prefactor = coulomb*fourpi_over_V*q[i]*damping*(C * s - D * c);

                    force[i] += prefactor * kx;
                    force[npart + i] += prefactor * ky;
                    force[2*npart + i] += prefactor * kz;
                }


            }
        }
    }

    // Remove self part
	double q2sum = 0.0;

    for (int i = 0; i < npart; ++i) q2sum += q[i]*q[i];

    *epot -=  coulomb * alpha / sqrt(M_PI) * q2sum;

	return retval;
}

/*******************************************************************
     * Correction for excluded intramolecular pairs.
     *
     * Since these pairs were omitted from the real-space calculation,
     * subtract
     *
     *     q_i q_j erf(alpha*r)/r
     *
     * from the reciprocal result.
     *******************************************************************/

double _ewald_remove_exclusion(double *force, double *pos, const double *charges, 
								const int *exclusion_list, const int maxexcl, const double coulomb, const double alpha,
								const double *lbox, const int npart){
	double dr[3];
	double epot = 0.0;


	for ( int n=0; n<npart; n++ ){

		for ( int j=0; j<maxexcl; j++ ){
	
			int  m = exclusion_list[n + j*npart];
			
			if ( m == -1 ) break;	
	
			for ( int k=0; k<3; k++ ){
				dr[k] = pos[k*npart + n] - pos[k*npart + m];
				_Wrap( dr[k], lbox[k] );
		  	}

		  	double r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
			double r = sqrt(r2);
			double qq = coulomb*charges[n]*charges[m];

			double ar   = alpha*r;
			double erfc_term = erf(ar);
    		double exp_term  = exp(-ar*ar);
			
		    /* Energy correction */
    		epot -= qq*erfc_term/r;

			double fscalar = qq * (erfc_term/(r2*r)- 2.0*alpha/sqrt(M_PI)*exp_term/r2);
			
			for ( int k=0; k<3; k++ ){
				force[n + k*npart] += fscalar*dr[k];
				force[m + k*npart] -= fscalar*dr[k];
			}

		}
	}

	return epot;
}	


