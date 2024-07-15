
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_NNEIGHB 500

#define _Wrap( x, y )                          \
{                                                \
if ( x > 0.5*y ) x -= y;                         \
else if  ( x < -0.5*y ) x += y;                  \
}


void _build_cell_list(int *list, double *pos, double *lbox, unsigned *nsubbox, double *lsubbox, unsigned npart);
void _build_neighb_list(int *nighb_list, double *pos, int *cell_list, unsigned *nsubbox, double cf, double *lbox,  
											double skin, unsigned npart, unsigned max_nneighb);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	// Input check
	if ( nrhs != 7 )
		mexErrMsgTxt ("Input wrong");

	// Get the input variables	
	int *neighb_list = (int *)mxGetPr(prhs[0]);
	double *r = mxGetPr(prhs[1]);
	double *r0 = mxGetPr(prhs[2]);
	double *lbox = mxGetPr(prhs[3]);
	double cf = mxGetScalar(prhs[4]);
	double skin = mxGetScalar(prhs[5]);
	unsigned int npart = (unsigned int)mxGetScalar(prhs[6]);

	// Calculate the cell grid
	unsigned int ncells[3]; double lcells[3];
	for ( int k=0; k<3; k++ ){
  		ncells[k] = (int)(lbox[k]/(cf + skin));
		if ( ncells[k] < 3 )
			fprintf(stderr, "Number of cells in %d direction less than three - BAILING OUT!", ncells[k]);
  		lcells[k] = lbox[k]/ncells[k];
	}

	// Allocate/prepare for cell list build
	int *cell_list = malloc(sizeof(int)*(npart + ncells[0]*ncells[1]*ncells[2]));
	if ( cell_list == NULL ) {
		fprintf(stderr, "Mem. allocation error\n");
		exit(EXIT_FAILURE);
	}

	_build_cell_list(cell_list, r, lbox, ncells, lcells, npart);
	_build_neighb_list(neighb_list, r, cell_list, ncells, cf, lbox, skin, npart, MAX_NNEIGHB);

	for ( unsigned n=0; n<npart; n++ ){
		for ( unsigned k=0; k<3; k++ ){
			unsigned idx = k*npart + n;  
			r0[idx] = r[idx];
		}
	}

	free(cell_list);
}


void _build_cell_list(int *cell_list, double *pos, double *lbox, unsigned *nsubbox, double *lsubbox, unsigned npart){

	const unsigned nsubbox2 = nsubbox[0]*nsubbox[1]; 
	const unsigned nsubbox3 = nsubbox2*nsubbox[2];
	const unsigned length = npart + nsubbox3;

	for ( unsigned n=0; n<length; n++ ) cell_list[n] = -1;

	for ( unsigned n=0; n<npart; n++){
		unsigned i = (unsigned)(pos[n]/lsubbox[0]) + (unsigned)(pos[npart+n]/lsubbox[1])*nsubbox[0] +
			(unsigned)(pos[2*npart+n]/lsubbox[2])*nsubbox2;

		if ( i > length - 1 ) {
			fprintf(stderr, "%s at %d: Index larger than array length", __func__, __LINE__);
			exit(EXIT_FAILURE);
		}

		cell_list[n+nsubbox3] = cell_list[i];
		cell_list[i] = n;  
	}

}


void _build_neighb_list(int *neighb_list, double *pos, int *cell_list, unsigned *nsubbox, double cf, double *lbox,  
											double skin, unsigned npart, unsigned max_nneighb){
	
	const int iofX[] = {0,1,1,0,-1,0,1,1,0,-1,-1,-1, 0, 1}; 
	const int iofY[] = {0,0,1,1, 1,0,0,1,1, 1, 0,-1,-1,-1};
	const int iofZ[] = {0,0,0,0, 0,1,1,1,1, 1, 1, 1, 1, 1}; 
	
	size_t nbytes = sizeof(int)*npart;
	int *index = malloc(nbytes); 
	memset(index, 0, nbytes);

	nbytes = sizeof(int)*(npart+1);
	int *icc = malloc(nbytes); 
	memset(index, 0, nbytes);
	
	unsigned nsubbox2 = nsubbox[1]*nsubbox[0]; 
	unsigned nsubbox3 = nsubbox[2]*nsubbox2;
	
	const double cf2 = (cf + skin)*(cf + skin);

	memset(neighb_list, -1, sizeof(int)*(npart*max_nneighb));

	int j1;
	for ( unsigned m1Z = 0; m1Z < nsubbox[2]; m1Z++ ){
		for ( unsigned m1Y = 0; m1Y < nsubbox[1]; m1Y++ ) {
			for ( unsigned m1X = 0; m1X < nsubbox[0]; m1X++ )  {

				// cell index
				unsigned m1 = m1Z*nsubbox2 + m1Y*nsubbox[0] + m1X;

				// Re-organize
				if ( cell_list[m1] == -1 ) continue;
				else j1 = cell_list[m1];

				int nccell = 0;
				while ( j1 != -1 ) {
					icc[nccell] = j1;
					nccell++;
					j1 = cell_list[j1+nsubbox3];
				}

//#pragma omp parallel for schedule(dynamic) private(offset, j1, m2X, m2Y, m2Z, m2, j2, r2, dr, k)
				for ( int i=0; i<nccell; i++ ){

					j1 = icc[i];

					// Over neighb. cells
					for ( unsigned offset = 0; offset < 14; offset++ ) { 

						int m2X = m1X + iofX[offset];    
						if ( m2X == nsubbox[0] ) m2X = 0;    
						else if ( m2X == -1 ) m2X = nsubbox[0]-1;    

						int m2Y = m1Y + iofY[offset];    
						if ( m2Y == nsubbox[1] )  m2Y = 0;    
						else if ( m2Y == -1 ) m2Y = nsubbox[1]-1;    

						int m2Z = m1Z + iofZ[offset];    
						if ( m2Z == nsubbox[2] ) m2Z = 0;    

						// Neighb. cell index
						int m2 = m2Z*nsubbox2 + m2Y*nsubbox[0] + m2X;

						// Head-index for nieghb. cell
						int j2 = cell_list[m2];

						// Loop over particles in neighb. cell
						while ( j2 != -1 ) {

							if ( m1 != m2 || j2 > j1 ) {
								double dr2 = 0.0f;
								for ( int k = 0; k<3; k++ ){
									double dr = pos[k*npart + j1] - pos[k*npart + j2];
									_Wrap( dr, lbox[k] );
									dr2 += dr*dr;
								}

								if ( dr2 < cf2 ) {
									neighb_list[index[j1]*npart + j1] = j2;
									index[j1]++;
								}	  

								if ( index[j1] == max_nneighb ){
									fprintf(stderr, "Found too many neighbours - Bailing out!");
									exit(EXIT_FAILURE);
								}
							}
							
							// Get next particle in list for cell m2 
							j2 =cell_list[j2+nsubbox3]; 
						}
					} 
				}
			} } }

	free(index);
	free(icc);
}


