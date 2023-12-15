
#ifndef __SEPCUDAMOL_H__
#define __SEPCUDAMOL_H__

#include "sepcudadefs.h"
#include "sepcudamisc.h"

typedef struct _sepcumol {
	unsigned nmols; 
	unsigned *hnuau, *dnuau;
	
	unsigned nbonds;    /**< Total number of bonds */
	unsigned *hblist, *dblist;    /**< Bond list: (the two bonded part. indicies + bond type)*num_bonds */
	unsigned nbondblocks;
	
	unsigned nangles;    /**< Total number of bonds */
	unsigned *halist, *dalist;    /**< Bond list: (the three bonded part. indicies + bond type)*num_angles */
	unsigned nangleblocks;
	
	unsigned ndihedrals;    /**< Total number of dihedrals */
	unsigned *hdlist, *ddlist;    /**< Bond list: (the four bonded part. indicies + bond type)*num_dihedrals */
	unsigned ndihedralblocks;  

	int *alist;  /** Atomic list - host only */

	float3 *hf, *df;  /** Molecular forces - currently device part disabled */
	float3 *hx, *dx;  /** Centre of mass  - currently device part disabled */
	float3 *hv, *dv;  /** Centre of mass velocity - currently device part disabled */
	float3 *hpel, *dpel; /** Mol. dipole - currently device part disabled */ 

	float *masses;    /** Molecular masses (host only) */
	float3 *hfij, *dfij; /** Config. part of pressure tensor */
} sepcumol;


sepcumol * sep_cuda_init_mol(void);

FILE *sep_cuda_set_file_pointer(FILE *fptr, const char *section);

void sep_cuda_read_bonds(sepcupart *pptr, sepcumol *mptr, const char *file, const char opt);
void sep_cuda_read_angles(sepcumol *mptr, const char *file, const char opt);
void sep_cuda_read_dihedrals(sepcumol *mptr, const char *file, const char opt);

void sep_cuda_free_bonds(sepcumol *mptr);
void sep_cuda_free_angles(sepcumol *mptr);
void sep_cuda_free_mols(sepcumol *mptr);
void sep_cuda_free_dihedrals(sepcumol *mptr);

void sep_cuda_make_atomslist(sepcumol *mptr, sepcupart *pptr);

// Calculators
void sep_cuda_mol_calc_cmprop(sepcupart *pptr, sepcumol *mptr);
void sep_cuda_mol_calc_molpress(double *P, sepcupart *pptr,  sepcumol *mptr);
void sep_cuda_mol_calc_dipoles(sepcupart *pptr, sepcumol *mptr);
double sep_cuda_mol_calc_avdipole(sepcumol *mptr);
void sep_cuda_mol_calc_forceonmol(sepcupart *pptr, sepcumol *mptr);

// Kernels
__global__ void sep_cuda_bond_harmonic(unsigned int *blist, unsigned nbonds, float3 bondspec, 
								  float4 *pos, float4 *force, float3 lbox);
__global__ void sep_cuda_angle(unsigned *alist, unsigned nangles, float3 anglespec, 
								  float4 *pos, float4 *force, float3 lbox);
__global__ void sep_cuda_ryckertbellemann(unsigned *dlist, unsigned ndihedrals, int type, float *params, 
											float4 *pos, float4 *force, float3 lbox);

// Device
__device__ float sep_cuda_mol_dot(float4 a);
__device__ float sep_cuda_mol_dot(float3 a, float3 b);
__device__ float sep_cuda_mol_wrap(float x, float lbox);

// Wrappers
void sep_cuda_force_harmonic(sepcupart *pptr, sepcumol *mptr, int type, float ks, float lbond);
void sep_cuda_force_angle(sepcupart *pptr, sepcumol *mptr, int type, float ktheta, float angle0);
void sep_cuda_force_dihedral(sepcupart *pptr, sepcumol *mptr, int type, float *params);

#endif
