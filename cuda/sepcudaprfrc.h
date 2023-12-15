
#ifndef __SEPCUDAPRFRC_H__
#define  __SEPCUDAPRFRC_H__

#include "sepcudadefs.h"
#include "sepcudamisc.h"
#include <stdlib.h>

void sep_cuda_check_neighblist(sepcupart *ptr, float skin);
void sep_cuda_set_exclusion(sepcupart *aptr, const char rule[]);

// Kernels
__global__ 
void sep_cuda_build_neighblist(int *neighlist, float4 *p, float *dist, float cf, float3 lbox, unsigned nneighmax, unsigned npart);

__global__ 
void sep_cuda_build_neighblist(int *neighlist, float *dist, float4 *p, int *molindex, float cf, float3 lbox, unsigned nneighmax, unsigned npart); 

/* Pair interactions - types specified */
__global__ 
void sep_cuda_lj(const char type1, const char type2, float4 params, int *neighblist, float4 *pos, float4 *force,
							int *molindex, float *epot, float4 *press, unsigned maxneighb, float3 lbox, const unsigned npart);

/* Pair interactions - all types have same interactions (faster) */
__global__ 
void sep_cuda_lj(int *neighblist, float4 *pos, float4 *force,
		float *epot, float4 *press, unsigned maxneighb, float3 lbox, const unsigned npart);

/* Shifted force */
__global__ 
void sep_cuda_lj_sf(const char type1, const char type2, float3 params, int *neighblist, float4 *pos, float4 *force,
								float *epot, float4 *press, unsigned maxneighb, float3 lbox, const unsigned npart);
/* Coulomb shifted force */
__global__
void sep_cuda_sf(float cf, int *neighblist, float4 *pos, float4 *vel, float4 *force,
							float *epot, float4 *press, unsigned maxneighb, float3 lbox, const unsigned npart);

/* Calculate the molecular force due to LJ interactions */
__global__ void sep_cuda_calc_molforce(float3 *mforce,  const char type1, const char type2, float3 params, float4 *pos, 
										int *neighblist,  unsigned maxneighb, float4 *force, 
										float3 lbox, int *molindex, unsigned nmol, const unsigned npart);

/* Calculate molecular force due to sf Coulomb interactions */
__global__ void sep_cuda_calc_molforce(float3 *mforce, float cf, int *neighblist, float4 *pos, float4 *vel, 
										 unsigned maxneighb, int *molindex, unsigned nmols, 
										 float3 lbox, const unsigned npart);

/* Lattice forces */
__global__ 
void sep_cuda_lattice_force(const char type, float springConstant, float4 *pos, float4 *pos0, float4 *force,
									   float3 lbox, const unsigned npart);


/* Wrapper interfaces*/

void sep_cuda_force_lj(sepcupart *pptr, const char types[], float params[4]);
void sep_cuda_force_lj(sepcupart *pptr);
void sep_cuda_force_lj_sf(sepcupart *pptr, const char types[], float params[3]);
void sep_cuda_force_sf(sepcupart *pptr, const float cf);
void sep_cuda_update_neighblist(sepcupart *pptr, float maxcutoff);
void sep_cuda_force_lattice(sepcupart *pptr, const char type, float springConstant);

#endif
