
#ifndef __SEPCUDAINTGR_H__
#define __SEPCUDAINTGR_H__

#include "sepcudadefs.h"
#include "sepcudamisc.h"
	
__global__ 
void sep_cuda_leapfrog(float4 *pos, float4 *vel, float4 *force, float *dist, int3 *crossing, float dt, float3 lbox, unsigned npart);

__global__ 
void sep_cuda_update_nosehoover(float *alpha, float3 *denergies, float temp0, float tau, float dt, unsigned int npart);

__global__
void sep_cuda_nosehoover(float *alpha, float4 *pos, float4 *vel, float4 *force, unsigned npart);

__global__ 
void sep_cuda_nosehoover(float *alpha, const char type, float4 *pos, float4 *vel, float4 *force, unsigned npart);

/* Wrapper interfaces */
void sep_cuda_thermostat_nh(sepcupart *pptr, float temp0, float tau);
void sep_cuda_thermostat_nh(const char type, sepcupart *pptr, float temp0, float tau);

void sep_cuda_integrate_leapfrog(sepcupart *pptr);

#endif
