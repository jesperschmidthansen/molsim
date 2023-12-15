
#ifndef __SEPCUDAMISC_H__
#define __SEPCUDAMISC_H__

#include "sepcudadefs.h"
#include "sepcudamol.h"


void sep_cuda_mem_error(void);
void sep_cuda_file_error(void);

void sep_cuda_copy(sepcupart *ptr, char opt_quantity, char opt_direction);
void sep_cuda_copy_energies(sepcusys *sptr);

void sep_cuda_load_lattice_positions(sepcupart *ptr, const char *xyzfile);

void sep_cuda_save_crossings(sepcupart *ptr, const char *filestr, float time);
void sep_cuda_save_xyz(sepcupart *ptr, const char *filestr);

void sep_cuda_compress_box(sepcupart *aptr, float rho0, float compressfactor[3]);

float sep_cuda_eval_momentum(float momentum[3], sepcupart *aptr);

void sep_cuda_reset_momentum(sepcupart *aptr);

bool sep_cuda_logrem(unsigned n, int base);

void sep_cuda_set_molforcecalc_on(sepcusys *sptr, unsigned isample);

float sep_cuda_wrap_host(float x, float lbox);
float sep_cuda_periodic_host(float x, float lbox);
float sep_cuda_dot_host(float3 a);

/* Kernel functions */
__global__ void sep_cuda_set_value(float *variable, float value);
__global__ void sep_cuda_set_value(int *variable, int value);
__global__ void sep_cuda_set_prevpos(float4 *p, float4 *pprev, unsigned npart);

__global__ void sep_cuda_print_value(float *value);

__global__ void sep_cuda_sum_distance(float *totalsum, float *dist, unsigned npart);
__global__ void sep_cuda_sum_energies(float3 *totalsum, float4* dx, float4 *dv, float4 *df, 
									 float dt, float *epot, unsigned npart);
__global__ void sep_cuda_sum_ekin(float3 *totalsum, const char type, float4* dx, float4 *dv, float4 *df, 
									 float dt, float *epot, unsigned npart);

__global__ void sep_cuda_reset_variable(float3 *a);
__global__ void sep_cuda_reset_iteration(float4 *force, float *epot, float4 *press, float4 *sumpress, float3 *energies, unsigned npart);
__global__ void sep_cuda_reset_mol_fij(float3 *force, unsigned nmol);

__global__ void sep_cuda_calc_dist(float *dist, float4 *p, float4 *pprev, float3 lbox, unsigned npart);
__global__ void sep_cuda_get_pressure(float4 *press, float4 *pos, float4 *vel, float4 *ppress, int npart);

/* Wrapper interface */
void sep_cuda_reset_iteration(sepcupart *pptr);
void sep_cuda_get_energies(sepcupart *ptr);
void sep_cuda_get_pressure(double *npress, double *shearpress, sepcupart *aptr);

#endif
