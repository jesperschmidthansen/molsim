
#ifndef ___SEPCUDA_H__
#define ___SEPCUDA_H__


#include "cuda.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#define SEP_CUDA_NTHREADS 32  
#define SEP_CUDA_MAXNEIGHBS 512

#define SEP_CUDA_PI 3.14159265

#define SEP_MAX_NUMB_EXCLUSION 20
#define SEP_CUDA_EXCL_NONE 0
#define SEP_CUDA_EXCL_BONDS 1
#define SEP_CUDA_EXCL_MOLECULE 2

#define SEP_CUDA_MAXNUAU 100

typedef struct _sepcupart {
	
	float4 *hx, *dx; //x,y,z,mass
	float4 *hv, *dv; //vx,vy,vz,charge 
	float4 *hf, *df; //fx,f,fz,type 
	float4 *hx0, *dx0; //Virtual lattice sites, x0, y0, z0 
	
	float4 *dxprev; //Positions at previous neighbour list update 

	unsigned maxneighb; 
	int *neighblist; // neighb indicies + trailing -1s
	int3 *hcrossings, *dcrossings; // Simulation box crossing

	float *epot;  // Potential energy on particle - on device
	float4 *press; // sumdiag,xy,xz,yz contrib pressures per particle - on device 
	float4 *sumpress; // sum of diag, xy, xz, yz
	
	char *ht;   // Type on host only
	
	float *ddist;  // distance travelled by atom - on device
	float dsumdist; // total distance travelled by all atoms - on device
	
	// A few additional members in order to reduce functions API argument list
	// Some of them also in the sys structure
	unsigned nthreads, nblocks;
	float3 lbox;
	unsigned npart, npart_padding;
	
	int *hmolindex, *dmolindex;	// Molecule index 

	// Pair exlusion rules 0 - no exclusion rule, 1 - exclude bonds, 2 - exclude mol.
	unsigned hexclusion_rule, dexclusion_rule; 

	// Data exchange
	struct _sepcusys *sptr;

} sepcupart;


typedef struct _sepcusys {
	// GPU infrastructure 
	unsigned nthreads, nblocks;

	// System size
	unsigned npart, npart_padding;
	float3 lbox; 

	// Simulation details
	float dt;
	float skin;
	long iteration; 
	bool neighbupdate;

	float *dalpha; // Nose-Hoover thermostat coupling Device only 

	// Thermodynamic state of the system
	float3 *henergies, *denergies;  // ekin, epot, momentum
	float ekin, epot, etot;
	float temp;
	
	bool molprop;
	unsigned molpropinterval;

	bool cmflag;

	// Data exchange
	struct _sepcupart *pptr;
	struct _sepcumol *mptr;

} sepcusys;



#endif

