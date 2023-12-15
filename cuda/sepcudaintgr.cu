
#include "sepcudaintgr.h"


// Again, I have no clue why these function must be defined here
// Octave's file scope versus CUDA...

#ifdef OCTAVE

__inline__ __device__ float sep_cuda_dot(float4 a){
	
	return (a.x*a.x + a.y*a.y + a.z*a.z);
	
}


__global__ void oct_sep_cuda_sumenergies(float3 *totalsum, float4* dx, float4 *dv, float4 *df, 
										 float dt, float *epot, unsigned npart){

	int id = blockIdx.x*blockDim.x + threadIdx.x;
	__shared__ float3 sums;

	if ( threadIdx.x==0 ) {
		sums.x = sums.y = sums.z = 0.0f;
	}
	__syncthreads();

	if ( id < npart ){
		float4 vel; 
		vel.x =  dv[id].x - 0.5*dt*df[id].x/dx[id].w;
		vel.y =  dv[id].y - 0.5*dt*df[id].y/dx[id].w;
		vel.z =  dv[id].z - 0.5*dt*df[id].z/dx[id].w;
		
		float mykin = 0.5*sep_cuda_dot(vel)*dx[id].w;
		float mymom = (dv[id].x + dv[id].y + dv[id].z)*dx[id].w;
		
		atomicAdd(&sums.x, mykin);
		atomicAdd(&sums.y, epot[id]);
		atomicAdd(&sums.z, mymom);
	}

	__syncthreads();
	
	if ( id < npart && threadIdx.x == 0 ) {
		atomicAdd(&(totalsum->x), sums.x);
		atomicAdd(&(totalsum->y), sums.y);
		atomicAdd(&(totalsum->z), sums.z);
	}
	
}

__global__ void oct_sep_cuda_sum_ekin(float3 *totalsum, const char type, float4* dx, float4 *dv, float4 *df, 
									  float dt, float *epot, unsigned npart){

	int id = blockIdx.x*blockDim.x + threadIdx.x;
	__shared__ float sumekin;
	__shared__ int numtype;

	if ( threadIdx.x==0 ){
	   	sumekin = 0.0f;
		numtype = 0;
	}
	__syncthreads();

	int itype = __float2int_rd(df[id].w);

	if ( id < npart && itype==(int)type ){
		float4 vel; 
		vel.x =  dv[id].x - 0.5*dt*df[id].x/dx[id].w;
		vel.y =  dv[id].y - 0.5*dt*df[id].y/dx[id].w;
		vel.z =  dv[id].z - 0.5*dt*df[id].z/dx[id].w;
		
		float myekin = 0.5*sep_cuda_dot(vel)*dx[id].w;
		
		atomicAdd(&sumekin, myekin);
		atomicAdd(&numtype, 1.0);
	}

	__syncthreads();
	
	if ( id < npart && threadIdx.x == 0 ){ 
		atomicAdd(&(totalsum->x), sumekin);
		atomicAdd(&(totalsum->y), numtype);
	}
}


__global__ void oct_sep_cuda_reset_variable(float3 *a){

	a->x = a->y = a->z = 0.0f;

}



#endif

__inline__ __device__ float sep_cuda_wrap(float x, float lbox){
	
	if ( x > 0.5*lbox ) 
		x -= lbox;
	else if  ( x < -0.5*lbox ) 
		x += lbox;
	
	return x;
}

__inline__ __device__ float sep_cuda_periodic(float x, float lbox, int *crossing){
	
	if ( x > lbox ) {
		x -= lbox;  
		*crossing = *crossing + 1;
	}
	else if  ( x < 0 ) {
		x += lbox;
		*crossing = *crossing - 1;
	}
	
	return x;
}


__global__ void sep_cuda_leapfrog(float4 *pos, float4 *vel, 
		  float4 *force, float *dist, int3 *crossing, float dt, float3 lbox, unsigned npart){

	int i = blockDim.x * blockIdx.x + threadIdx.x;

	// Local variables not even necessary for speed up	
	float4 mypos = make_float4(pos[i].x, pos[i].y, pos[i].z, pos[i].w);
	float4 myvel = make_float4(vel[i].x, vel[i].y, vel[i].z, vel[i].w);

	if ( i < npart ) {
		float imass = 1.0/mypos.w;
		
		myvel.x += force[i].x*imass*dt; 
		myvel.y += force[i].y*imass*dt;
		myvel.z += force[i].z*imass*dt;
		
		mypos.x += myvel.x*dt;
		mypos.x = sep_cuda_periodic(mypos.x, lbox.x, &(crossing[i].x));
		
		mypos.y += myvel.y*dt;
		mypos.y = sep_cuda_periodic(mypos.y, lbox.y, &(crossing[i].y));
	
		mypos.z += myvel.z*dt;
		mypos.z = sep_cuda_periodic(mypos.z, lbox.z, &(crossing[i].z));
	
		pos[i].x = mypos.x; pos[i].y = mypos.y; pos[i].z = mypos.z;
		vel[i].x = myvel.x; vel[i].y = myvel.y; vel[i].z = myvel.z; 
	}
	
}


__global__ void sep_cuda_update_nosehoover(float *alpha, float3 *denergies, float temp0, 
										   float tau, float dt, unsigned int npart){

	float temp = (2.0/3.0)*denergies->x/npart; 

	*alpha = *alpha + dt/(tau*tau)*(temp/temp0 - 1.0);

}


__global__ void sep_cuda_update_nosehoover(float *alpha, float3 *denergies, float temp0, 
										   float tau, float dt){

	float temp = (2.0/3.0)*denergies->x/denergies->y; 

	*alpha = *alpha + dt/(tau*tau)*(temp/temp0 - 1.0);

}


__global__ void sep_cuda_nosehoover(float *alpha, float4 *pos, float4 *vel, float4 *force, unsigned npart){
	
	unsigned id = blockIdx.x*blockDim.x + threadIdx.x;

	if ( id < npart ){	
		float fac = (*alpha)*pos[id].w;
		force[id].x -= fac*vel[id].x; 
		force[id].y -= fac*vel[id].y; 
		force[id].z -= fac*vel[id].z;		
	}
}

__global__ void sep_cuda_nosehoover(float *alpha, const char type, float4 *pos, float4 *vel, float4 *force, unsigned npart){
	
	unsigned id = blockIdx.x*blockDim.x + threadIdx.x;
	int itype = __float2int_rd(force[id].w);

	if ( id < npart && itype == (int)type ){	
		float fac = (*alpha)*pos[id].w;
		force[id].x -= fac*vel[id].x; 
		force[id].y -= fac*vel[id].y; 
		force[id].z -= fac*vel[id].z;		
	}
}


void sep_cuda_thermostat_nh(sepcupart *pptr, float temp0, float tau){
	const int nb = pptr->sptr->nblocks; 
	const int nt = pptr->sptr->nthreads;
	
	// Get current system kinetic energy
#ifdef OCTAVE
	oct_sep_cuda_sumenergies<<<nb,nt>>>
		(pptr->sptr->denergies, pptr->dx, pptr->dv, pptr->df, pptr->sptr->dt, pptr->epot, pptr->sptr->npart);
#else 
	sep_cuda_sum_energies<<<nb,nt>>>
		(pptr->sptr->denergies, pptr->dx, pptr->dv, pptr->df, pptr->sptr->dt, pptr->epot, pptr->sptr->npart);
#endif

	//cudaDeviceSynchronize();
	
	// Update nh-alpha dynamics (single thread)
	sep_cuda_update_nosehoover<<<1,1>>>
		(pptr->sptr->dalpha, pptr->sptr->denergies, temp0, tau, pptr->sptr->dt, pptr->sptr->npart);
	//cudaDeviceSynchronize();

	// Add thermostat force
	sep_cuda_nosehoover<<<nb, nt>>>
		(pptr->sptr->dalpha, pptr->dx, pptr->dv, pptr->df, pptr->sptr->npart);
	//cudaDeviceSynchronize();		

	sep_cuda_reset_variable<<<1,1>>>(pptr->sptr->denergies);	
}

 
void sep_cuda_thermostat_nh(const char type, sepcupart *pptr, float temp0, float tau){
	const int nb = pptr->sptr->nblocks; 
	const int nt = pptr->sptr->nthreads;
	
	// Get current system kinetic energy
#ifdef OCTAVE
	oct_sep_cuda_sum_ekin<<<nb,nt>>>
		(pptr->sptr->denergies, type, pptr->dx, pptr->dv, pptr->df, pptr->sptr->dt, pptr->epot, pptr->sptr->npart);
#else 
	sep_cuda_sum_ekin<<<nb,nt>>>
		(pptr->sptr->denergies, type, pptr->dx, pptr->dv, pptr->df, pptr->sptr->dt, pptr->epot, pptr->sptr->npart);
#endif

	//cudaDeviceSynchronize();
	
	// Update nh-alpha dynamics (single thread)
	sep_cuda_update_nosehoover<<<1,1>>>
		(pptr->sptr->dalpha, pptr->sptr->denergies, temp0, tau, pptr->sptr->dt);
	//cudaDeviceSynchronize();

	// Add thermostat force
	sep_cuda_nosehoover<<<nb, nt>>>
		(pptr->sptr->dalpha, type,  pptr->dx, pptr->dv, pptr->df, pptr->sptr->npart);
	cudaDeviceSynchronize();		

#ifdef OCTAVE
	oct_sep_cuda_reset_variable<<<1,1>>>(pptr->sptr->denergies);
#else
	sep_cuda_reset_variable<<<1,1>>>(pptr->sptr->denergies);	
#endif

}

  
void sep_cuda_integrate_leapfrog(sepcupart *pptr){
	const int nb = pptr->sptr->nblocks; 
	const int nt = pptr->sptr->nthreads;

	sep_cuda_leapfrog<<<nb, nt>>>
		(pptr->dx, pptr->dv, pptr->df, pptr->ddist, pptr->dcrossings, pptr->sptr->dt, pptr->lbox, pptr->npart);
	//cudaDeviceSynchronize();
	
}


