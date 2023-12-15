
#include "sepcudaprfrc.h"

__inline__  __device__ float sep_cuda_wrap(float x, float lbox){
	
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

// I have no clue why I need to define these kernels in this particular file scope
// It will not work if I just use the kernels declared and defined in sepcudamisc
// In sepcudaintgr calling kernles from sepcudamisc works
// THIS MUST BE SOLVED
#ifdef OCTAVE
__global__ void oct_sep_cuda_calc_dist(float *dist, float4 *p, float4 *pprev, float3 lbox, unsigned npart){

	int pidx = blockDim.x * blockIdx.x + threadIdx.x;
	
	if ( pidx < npart ){
		float dx = pprev[pidx].x - p[pidx].x; dx = sep_cuda_wrap(dx, lbox.x);
		float dy = pprev[pidx].y - p[pidx].y; dy = sep_cuda_wrap(dy, lbox.y);
		float dz = pprev[pidx].z - p[pidx].z; dz = sep_cuda_wrap(dz, lbox.z);

		dist[pidx] = sqrtf(dx*dx + dy*dy + dz*dz);
	}

}

__global__ void oct_sep_cuda_sumdistance(float *totalsum, float *dist, unsigned npart){
	
	__shared__ float sum;
	if (threadIdx.x==0) sum=.0f;
	__syncthreads();
	
	int id = blockIdx.x*blockDim.x + threadIdx.x;
	if ( id < npart ) atomicAdd(&sum, dist[id]);
	__syncthreads();
	
	if ( threadIdx.x == 0 ) {
		atomicAdd(totalsum, sum);
	}
}

__global__ void oct_sep_cuda_set_prevpos(float4 *p, float4 *pprev, unsigned npart){
        int pidx = blockDim.x * blockIdx.x + threadIdx.x;
                
        if ( pidx < npart ){
                pprev[pidx].x = p[pidx].x; pprev[pidx].y = p[pidx].y; pprev[pidx].z = p[pidx].z;
        }

}


__global__ void oct_sep_cuda_setvalue(float *variable, float value){*variable = value;}


#endif



// Host functions

void sep_cuda_check_neighblist(sepcupart *ptr, float skin){

	const int nb = ptr->nblocks;
	const int nt = ptr->nthreads;

#ifdef OCTAVE
	oct_sep_cuda_calc_dist<<<nb, nt>>>(ptr->ddist, ptr->dx, ptr->dxprev, ptr->lbox, ptr->npart);
	oct_sep_cuda_sumdistance<<<nb, nt>>>(&(ptr->dsumdist), ptr->ddist, ptr->npart);
#else
	sep_cuda_calc_dist<<<nb, nt>>>(ptr->ddist, ptr->dx, ptr->dxprev, ptr->lbox, ptr->npart);
	sep_cuda_sum_distance<<<nb, nt>>>(&(ptr->dsumdist), ptr->ddist, ptr->npart);
#endif

	//cudaDeviceSynchronize();
	
	float sumdr=0.0f;
	cudaMemcpy(&sumdr, &(ptr->dsumdist), sizeof(float), cudaMemcpyDeviceToHost);
	
	float avsumdr = sumdr/ptr->npart;
		
	if ( avsumdr > skin ){
#ifdef OCTAVE
		oct_sep_cuda_setvalue<<<1,1>>>(&(ptr->dsumdist), 0);
		oct_sep_cuda_set_prevpos<<<nb, nt>>>(ptr->dx, ptr->dxprev, ptr->npart);
#else 
		sep_cuda_set_value<<<1,1>>>(&(ptr->dsumdist), 0.0);
		sep_cuda_set_prevpos<<<nb, nt>>>(ptr->dx, ptr->dxprev, ptr->npart);
#endif

		//cudaDeviceSynchronize();
		ptr->sptr->neighbupdate = true;
	}	
	else 
		ptr->sptr->neighbupdate = false;


}


void sep_cuda_set_exclusion(sepcupart *aptr, const char rule[]){
	
	if ( strcmp(rule, "bonds")==0 ){
		aptr->hexclusion_rule = SEP_CUDA_EXCL_BONDS;
	}
	else if (strcmp(rule, "molecule")==0 ){
		aptr->hexclusion_rule = SEP_CUDA_EXCL_MOLECULE;
	}
	else {
		fprintf(stderr, "Not valid exclusion rule\n");
	}
	
	size_t nbytes = sizeof(unsigned);
	cudaMemcpy(&(aptr->dexclusion_rule), &(aptr->hexclusion_rule), nbytes, cudaMemcpyHostToDevice);
	
}

// Kernel functions

/* Neighbourlist for particles - no exclusion */
__global__ void sep_cuda_build_neighblist(int *neighlist, float4 *p, float *dist, float cf, 
										  float3 lbox, unsigned nneighmax, unsigned npart) {

	int pidx = blockDim.x * blockIdx.x + threadIdx.x;
		
	if ( pidx < npart ){
		float cfsqr = cf*cf; 
		int arrayOffset = pidx*nneighmax;
	
		float mpx = __ldg(&p[pidx].x); float mpy = __ldg(&p[pidx].y); float mpz = __ldg(&p[pidx].z);

		#pragma unroll	
		for ( int n=0; n<nneighmax; n++ ) neighlist[arrayOffset + n] = -1; //<- this should be optimized 
		
		dist[pidx] = 0.0f;
		
		int shift = 0;
		for ( int tile = 0; tile < gridDim.x; tile++ ) {
		
			for ( int j = 0; j < SEP_CUDA_NTHREADS; j++ ) {
				int idxj = tile*blockDim.x + j;
				
				if ( idxj >= npart )  break;
		
				float dx = mpx - p[idxj].x;  dx = sep_cuda_wrap(dx, lbox.x);
				float dy = mpy - p[idxj].y;  dy = sep_cuda_wrap(dy, lbox.y);
				float dz = mpz - p[idxj].z;  dz = sep_cuda_wrap(dz, lbox.z);
	
				float distSqr = dx*dx + dy*dy + dz*dz;

				if ( distSqr < 2.0*FLT_EPSILON ) continue; // Self contribution
				
				if ( distSqr < cfsqr ) {
						
					if ( shift < nneighmax )
						neighlist[arrayOffset + shift] = idxj;
					else if ( shift >= nneighmax ) {
						printf("Neighbour list generation failed\n");
						return;
					}	
					
					shift++;
				}
			}

			__syncthreads();
			
		}
	}
}
	
/* Neighbourlist for particles excluding particles in same molecule */
__global__ void sep_cuda_build_neighblist(int *neighlist, float *dist, float4 *p, int *molindex, 
										  float cf, float3 lbox, unsigned nneighmax, unsigned npart) {

	int pidx = blockDim.x * blockIdx.x + threadIdx.x;
		
	if ( pidx < npart ){
		float cfsqr = cf*cf; 
		int arrayOffset = pidx*nneighmax;
		int moli = molindex[pidx];
		float mpx =p[pidx].x; float mpy = p[pidx].y; float mpz = p[pidx].z;
	
		#pragma unroll	
		for ( int n=0; n<nneighmax; n++ ) neighlist[arrayOffset + n] = -1; //<- this should be optimized 
		
		// Reset the distance travelled since last update - in danger zone
		dist[pidx] = 0.0f;
		
		int shift = 0;
		for ( int tile = 0; tile < gridDim.x; tile++ ) {

			for ( int j = 0; j < SEP_CUDA_NTHREADS; j++ ) {
				int idxj = tile*blockDim.x + j;
				
				if ( idxj >= npart )  break;
				
				if ( moli == molindex[idxj] ) continue;
	
				float dx = mpx - p[idxj].x; dx = sep_cuda_wrap(dx, lbox.x);
				float dy = mpy - p[idxj].y;  dy = sep_cuda_wrap(dy, lbox.y);
				float dz = mpz - p[idxj].z;  dz = sep_cuda_wrap(dz, lbox.z);
			
				float distSqr = dx*dx + dy*dy + dz*dz;

				if ( distSqr < 2.0*FLT_EPSILON ) continue; // Self contribution
				
				if ( distSqr < cfsqr ) {
						
					if ( shift < nneighmax )
						neighlist[arrayOffset + shift] = idxj;
					else if ( shift >= nneighmax ) {
						printf("Neighbour list generation failed\n");
						return;
					}	
					
					shift++;
				}
			}
		}
	}

}

/* Pair interactions - types specified */
__global__ void sep_cuda_lj(const char type1, const char type2, float4 params, int *neighblist, float4 *pos, float4 *force,
							float *epot, float4 *press, unsigned maxneighb, float3 lbox, const unsigned npart){


	int pidx = blockDim.x * blockIdx.x + threadIdx.x;

	if ( pidx < npart ) {
		
		int itype = __float2int_rd(force[pidx].w);
		int atype = (int)type1; int btype = (int)type2; //cast is stupid!
		
		if ( itype != atype && itype != btype ) return;
		
		float sigma = params.x; 
		float epsilon = params.y; 
		float cf = params.z; //__ldg does not work..?
		float aw = params.w; 

		float cfsqr = cf*cf;
		float awh = 0.5*aw;
		float Epot_shift = 4.0*epsilon*(powf(sigma/cf, 12.) - aw*powf(sigma/cf,6.));
		
		int offset = pidx*maxneighb;
			
		float mpx = __ldg(&pos[pidx].x); float mpy = __ldg(&pos[pidx].y); float mpz = __ldg(&pos[pidx].z);
				
		float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f; 
		float Epot = 0.0f; 
		float4 mpress; mpress.x = mpress.y = mpress.z = mpress.w = 0.0f;

		int n = 0;
		while ( neighblist[n+offset] != -1 ){
			int pjdx = neighblist[n+offset];
			int jtype = __float2int_rd(force[pjdx].w);
			
			if ( (itype == atype && jtype == btype) || (itype == btype && jtype == atype) ){
			
				float dx = mpx - pos[pjdx].x; dx = sep_cuda_wrap(dx, lbox.x);
				float dy = mpy - pos[pjdx].y; dy = sep_cuda_wrap(dy, lbox.y);
				float dz = mpz - pos[pjdx].z; dz = sep_cuda_wrap(dz, lbox.z);
		
				float distSqr = dx*dx + dy*dy + dz*dz;

				if ( distSqr < cfsqr ) {
					float rri = sigma*sigma/distSqr; 
					float rri3 = rri*rri*rri;
					float ft = 48.0*epsilon*rri3*(rri3 - awh)*rri;
				
					Fx += ft*dx; Fy += ft*dy; Fz += ft*dz;
					Epot += 0.5*(4.0*epsilon*rri3*(rri3 - aw) - Epot_shift);
					
					mpress.x += dx*ft*dx + dy*ft*dy + dz*ft*dz; 
					mpress.y += dx*ft*dy; mpress.z += dx*ft*dz; mpress.w += dy*ft*dz;
				}
			}
			
			n++;
		}
		
		force[pidx].x += Fx; force[pidx].y += Fy; force[pidx].z += Fz; 
		epot[pidx] += Epot; 
		
		press[pidx].x += mpress.x;
		press[pidx].y += mpress.y; press[pidx].z += mpress.z; press[pidx].w += mpress.w; 
	}
		
}


/* Pair interactions - all types have same interactions (faster) */
__global__ void sep_cuda_lj(int *neighblist, float4 *pos, float4 *force,
				float *epot, float4 *press, unsigned maxneighb, float3 lbox, const unsigned npart){


	int pidx = blockDim.x * blockIdx.x + threadIdx.x;

	if ( pidx < npart ) {

		float cf = 2.5; 
		float cfsqr = cf*cf;
		float Epot_shift = 4.0*(powf(1.0/cf, 12.) - powf(1.0/cf,6.));

		int offset = pidx*maxneighb;
			
		float mpx = pos[pidx].x; float mpy = pos[pidx].y; float mpz = pos[pidx].z;
	
		float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f; 
		float Epot = 0.0f; 
		float4 mpress; mpress.x = mpress.y = mpress.z = mpress.w = 0.0f;
		
		register int n = 0;
		while ( neighblist[n+offset] != -1 ){
			int pjdx = neighblist[n+offset];
				
			float dx = mpx - pos[pjdx].x; dx = sep_cuda_wrap(dx, lbox.x);
			float dy = mpy - pos[pjdx].y; dy = sep_cuda_wrap(dy, lbox.y);
			float dz = mpz - pos[pjdx].z; dz = sep_cuda_wrap(dz, lbox.z);
	
			float distSqr = dx*dx + dy*dy + dz*dz;

			if ( distSqr < cfsqr ) {
				float rri = 1.0/distSqr; 
				float rri3 = rri*rri*rri;
				float ft =  48.0*rri3*(rri3 - 0.5)*rri; 
				
				Fx += ft*dx; Fy += ft*dy; Fz += ft*dz;
			
			   	Epot += 0.5*(4.0*rri3*(rri3 - 1.0) - Epot_shift);
				mpress.x += dx*ft*dx + dy*ft*dy + dz*ft*dz; 
				mpress.y += dx*ft*dy; mpress.z += dx*ft*dz; mpress.w += dy*ft*dz;
			}
			
			n++;
		}
			
		force[pidx].x += Fx; force[pidx].y += Fy; force[pidx].z += Fz;
		
		epot[pidx] += Epot; 
		press[pidx].x += mpress.x;
		press[pidx].y += mpress.y; press[pidx].z += mpress.z; press[pidx].w += mpress.w; 
	
	}
}



__global__ void sep_cuda_lj_sf(const char type1, const char type2, float3 params, int *neighblist, float4 *pos, float4 *force,
							float *epot, float4 *press, unsigned maxneighb, float3 lbox, const unsigned npart){


	int pidx = blockDim.x * blockIdx.x + threadIdx.x;

	if ( pidx < npart ) {
		
		int itype = __float2int_rd(force[pidx].w);
		int atype = (int)type1; int btype = (int)type2; //cast stupid
		
		if ( itype != atype && itype != btype ) return;
		
		float sigma = params.x; 
		float epsilon = params.y; 
		float cf = params.z; //__ldg does not work..?
		float cfsqr = cf*cf; 
		float Epot_shift = 4.0*epsilon*(powf(sigma/cf, 12.) - powf(sigma/cf,6.));
		float force_shift = 48.0*epsilon*powf(sigma/cf,6.0)*(powf(sigma/cf,3.0) - 0.5)*pow(sigma/cf, 2.0);
		
		int offset = pidx*maxneighb;
			
		float mpx = __ldg(&pos[pidx].x); float mpy = __ldg(&pos[pidx].y); float mpz = __ldg(&pos[pidx].z);
				
		float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f; 
		float Epot = 0.0f; 
		float4 mpress; mpress.x = mpress.y = mpress.z = mpress.w = 0.0f;
		int n = 0;
		while ( neighblist[n+offset] != -1 ){
			int pjdx = neighblist[n+offset];
			int jtype = __float2int_rd(force[pjdx].w);
			
			if ( (itype == atype && jtype == btype) || (itype == btype && jtype == atype) ){
				
				float dx = mpx - pos[pjdx].x; dx = sep_cuda_wrap(dx, lbox.x);
				float dy = mpy - pos[pjdx].y; dy = sep_cuda_wrap(dy, lbox.y);
				float dz = mpz - pos[pjdx].z; dz = sep_cuda_wrap(dz, lbox.z);

				float distSqr = dx*dx + dy*dy + dz*dz;

				if ( distSqr < cfsqr ) {
					float rri = sigma*sigma/distSqr; 
					float rri3 = rri*rri*rri;
					float ft = 48.0*epsilon*rri3*(rri3 - 0.5)*rri + force_shift;
				
					Fx += ft*dx; Fy += ft*dy; Fz += ft*dz;
					Epot += 0.5*(4.0*epsilon*rri3*(rri3 - 1.0) - Epot_shift);
					mpress.x += dx*ft*dx + dy*ft*dy + dz*ft*dz; 
					mpress.y += dx*ft*dy; mpress.z += dx*ft*dz; mpress.w += dy*ft*dz;
				}
			}
			
			n++;
		}
		
		force[pidx].x += Fx; force[pidx].y += Fy; force[pidx].z += Fz;
		epot[pidx] += Epot; 
		press[pidx].x += mpress.x;
		press[pidx].y += mpress.y; press[pidx].z += mpress.z; press[pidx].w += mpress.w; 
	}
	
}



__global__ void sep_cuda_sf(float cf, int *neighblist, float4 *pos, float4 *vel, float4 *force,
						float *epot, float4 *press, unsigned maxneighb, float3 lbox, const unsigned npart){

	__const__ int pidx = blockDim.x * blockIdx.x + threadIdx.x;

	if ( pidx < npart ) {
		
		float cfsqr = cf*cf;
		float icf = 1.0/cf;
		float icf2 = 1.0/cfsqr;
		
		int offset = pidx*maxneighb;
			
		float mpx = __ldg(&pos[pidx].x); 
		float mpy = __ldg(&pos[pidx].y); 
		float mpz = __ldg(&pos[pidx].z);
				
		float Fx = 0.0; float Fy = 0.0; float Fz = 0.0; float Epot = 0.0;		
		float4 mpress; mpress.x = mpress.y = mpress.z = mpress.w = 0.0f;
		
		int n = 0;
		while ( neighblist[n+offset] != -1 ){
			int pjdx = neighblist[n+offset];
				
			float dx = mpx - pos[pjdx].x; dx = sep_cuda_wrap(dx, lbox.x);
			float dy = mpy - pos[pjdx].y; dy = sep_cuda_wrap(dy, lbox.y);
			float dz = mpz - pos[pjdx].z; dz = sep_cuda_wrap(dz, lbox.z);

			float distSqr = dx*dx + dy*dy + dz*dz;

			if ( distSqr < cfsqr ) {
				float zizj = vel[pidx].w*vel[pjdx].w;
				float dist = sqrtf(distSqr); 
				float ft = zizj*(1.0/distSqr - icf2)/dist; 
				
				Fx += ft*dx; Fy += ft*dy; Fz += ft*dz;
				
				Epot += 0.5*zizj*(1.0/dist + (dist-cf)*icf2 - icf);
				mpress.x += dx*ft*dx + dy*ft*dy + dz*ft*dz; 
				mpress.y += dx*ft*dy; mpress.z += dx*ft*dz; mpress.w += dy*ft*dz;
			}

			n++;
		}
		
		force[pidx].x += Fx; force[pidx].y += Fy; force[pidx].z += Fz;
		epot[pidx] += Epot;	
		press[pidx].x += mpress.x;
		press[pidx].y += mpress.y; press[pidx].z += mpress.z; press[pidx].w += mpress.w;
	}	
		
}

// Molecular force calculations (for pressure etc)
__global__ void sep_cuda_calc_molforce(float3 *mforce,  const char type1, const char type2, float4 params, float4 *pos, 
										int *neighblist,  unsigned maxneighb, float4 *force, 
										float3 lbox, int *molindex, unsigned nmol, const unsigned npart) {


	int pidx = blockDim.x * blockIdx.x + threadIdx.x;

	if ( pidx < npart ) {
		
		int itype = __float2int_rd(force[pidx].w);
		int atype = (int)type1; int btype = (int)type2; //cast is stupid!
		
		if ( itype != atype && itype != btype ) return;
		
		float sigma = params.x;	float epsilon = params.y; float cf = params.z; float aw = params.w;
		float cfsqr = cf*cf; float awh = 0.5*aw;
		int molidx = molindex[pidx];

		int offset = pidx*maxneighb;
			
		float mpx = __ldg(&pos[pidx].x); float mpy = __ldg(&pos[pidx].y); float mpz = __ldg(&pos[pidx].z);
				
		int n = 0;
		while ( neighblist[n+offset] != -1 ){
			int pjdx = neighblist[n+offset];
			int jmolidx = molindex[pjdx]; 

			if ( jmolidx == molidx ) { n++; continue; }

			int jtype = __float2int_rd(force[pjdx].w);
				
			if ( (itype == atype && jtype == btype) || (itype == btype && jtype == atype) ){
				
				float dx = mpx - pos[pjdx].x; dx = sep_cuda_wrap(dx, lbox.x);
				float dy = mpy - pos[pjdx].y; dy = sep_cuda_wrap(dy, lbox.y);
				float dz = mpz - pos[pjdx].z; dz = sep_cuda_wrap(dz, lbox.z);

				float distSqr = dx*dx + dy*dy + dz*dz;

				if ( distSqr < cfsqr ) {
					float rri = sigma*sigma/distSqr; 
					float rri3 = rri*rri*rri;
					float ft = 48.0*epsilon*rri3*(rri3 - awh)*rri;
			
					unsigned offset = molidx*nmol + jmolidx;
					float Fx = ft*dx; float Fy = ft*dy; float Fz = ft*dz;
					
					atomicAdd(&mforce[offset].x, Fx);
					atomicAdd(&mforce[offset].y, Fy);
					atomicAdd(&mforce[offset].z, Fz);
				}
			}
			
			n++;
		}
	}
		
}

__global__ void sep_cuda_calc_molforce(float3 *mforce, float cf, int *neighblist, float4 *pos, float4 *vel, 
										 unsigned maxneighb, int *molindex, unsigned nmols, 
										 float3 lbox, const unsigned npart){

	__const__ int pidx = blockDim.x * blockIdx.x + threadIdx.x;

	if ( pidx < npart ) {
		
		float cfsqr = cf*cf;
		float icf2 = 1.0/cfsqr;
		
		int offset = pidx*maxneighb;
			
		int molidx = molindex[pidx];

		float mpx = __ldg(&pos[pidx].x); 
		float mpy = __ldg(&pos[pidx].y); 
		float mpz = __ldg(&pos[pidx].z);
		
		int n = 0;
		while ( neighblist[n+offset] != -1 ){
			int pjdx = neighblist[n+offset];
			int jmolidx = molindex[pjdx]; 
	
			if ( jmolidx == molidx ) { n++; continue; }

			float dx = mpx - pos[pjdx].x; dx = sep_cuda_wrap(dx, lbox.x);
			float dy = mpy - pos[pjdx].y; dy = sep_cuda_wrap(dy, lbox.y);
			float dz = mpz - pos[pjdx].z; dz = sep_cuda_wrap(dz, lbox.z);

			float distSqr = dx*dx + dy*dy + dz*dz;

			if ( distSqr < cfsqr ) {
				float zizj = vel[pidx].w*vel[pjdx].w;
				float dist = sqrtf(distSqr); 
				float ft = zizj*(1.0/distSqr - icf2)/dist; 
				
				unsigned offset = nmols*molidx + jmolidx;		
				float Fx = ft*dx; float Fy = ft*dy; float Fz = ft*dz;

				atomicAdd(&mforce[offset].x, Fx);
				atomicAdd(&mforce[offset].y, Fy);
				atomicAdd(&mforce[offset].z, Fz);
			}

			n++;
		}
	}	
		
}

__global__ void sep_cuda_lattice_force(const char type, float springConstant, float4 *pos, float4 *pos0, float4 *force,
									   float3 lbox, const unsigned npart){

	
	unsigned pidx = blockDim.x * blockIdx.x + threadIdx.x;
	int itype = __float2int_rd(force[pidx].w);
		
	if ( pidx < npart && itype == (int)type ){
		
		float dx = pos[pidx].x - pos0[pidx].x; dx = sep_cuda_wrap(dx, lbox.x);
		float dy = pos[pidx].y - pos0[pidx].y; dy = sep_cuda_wrap(dy, lbox.y);
		float dz = pos[pidx].z - pos0[pidx].z; dz = sep_cuda_wrap(dz, lbox.z);

		force[pidx].x = - springConstant*dx;
		force[pidx].y = - springConstant*dy;
		force[pidx].z = - springConstant*dz;
		
	}
	
}


// Wrapper section

void sep_cuda_force_lj(sepcupart *pptr, const char types[], float params[4]){
	const int nb = pptr->nblocks; 
	const int nt = pptr->nthreads;

	float4 ljparams = make_float4(params[0],params[1],params[2],params[3]);

	sep_cuda_lj<<<nb, nt>>>
		(types[0], types[1], ljparams, pptr->neighblist, pptr->dx, pptr->df,  
				pptr->epot, pptr->press, pptr->maxneighb, pptr->lbox, pptr->npart);

	// Calculates the molecular forces for molecular stress tensor
	if ( pptr->sptr->molprop && pptr->sptr->iteration%pptr->sptr->molpropinterval==0 ){
		sep_cuda_calc_molforce<<<nb,nt>>>(pptr->sptr->mptr->dfij, types[0], types[1], ljparams, 
					pptr->dx, pptr->neighblist, pptr->maxneighb, pptr->df, 
					pptr->sptr->lbox, pptr->dmolindex, pptr->sptr->mptr->nmols, pptr->sptr->npart);
	}

	//cudaDeviceSynchronize();

}

void sep_cuda_force_lj(sepcupart *pptr){
	const int nb = pptr->nblocks; 
	const int nt = pptr->nthreads;
	
	sep_cuda_lj<<<nb, nt>>>
		(pptr->neighblist, pptr->dx, pptr->df, pptr->epot, pptr->press, pptr->maxneighb, pptr->lbox, pptr->npart);
		
	//cudaDeviceSynchronize();

}


void sep_cuda_force_lj_sf(sepcupart *pptr, const char types[], float params[3]){
	const int nb = pptr->nblocks; 
	const int nt = pptr->nthreads;
	
	float3 ljparams = make_float3(params[0],params[1],params[2]);
	
	sep_cuda_lj_sf<<<nb, nt>>>
		(types[0], types[1], ljparams, pptr->neighblist, pptr->dx, pptr->df, pptr->epot, 
									pptr->press, pptr->maxneighb, pptr->lbox, pptr->npart);
	//cudaDeviceSynchronize();

}



void sep_cuda_force_sf(sepcupart *pptr, const float cf){
	const int nb = pptr->nblocks; 
	const int nt = pptr->nthreads;
	
	sep_cuda_sf<<<nb,nt>>>
		(cf, pptr->neighblist, pptr->dx, pptr->dv, pptr->df, pptr->epot, 
											pptr->press, pptr->maxneighb, pptr->lbox, pptr->npart);

	if ( pptr->sptr->molprop && pptr->sptr->iteration%pptr->sptr->molpropinterval==0 ){
		sep_cuda_calc_molforce<<<nb, nt>>>(pptr->sptr->mptr->dfij, cf, pptr->neighblist, pptr->dx, pptr->dv, 
										 pptr->maxneighb, pptr->dmolindex, pptr->sptr->mptr->nmols, 
										 pptr->lbox, pptr->npart);
	}

	//cudaDeviceSynchronize();

}

void sep_cuda_update_neighblist(sepcupart *pptr, float maxcutoff){

	if ( pptr->sptr->neighbupdate ){
	
		const int nb = pptr->sptr->nblocks; 
		const int nt = pptr->sptr->nthreads;

		if ( pptr->hexclusion_rule == SEP_CUDA_EXCL_NONE ) {
			sep_cuda_build_neighblist<<<nb, nt>>>
				(pptr->neighblist, pptr->dx, pptr->ddist, pptr->sptr->skin+maxcutoff, pptr->lbox, pptr->maxneighb, pptr->npart);
		}
		else if ( pptr->hexclusion_rule == SEP_CUDA_EXCL_MOLECULE ) {
			sep_cuda_build_neighblist<<<nb, nt>>>
				(pptr->neighblist, pptr->ddist, pptr->dx, pptr->dmolindex, pptr->sptr->skin+maxcutoff, pptr->lbox, pptr->maxneighb,pptr->npart);
		}
		else {
			fprintf(stderr, "Exclusion rule invalid");
		}

		pptr->sptr->neighbupdate = false;

		//cudaDeviceSynchronize();
	}
}

void sep_cuda_force_lattice(sepcupart *pptr, const char type, float springConstant){
	const int nb = pptr->nblocks; 
	const int nt = pptr->nthreads;
	
	sep_cuda_lattice_force<<<nb, nt>>>
		(type, springConstant, pptr->dx, pptr->dx0, pptr->df, pptr->lbox, pptr->npart);
		
	//cudaDeviceSynchronize();

}


