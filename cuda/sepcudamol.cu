#include "sepcudamol.h"

// Device functions

__device__ float sep_cuda_mol_dot(float4 a){
	
	return (a.x*a.x + a.y*a.y + a.z*a.z);
	
}

__device__ float sep_cuda_mol_dot(float3 a, float3 b){
	
	return (a.x*b.x + a.y*b.y + a.z*b.z);
	
}

__device__ float sep_cuda_mol_wrap(float x, float lbox){
	
	if ( x > 0.5*lbox ) 
		x -= lbox;
	else if  ( x < -0.5*lbox ) 
		x += lbox;
	
	return x;
}



sepcumol * sep_cuda_init_mol(void){
	
	sepcumol *mptr = (sepcumol *)malloc(sizeof(sepcumol));
	if ( mptr==NULL ) sep_cuda_mem_error();
	
	mptr->nmols = 0; 
	
	return mptr;
}



FILE *sep_cuda_set_file_pointer(FILE *fptr, const char *section){
	char line[256];

	do {

		if ( fgets(line, 256, fptr) == NULL ) break;

		if ( strcmp(line, section) == 0 ){
			if ( fgets(line, 256, fptr)==NULL )
				sep_cuda_file_error();
			return fptr;
		}

	}  while ( !feof(fptr) );

	return NULL;
}

void sep_cuda_make_atomslist(sepcumol *mptr, sepcupart *pptr){

	unsigned nmols = mptr->nmols;
	unsigned npart = pptr->sptr->npart;

	int *offset = (int *)malloc(sizeof(int)*nmols);
	if ( offset==NULL ) sep_cuda_mem_error();

	// Initializing
	for ( unsigned m=0; m<nmols; m++ ) offset[m]=0;
	for ( unsigned n=0; n<nmols*SEP_CUDA_MAXNUAU; n++ ) mptr->alist[n] = -1;

	for ( unsigned n=0; n<npart; n++ ){
		int molidx = pptr->hmolindex[n];
		if ( molidx > -1 ){
	    	unsigned idx = SEP_CUDA_MAXNUAU*molidx + offset[molidx];
			mptr->alist[idx]=n;
			offset[molidx]++;		
		}
	}

	free(offset);
}

void sep_cuda_read_bonds(sepcupart *pptr, sepcumol *mptr, const char *file, const char opt){
	const char section[] = {'[', ' ', 'b', 'o', 'n', 'd', 's', ' ', ']','\n', '\0'};
	char line[256];
	fpos_t pos_file;
	unsigned moli, a, b, type;

	mptr->nbonds = 0; 
	
	FILE *fptr = fopen(file, "r");
	if ( fptr == NULL ) sep_cuda_file_error();
	 
	// Find the 'bonds' section 
	fptr = sep_cuda_set_file_pointer(fptr, section);
	if ( fptr == NULL ) {
		if ( opt=='v' ){
		   fprintf(stdout, "No bonds found\n");
	   	}
 		return;
	}		
	
	// We *must* init the pointer since it will be free no matter if the read is sucessful or not
	mptr->hblist = (unsigned *)malloc(0);  // stupid-cast
	if ( mptr->hblist == NULL ) sep_cuda_mem_error();

	do {

		fgetpos(fptr, &pos_file); 
		if ( fgets(line, 256, fptr) == NULL ) sep_cuda_file_error();
		
		if ( line[0] == '[' ) {
			break;
		}
		else {
      
			fsetpos(fptr, &pos_file); 
      
			int sc = fscanf(fptr, "%u%u%u%u\n", &moli, &a, &b, &type);
			if ( sc != 4 ) sep_cuda_file_error();
     
			(mptr->nbonds) ++;
    
			mptr->hblist = (unsigned *)realloc((mptr->hblist), sizeof(unsigned)*3*mptr->nbonds);
			if ( mptr->hblist == NULL ) sep_cuda_mem_error();
      
			int index0 = (mptr->nbonds-1)*3;
      
			pptr->hmolindex[a]=moli;
			pptr->hmolindex[b]=moli;
			
			mptr->hblist[index0] = a;
			mptr->hblist[index0+1] = b;
			mptr->hblist[index0+2] = type;
			
			/* Broken
			sep_cuda_set_hexclusion(pptr, a, b); 
			sep_cuda_set_hexclusion(pptr, b, a);
			*/
			
			if ( moli > mptr->nmols ) mptr->nmols = moli;
		}
	} while ( !feof(fptr) ); 
	
	fclose(fptr);

	// Since number of mols is one bigger than the index
	(mptr->nmols)++; 
	mptr->nbondblocks = mptr->nbonds/SEP_CUDA_NTHREADS + 1 ;
	
	if ( opt=='v' ){	
		fprintf(stdout, "Succesfully read 'bond' section in file %s -> ", file);
    	fprintf(stdout, "Found %d molecule(s) and %d bond(s)\n", mptr->nmols, mptr->nbonds);
		fprintf(stdout, "Copying to device\n");
	}

	size_t nbytes =  3*(mptr->nbonds)*sizeof(unsigned int);
	if ( cudaMalloc((void **)&(mptr->dblist),nbytes) == cudaErrorMemoryAllocation )
		sep_cuda_mem_error();
	
	cudaMemcpy(mptr->dblist, mptr->hblist, nbytes, cudaMemcpyHostToDevice);

	nbytes = pptr->npart_padding*sizeof(int);
	cudaMemcpy(pptr->dmolindex, pptr->hmolindex, nbytes, cudaMemcpyHostToDevice);

	nbytes =  SEP_CUDA_MAXNUAU*(mptr->nmols)*sizeof(int);
	if ( cudaMallocHost((void **)&(mptr->alist),nbytes) == cudaErrorMemoryAllocation )
		sep_cuda_mem_error();

	sep_cuda_make_atomslist(mptr, pptr);

	// Allocating for molecular force, positions, and velocities
	nbytes = mptr->nmols*sizeof(float3);
	if ( cudaMalloc((void **)&(mptr->df),nbytes) == cudaErrorMemoryAllocation || 
		 cudaMallocHost((void **)&(mptr->hf),nbytes) == cudaErrorMemoryAllocation ||
		// cudaMalloc((void **)&(mptr->dx),nbytes) == cudaErrorMemoryAllocation || 
		 cudaMallocHost((void **)&(mptr->hx),nbytes) == cudaErrorMemoryAllocation ||
		// cudaMalloc((void **)&(mptr->dv),nbytes) == cudaErrorMemoryAllocation || 
		 cudaMallocHost((void **)&(mptr->hv),nbytes) == cudaErrorMemoryAllocation ||
	 	// cudaMalloc((void **)&(mptr->dpel),nbytes) == cudaErrorMemoryAllocation || 
		 cudaMallocHost((void **)&(mptr->hpel),nbytes) == cudaErrorMemoryAllocation )
		sep_cuda_mem_error();

	if ( cudaMallocHost((void **)&(mptr->masses), mptr->nmols*sizeof(float)) == cudaErrorMemoryAllocation )
		sep_cuda_mem_error();
	
	nbytes = mptr->nmols*mptr->nmols*sizeof(float3);
	if ( cudaMalloc((void **)&(mptr->dfij),nbytes) == cudaErrorMemoryAllocation || 
		 cudaMallocHost((void **)&(mptr->hfij),nbytes) == cudaErrorMemoryAllocation )
		sep_cuda_mem_error();

	// Now sys structure can access molecular info
	pptr->sptr->mptr=mptr;
}


void sep_cuda_free_bonds(sepcumol *mptr){

	if ( mptr->nbonds == 0 ) return;

	cudaFreeHost(mptr->hblist);	cudaFree(mptr->dblist);

	cudaFreeHost(mptr->hf);	cudaFree(mptr->df);
	cudaFreeHost(mptr->hv);	//cudaFree(mptr->dv);
	cudaFreeHost(mptr->hx);	//cudaFree(mptr->dx);
	cudaFreeHost(mptr->hpel);	//cudaFree(mptr->dpel);

	cudaFreeHost(mptr->hfij); cudaFree(mptr->dfij);

	cudaFreeHost(mptr->alist);

	cudaFreeHost(mptr->masses);
}


void sep_cuda_read_angles(sepcumol *mptr, const char *file, const char opt){
	const char section[] = {'[', ' ', 'a', 'n', 'g', 'l', 'e', 's', ' ', ']','\n', '\0'};
	char line[256];
	fpos_t pos_file;
	unsigned moli, a, b, c, type;

	if ( mptr->nmols == 0 ) {
		fprintf(stderr, "Bond section must be read before angle section");
		exit(EXIT_FAILURE);
	}
	
	mptr->nangles = 0; 

	FILE *fptr = fopen(file, "r");
	if ( fptr == NULL ) sep_cuda_file_error();
	 
	// Find the 'angles' section 
	fptr = sep_cuda_set_file_pointer(fptr, section);
	if ( fptr == NULL ) {
		if ( opt=='v' ){
		   fprintf(stdout, "No angles found");
	   	}
 		return;
	}		
	
	// We *must* init the pointer since it will be free no matter if the read is sucessful or not
	mptr->halist = (unsigned *)malloc(0);
	if ( mptr->halist == NULL ) sep_cuda_mem_error();
 
	do {

		fgetpos(fptr, &pos_file); 
		if ( fgets(line, 256, fptr) == NULL ) sep_cuda_file_error();
		
		if ( line[0] == '[' ) {
			break;
		}
		else {
      
			fsetpos(fptr, &pos_file); 
      
			int sc = fscanf(fptr, "%u%u%u%u%u\n", &moli, &a, &b, &c, &type);
			if ( sc != 5 ) sep_cuda_file_error();
     
			(mptr->nangles) ++;
    
			mptr->halist = (unsigned *)realloc((mptr->halist), sizeof(unsigned)*4*mptr->nangles);
			if ( mptr->halist == NULL ) sep_cuda_mem_error();
      
			int index0 = (mptr->nangles-1)*4;
      
			mptr->halist[index0] = a;
			mptr->halist[index0+1] = b;
			mptr->halist[index0+2] = c;
			mptr->halist[index0+3] = type;
			
			/* Broken
			sep_cuda_set_hexclusion(pptr, a, b); sep_cuda_set_hexclusion(pptr, a, c); 
			sep_cuda_set_hexclusion(pptr, b, a); sep_cuda_set_hexclusion(pptr, b, c); 
			sep_cuda_set_hexclusion(pptr, c, a); sep_cuda_set_hexclusion(pptr, c, b);			
			*/
		}
	} while ( !feof(fptr) ); 
	
	fclose(fptr);

	// Since number of mols is one bigger than the index
	mptr->nangleblocks = mptr->nangles/SEP_CUDA_NTHREADS + 1 ;
	
	if ( opt=='v') {
		fprintf(stdout, "Succesfully read 'angle' section in file %s -> ", file);
    	fprintf(stdout, "Found %d angles(s)\n", mptr->nangles);
	}

	if ( mptr->nangles > 0 ){
		fprintf(stdout, "Copying to device\n");
	
		size_t nbytes =  4*(mptr->nangles)*sizeof(unsigned int);
		if ( cudaMalloc((void **)&(mptr->dalist),nbytes) == cudaErrorMemoryAllocation )
			sep_cuda_mem_error();
	
		cudaMemcpy(mptr->dalist, mptr->halist, nbytes, cudaMemcpyHostToDevice);
		
		// Broken
		//sep_cuda_copy_exclusion(pptr);
	}
	
}


void sep_cuda_read_dihedrals(sepcumol *mptr, const char *file, const char opt){
	const char section[] = {'[', ' ', 'd', 'i', 'h', 'e', 'd', 'r', 'a', 'l', 's', ' ', ']','\n', '\0'};
	char line[256];
	fpos_t pos_file;
	unsigned moli, a, b, c, d, type;


	if ( mptr->nmols == 0 ) {
		fprintf(stderr, "Bond section must be read before dihedral section");
		exit(EXIT_FAILURE);
	}
	
	mptr->ndihedrals = 0; 
	FILE *fptr = fopen(file, "r");
	if ( fptr == NULL ) {printf("Here\n"); sep_cuda_file_error(); printf("and here\n");}

	fptr = sep_cuda_set_file_pointer(fptr, section);
	if ( fptr == NULL ) {
		if ( opt=='v' ){
		   fprintf(stdout, "No dihedrals found\n");
	   	}
 		return;
	}		
	
	// We *must* init the pointer since it will be free no matter if the read is sucessful or not
	mptr->hdlist = (unsigned *)malloc(0);
	if ( mptr->hdlist == NULL ) sep_cuda_mem_error();
 
	do {

		fgetpos(fptr, &pos_file); 
		if ( fgets(line, 256, fptr) == NULL ) sep_cuda_file_error();
		
		if ( line[0] == '[' ) {
			break;
		}
		else {
      
			fsetpos(fptr, &pos_file); 
      
			int sc = fscanf(fptr, "%u%u%u%u%u%u\n", &moli, &a, &b, &c, &d, &type);
			if ( sc != 6 ) sep_cuda_file_error();
     
			(mptr->ndihedrals) ++;
    
			mptr->hdlist = (unsigned *)realloc((mptr->hdlist), sizeof(unsigned)*5*mptr->ndihedrals);
			if ( mptr->hdlist == NULL ) sep_cuda_mem_error();
      
			int index0 = (mptr->ndihedrals-1)*5;
      
			mptr->hdlist[index0] = a;
			mptr->hdlist[index0+1] = b;
			mptr->hdlist[index0+2] = c;
			mptr->hdlist[index0+3] = d;
			mptr->hdlist[index0+4] = type;
			
			/* Broken
			sep_cuda_set_hexclusion(pptr, a, b); sep_cuda_set_hexclusion(pptr, a, c); sep_cuda_set_hexclusion(pptr, a, d);  
			sep_cuda_set_hexclusion(pptr, b, a); sep_cuda_set_hexclusion(pptr, b, c); sep_cuda_set_hexclusion(pptr, b, d); 
			sep_cuda_set_hexclusion(pptr, c, a); sep_cuda_set_hexclusion(pptr, c, b); sep_cuda_set_hexclusion(pptr, c, d);
			sep_cuda_set_hexclusion(pptr, d, a); sep_cuda_set_hexclusion(pptr, d, b); sep_cuda_set_hexclusion(pptr, d, c);
			*/
		}
	} while ( !feof(fptr) ); 
	
	fclose(fptr);

	// Since number of mols is one bigger than the index
	mptr->ndihedralblocks = mptr->ndihedrals/SEP_CUDA_NTHREADS + 1 ;
	
	fprintf(stdout, "Succesfully read 'dihedral' section in file %s -> ", file);
    fprintf(stdout, "Found %d angles(s)\n", mptr->ndihedrals);
	
	if ( mptr->ndihedrals > 0 ){
		fprintf(stdout, "Copying to device\n");
	
		size_t nbytes =  5*(mptr->nangles)*sizeof(unsigned int);
		if ( cudaMalloc((void **)&(mptr->ddlist),nbytes) == cudaErrorMemoryAllocation )
			sep_cuda_mem_error();
	
		cudaMemcpy(mptr->ddlist, mptr->hdlist, nbytes, cudaMemcpyHostToDevice);
		
		// Broken
		//sep_cuda_copy_exclusion(pptr);
		
		cudaDeviceSynchronize();
	}
	
}

void sep_cuda_free_angles(sepcumol *mptr){

	if ( mptr->nangles==0 ) return;

	cudaFreeHost(mptr->halist);
	cudaFree(mptr->dalist);
}

void sep_cuda_free_dihedrals(sepcumol *mptr){

	if ( mptr->ndihedrals == 0 ) return; 

	cudaFreeHost(mptr->hdlist);
	cudaFree(mptr->ddlist);
}

/*
void sep_cuda_free_mols(sepcumol *mptr){
	
	sep_cuda_free_bonds(mptr);
	sep_cuda_free_angles(mptr);
	sep_cuda_free_dihedrals(mptr);
}
*/

__global__ void sep_cuda_bond_harmonic(unsigned *blist, unsigned nbonds, float3 bondspec, 
								  float4 *pos, float4 *force, float3 lbox){
	
	unsigned i = blockDim.x*blockIdx.x + threadIdx.x;
	
	if ( i<nbonds ){
		
		int type =  __float2int_rd(bondspec.z);
		unsigned offset = i*3;
		
		if ( blist[offset+2] == type ) {
			
			unsigned a = blist[offset]; unsigned b = blist[offset+1]; 
			
			float dx = pos[a].x - pos[b].x; dx = sep_cuda_mol_wrap(dx, lbox.x);
			float dy = pos[a].y - pos[b].y; dy = sep_cuda_mol_wrap(dy, lbox.y);
			float dz = pos[a].z - pos[b].z; dz = sep_cuda_mol_wrap(dz, lbox.z);

			float dist = sqrtf(dx*dx + dy*dy + dz*dz);
						
			float ft = -bondspec.x*(dist - bondspec.y)/dist;
			
			//ACHTUNG slow perhaps
			atomicAdd(&(force[a].x), ft*dx); 
			atomicAdd(&(force[a].y), ft*dy); 
			atomicAdd(&(force[a].z), ft*dz);
			
			atomicAdd(&(force[b].x), -ft*dx); 
			atomicAdd(&(force[b].y), -ft*dy); 
			atomicAdd(&(force[b].z), -ft*dz);
			
		}
	}
	
}

__global__ void sep_cuda_angle(unsigned *alist, unsigned nangles, float3 anglespec, 
								  float4 *pos, float4 *force, float3 lbox){
	
	unsigned i = blockDim.x*blockIdx.x + threadIdx.x;
	
	if ( i<nangles ){
		
		int type =  __float2int_rd(anglespec.z);
		unsigned offset = i*4;
		float cCon = cos(SEP_CUDA_PI - anglespec.y);
		
		if ( alist[offset+3] == type ) {
			
			unsigned a = alist[offset]; 
			unsigned b = alist[offset+1];
			unsigned c = alist[offset+2];
			 
			float3 dr1, dr2;
			
			dr1.x = pos[b].x - pos[a].x; dr1.x = sep_cuda_mol_wrap(dr1.x, lbox.x);
			dr1.y = pos[b].y - pos[a].y; dr1.y = sep_cuda_mol_wrap(dr1.y, lbox.y);
			dr1.z = pos[b].z - pos[a].z; dr1.z = sep_cuda_mol_wrap(dr1.z, lbox.z);

			dr2.x = pos[c].x - pos[b].x; dr2.x = sep_cuda_mol_wrap(dr2.x, lbox.x);
			dr2.y = pos[c].y - pos[b].y; dr2.y = sep_cuda_mol_wrap(dr2.y, lbox.y);
			dr2.z = pos[c].z - pos[b].z; dr2.z = sep_cuda_mol_wrap(dr2.z, lbox.z);
	
			float c11 = sep_cuda_mol_dot(dr1, dr1);
			float c12 = sep_cuda_mol_dot(dr1, dr2);
			float c22 = sep_cuda_mol_dot(dr2, dr2);
      
			float cD = sqrtf(c11*c22); float cc = c12/cD; 

			float f = -anglespec.x*(cc - cCon);
      
			float3 f1, f2;
			
			f1.x = f*((c12/c11)*dr1.x - dr2.x)/cD; 
			f1.y = f*((c12/c11)*dr1.y - dr2.y)/cD;
			f1.z = f*((c12/c11)*dr1.z - dr2.z)/cD;
			
			f2.x = f*(dr1.x - (c12/c22)*dr2.x)/cD;
			f2.y = f*(dr1.y - (c12/c22)*dr2.y)/cD;
			f2.z = f*(dr1.z - (c12/c22)*dr2.z)/cD;
				
			//ACHTUNG slow perhaps
			atomicAdd(&(force[a].x), f1.x); 
			atomicAdd(&(force[a].y), f1.y); 
			atomicAdd(&(force[a].z), f1.z);
			
			atomicAdd(&(force[b].x), -f1.x-f2.x); 
			atomicAdd(&(force[b].y), -f1.y-f2.y); 
			atomicAdd(&(force[b].z), -f1.z-f2.z); 
			
			atomicAdd(&(force[c].x), f2.x); 
			atomicAdd(&(force[c].y), f2.y); 
			atomicAdd(&(force[c].z), f2.z);
			
		}
	}
	
}



__global__ void sep_cuda_ryckertbellemann(unsigned *dlist, unsigned ndihedrals, int type, 
										  float params1, float params2, float params3, float params4, float params5,
											float4 *pos, float4 *force, float3 lbox){
	
	unsigned i = blockDim.x*blockIdx.x + threadIdx.x;
	
	if ( i<ndihedrals ){
		
		unsigned offset = i*5;
		
		if ( dlist[offset+4] == type ) {
			
			unsigned a = dlist[offset]; 
			unsigned b = dlist[offset+1];
			unsigned c = dlist[offset+2];
			unsigned d = dlist[offset+3];
			
			float3 dr1, dr2, dr3;

			dr1.x = pos[b].x - pos[a].x; dr1.x = sep_cuda_mol_wrap(dr1.x, lbox.x);
			dr1.y = pos[b].y - pos[a].y; dr1.y = sep_cuda_mol_wrap(dr1.y, lbox.y);
			dr1.z = pos[b].z - pos[a].z; dr1.z = sep_cuda_mol_wrap(dr1.z, lbox.z);

			dr2.x = pos[c].x - pos[b].x; dr2.x = sep_cuda_mol_wrap(dr2.x, lbox.x);
			dr2.y = pos[c].y - pos[b].y; dr2.y = sep_cuda_mol_wrap(dr2.y, lbox.y);
			dr2.z = pos[c].z - pos[b].z; dr2.z = sep_cuda_mol_wrap(dr2.z, lbox.z);
	
			dr3.x = pos[d].x - pos[c].x; dr3.x = sep_cuda_mol_wrap(dr3.x, lbox.x);
			dr3.y = pos[d].y - pos[c].y; dr3.y = sep_cuda_mol_wrap(dr3.y, lbox.y);
			dr3.z = pos[d].z - pos[c].z; dr3.z = sep_cuda_mol_wrap(dr3.z, lbox.z);
	
			float c11 = sep_cuda_mol_dot(dr1, dr1); 
			float c12 = sep_cuda_mol_dot(dr1, dr2);
			float c13 = sep_cuda_mol_dot(dr1, dr3);
			float c22 = sep_cuda_mol_dot(dr2, dr2);
			float c23 = sep_cuda_mol_dot(dr2, dr3);
			float c33 = sep_cuda_mol_dot(dr3, dr3);
			 
			float cA = c13*c22 - c12*c23;
			float cB1 = c11*c22 - c12*c12;
			float cB2 = c22*c33 - c23*c23;
			float cD = sqrt(cB1*cB2); 
			float cc = cA/cD;

			float f = -(params1+(2.*params2+(3.*params3+(4.*params4+5.*params5*cc)*cc)*cc)*cc);
			
			float t1 = cA; 
			float t2 = c11*c23 - c12*c13;
			float t3 = -cB1; 
			float t4 = cB2;
			float t5 = c13*c23 - c12*c33; 
			float t6 = -cA;
			float cR1 = c12/c22; 
			float cR2 = c23/c22;
			
			float3 f1;
			f1.x = f*c22*(t1*dr1.x + t2*dr2.x + t3*dr3.x)/(cD*cB1);
			f1.y = f*c22*(t1*dr1.y + t2*dr2.y + t3*dr3.y)/(cD*cB1);
			f1.z = f*c22*(t1*dr1.z + t2*dr2.z + t3*dr3.z)/(cD*cB1);
			
			float3 f2;
			f2.x = f*c22*(t4*dr1.x + t5*dr2.x + t6*dr3.x)/(cD*cB2);
			f2.y = f*c22*(t4*dr1.y + t5*dr2.y + t6*dr3.y)/(cD*cB2);
			f2.z = f*c22*(t4*dr1.z + t5*dr2.z + t6*dr3.z)/(cD*cB2);
			
			//ACHTUNG slow perhaps
			atomicAdd(&(force[a].x), f1.x); 
			atomicAdd(&(force[a].y), f1.y);
			atomicAdd(&(force[a].z), f1.z);
			
			atomicAdd(&(force[b].x), (-(1.0 + cR1)*f1.x + cR2*f2.x)); 
			atomicAdd(&(force[b].y), (-(1.0 + cR1)*f1.y + cR2*f2.y));
			atomicAdd(&(force[b].z), (-(1.0 + cR1)*f1.z + cR2*f2.z));
			
			atomicAdd(&(force[c].x), (cR1*f1.x - (1.0 + cR2)*f2.x));
			atomicAdd(&(force[c].y), (cR1*f1.y - (1.0 + cR2)*f2.y));
			atomicAdd(&(force[c].z), (cR1*f1.z - (1.0 + cR2)*f2.z));
			
			atomicAdd(&(force[d].x), f2.x); 
			atomicAdd(&(force[d].y), f2.y);
			atomicAdd(&(force[d].z), f2.z);
			
      }
	}
		
}


__global__ void sep_cuda_calc_forceonmol(float3 *df, float3 *dfij, unsigned nmols){

	unsigned molidx = blockDim.x * blockIdx.x + threadIdx.x;

	if ( molidx < nmols ){
		df[molidx].x = df[molidx].y = df[molidx].z = 0.0f;
		for ( unsigned n=0; n<nmols; n++ ){
			df[molidx].x += dfij[nmols*molidx + n].x;
			df[molidx].y += dfij[nmols*molidx + n].y;
			df[molidx].z += dfij[nmols*molidx + n].z;
		}
	}


}

void sep_cuda_force_harmonic(sepcupart *pptr, sepcumol *mptr, int type, float ks, float lbond){
	int nb = mptr->nbondblocks; 
	int nt = pptr->nthreads;
	
	// Notice the change in sequence
	float3 bondinfo = make_float3(ks, lbond, type);
	
	sep_cuda_bond_harmonic<<<nb,nt>>>
		(mptr->dblist, mptr->nbonds, bondinfo, pptr->dx, pptr->df, pptr->lbox);
		
	cudaDeviceSynchronize();
}


void sep_cuda_force_angle(sepcupart *pptr, sepcumol *mptr, int type, float ktheta, float angle0){
	int nb = mptr->nangleblocks; 
	int nt = pptr->nthreads;
	
	// Notice the change in sequence
	float3 angleinfo = make_float3(ktheta, angle0, type);
	
	sep_cuda_angle<<<nb,nt>>>
		(mptr->dalist, mptr->nangles, angleinfo, pptr->dx, pptr->df, pptr->lbox);
		
	cudaDeviceSynchronize();
}

void sep_cuda_force_dihedral(sepcupart *pptr, sepcumol *mptr, int type, float params[6]){
	int nb = mptr->ndihedralblocks; 
	int nt = pptr->nthreads;
	
	
	sep_cuda_ryckertbellemann<<<nb, nt>>>
		(mptr->ddlist, mptr->ndihedrals, type, 
		 params[1], params[2], params[3], params[4], params[5], pptr->dx, pptr->df, pptr->lbox);

	cudaDeviceSynchronize();

}
/*
float sep_cuda_mol_translate_tobox(float x, float L){
	int nL;

	if ( x>0 ){
		nL = floor(x/L);
	   	x = x - nL*L;
	}
	else if ( x<0 ){
		nL = ceil(x/L);
		x = x + nL*L;
	}
	
	return x;
}	
*/
void sep_cuda_mol_calc_cmprop(sepcupart *pptr, sepcumol *mptr){

	unsigned nmols = mptr->nmols;
	float3 lbox = pptr->sptr->lbox;
	float3 dr;

	sep_cuda_copy(pptr, 'x', 'h'); 
	sep_cuda_copy(pptr, 'v', 'h'); 

	double *x = (double *)malloc(sizeof(double)*SEP_CUDA_MAXNUAU);
	double *y = (double *)malloc(sizeof(double)*SEP_CUDA_MAXNUAU);
	double *z = (double *)malloc(sizeof(double)*SEP_CUDA_MAXNUAU);
	
	if ( x==NULL || y==NULL || z==NULL ) sep_cuda_mem_error();

	for ( unsigned m=0; m<nmols; m++) {
        // First particle in molecule
		int i = mptr->alist[m*SEP_CUDA_MAXNUAU]; 	

		if ( i > -1 ){
			
			x[0] = pptr->hx[i].x; y[0] = pptr->hx[i].y; z[0] = pptr->hx[i].z;	
		
			int n=1; 
			while ( mptr->alist[m*SEP_CUDA_MAXNUAU+n] != -1 ){
				int aidx = mptr->alist[m*SEP_CUDA_MAXNUAU+n];
				
				dr.x = pptr->hx[aidx].x - pptr->hx[aidx-1].x; dr.x = sep_cuda_wrap_host(dr.x, lbox.x);
				dr.y = pptr->hx[aidx].y - pptr->hx[aidx-1].y; dr.y = sep_cuda_wrap_host(dr.y, lbox.y);
				dr.z = pptr->hx[aidx].z - pptr->hx[aidx-1].z; dr.z = sep_cuda_wrap_host(dr.z, lbox.z);

				x[n] = x[n-1] + dr.x; y[n] = y[n-1] + dr.y; z[n] = z[n-1]+dr.z;			

				n++;
			}

			mptr->hx[m].x=0.0f;  mptr->hx[m].y=0.0f; mptr->hx[m].z=0.0f;
			mptr->hv[m].x=0.0f;  mptr->hv[m].y=0.0f; mptr->hv[m].z=0.0f;
			mptr->masses[m] = 0.0f; 

			for ( int ni=0; ni<n; ni++ ){
				double mass = pptr->hx[ni+i].w;
				mptr->hx[m].x +=  mass*x[ni]; mptr->hx[m].y += mass*y[ni]; mptr->hx[m].z += mass*z[ni];
				
				mptr->hv[m].x += mass*pptr->hv[ni+i].x; 
				mptr->hv[m].y += mass*pptr->hv[ni+i].y; 
				mptr->hv[m].z += mass*pptr->hv[ni+i].z;

				mptr->masses[m] += mass;
			}
		
			mptr->hx[m].x = mptr->hx[m].x/mptr->masses[m]; 
			mptr->hx[m].x = sep_cuda_periodic_host(mptr->hx[m].x, lbox.x);
			
			mptr->hx[m].y = mptr->hx[m].y/mptr->masses[m]; 
			mptr->hx[m].y = sep_cuda_periodic_host(mptr->hx[m].y, lbox.y);

			mptr->hx[m].z = mptr->hx[m].z/mptr->masses[m]; 
			mptr->hx[m].z = sep_cuda_periodic_host(mptr->hx[m].z, lbox.z);

			mptr->hv[m].x = mptr->hv[m].x/mptr->masses[m];
			mptr->hv[m].y = mptr->hv[m].y/mptr->masses[m];
			mptr->hv[m].z = mptr->hv[m].z/mptr->masses[m];
		}
	}

	free(x); free(y); free(z); 

	pptr->sptr->cmflag = true;
}

void sep_cuda_mol_calc_dipoles(sepcupart *pptr, sepcumol *mptr){

	if ( !pptr->sptr->cmflag )
		sep_cuda_mol_calc_cmprop(pptr, mptr); // Cm prop. done and saved on host

	const unsigned nmols = mptr->nmols;
	float3 lbox = pptr->sptr->lbox;

	for ( unsigned m=0; m<nmols; m++) {
    	double sumz = 0.0;
    	double rpos[3]={0.0, 0.0, 0.0}, rneg[3]={0.0, 0.0, 0.0};
    	int pos_counter=0, neg_counter=0;

		int n=0; 
		do {
        	// First particle in molecule
			int i = mptr->alist[m*SEP_CUDA_MAXNUAU+n]; 	
			if ( i<0 ) break;

			double offset[3]= {0.0, 0.0, 0.0};
			if ( fabs(pptr->hx[i].x - mptr->hx[m].x) > 0.5*lbox.x ){
	  			if ( pptr->hx[i].x - mptr->hx[m].x > 0.0 ) offset[0] = -lbox.x;
	  			else offset[0] = lbox.x;
			}
      		if ( fabs(pptr->hx[i].y - mptr->hx[m].y) > 0.5*lbox.y ){
	  			if ( pptr->hx[i].y - mptr->hx[m].y > 0.0 ) offset[1] = -lbox.y;
	  			else offset[1] = lbox.y;
			}
      		if ( fabs(pptr->hx[i].z - mptr->hx[m].z) > 0.5*lbox.z ){
	  			if ( pptr->hx[i].z - mptr->hx[m].z > 0.0 ) offset[2] = -lbox.z;
	  			else offset[2] = lbox.z;
			}
      		
      		if ( pptr->hv[i].w > 0.0 ){
				rpos[0] += pptr->hx[i].x + offset[0];
				rpos[1] += pptr->hx[i].y + offset[1];
				rpos[2] += pptr->hx[i].z + offset[2];
				sumz += pptr->hv[i].w;
				pos_counter++;
      		}		
      		else if ( pptr->hv[i].w < 0.0 ){
				rneg[0] += pptr->hx[i].x + offset[0];
				rneg[1] += pptr->hx[i].y + offset[1];
				rneg[2] += pptr->hx[i].z + offset[2];

				neg_counter++;
		    }
		  	n++;
		} while (1); // end do
	 	
		double d[3]={0.0};
 		if ( neg_counter > 0 && pos_counter > 0 )
      		for ( int k=0; k<3; k++ ) d[k] = rpos[k]/pos_counter - rneg[k]/neg_counter;

		mptr->hpel[m].x = sumz*d[0]; mptr->hpel[m].y = sumz*d[1]; mptr->hpel[m].z = sumz*d[2];

	} // end nmol
	

}

void sep_cuda_mol_calc_molpress(double *P, sepcupart *pptr, sepcumol *mptr){

	if ( !pptr->sptr->molprop ) {
		fprintf(stderr, "Mol. properties flag not set to 'on' - mol. pressure not calculated\n");
		return;
	}

	if ( !pptr->sptr->cmflag )
		sep_cuda_mol_calc_cmprop(pptr, mptr); // Cm prop. done and saved on host
	
	sep_cuda_copy(pptr, 'm', 'h');        // Copy molecular forces to host *SLOW*

	unsigned nmols = mptr->nmols;

	double Pkin[3][3]={0.0f};
	for ( unsigned m=0; m<nmols; m++ ){
		Pkin[0][0] += mptr->masses[m]*mptr->hv[m].x*mptr->hv[m].x;
		Pkin[0][1] += mptr->masses[m]*mptr->hv[m].x*mptr->hv[m].y;
		Pkin[0][2] += mptr->masses[m]*mptr->hv[m].x*mptr->hv[m].z;

		Pkin[1][0] += mptr->masses[m]*mptr->hv[m].y*mptr->hv[m].x;
		Pkin[1][1] += mptr->masses[m]*mptr->hv[m].y*mptr->hv[m].y;
		Pkin[1][2] += mptr->masses[m]*mptr->hv[m].y*mptr->hv[m].z;

		Pkin[2][0] += mptr->masses[m]*mptr->hv[m].z*mptr->hv[m].x;
		Pkin[2][1] += mptr->masses[m]*mptr->hv[m].z*mptr->hv[m].y;
		Pkin[2][2] += mptr->masses[m]*mptr->hv[m].z*mptr->hv[m].z;
	}


	double Ppot[3][3]={0.0f};
	for ( unsigned m=0; m<nmols-1; m++ ){ 	
		for ( unsigned n=m+1; n<nmols; n++ ){
			float x = mptr->hx[m].x - mptr->hx[n].x; x = sep_cuda_wrap_host(x, pptr->sptr->lbox.x);
			float y = mptr->hx[m].y - mptr->hx[n].y; y = sep_cuda_wrap_host(y, pptr->sptr->lbox.y);
			float z = mptr->hx[m].z - mptr->hx[n].z; z = sep_cuda_wrap_host(z, pptr->sptr->lbox.z);

			Ppot[0][0] += x*mptr->hfij[m*nmols+n].x; Ppot[0][1] += x*mptr->hfij[m*nmols+n].y; Ppot[0][2] += x*mptr->hfij[m*nmols+n].z; 
			Ppot[1][0] += y*mptr->hfij[m*nmols+n].x; Ppot[1][1] += y*mptr->hfij[m*nmols+n].y; Ppot[1][2] += y*mptr->hfij[m*nmols+n].z;
			Ppot[2][0] += z*mptr->hfij[m*nmols+n].x; Ppot[2][1] += z*mptr->hfij[m*nmols+n].y; Ppot[2][2] += z*mptr->hfij[m*nmols+n].z;
		}
	}
	

	double ivolume = 1.0/(pptr->sptr->lbox.x*pptr->sptr->lbox.y*pptr->sptr->lbox.z);
	P[0] = (Pkin[0][0]+Ppot[0][0])*ivolume;  P[1] = (Pkin[0][1]+Ppot[0][1])*ivolume; P[2] = (Pkin[0][2]+Ppot[0][2])*ivolume;
	P[3] = (Pkin[1][0]+Ppot[1][0])*ivolume;  P[4] = (Pkin[1][1]+Ppot[1][1])*ivolume; P[5] = (Pkin[1][2]+Ppot[1][2])*ivolume;
	P[6] = (Pkin[2][0]+Ppot[2][0])*ivolume;  P[7] = (Pkin[2][1]+Ppot[2][1])*ivolume; P[8] = (Pkin[2][2]+Ppot[2][2])*ivolume;

}

double sep_cuda_mol_calc_avdipole(sepcumol *mptr){

	const unsigned nmols = mptr->nmols;

	double mu=0.0f;
	for (unsigned m=0; m<nmols; m++ )
		mu += sqrt(sep_cuda_dot_host(mptr->hpel[m]));

	return mu/nmols;

}

void sep_cuda_mol_calc_forceonmol(sepcupart *pptr, sepcumol *mptr){
	
	const int nb = mptr->nmols/pptr->nthreads+1; 
	const int nt = pptr->nthreads;

	sep_cuda_calc_forceonmol<<<nb, nt>>>(mptr->df, mptr->dfij, mptr->nmols);

	cudaDeviceSynchronize();
}

/*
void sep_cuda_mol_calc_forceonmol(sepcupart *pptr, sepcumol *mptr){
	unsigned n,m;

	const unsigned nmols = mptr->nmols;
	sep_cuda_copy(pptr, 'm', 'h');        // Copy molecular forces to host

	for ( m=0; m<nmols; m++ ){
		mptr->hf[m].x = mptr->hf[m].y = mptr->hf[m].z = 0.0f;
		for ( n=0; n<nmols; n++ ) {
			mptr->hf[m].x += mptr->hfij[n*nmols + m].x;
			mptr->hf[m].y += mptr->hfij[n*nmols + m].y;
			mptr->hf[m].z += mptr->hfij[n*nmols + m].z;
		}
	}

}
*/
