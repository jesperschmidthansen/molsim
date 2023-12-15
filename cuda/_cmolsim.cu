
/********************************
 * 
 * A wrapper for sepcuda
 * 
 * CUDA headers and libraries are stored different places
 * depending on platform. This wrapper will hopefully make 
 * the installation procedure less cumbersom
*
 ********************************/
 
#include "cmolsim.h"

#include "sepcuda.h"

#include <stdarg.h>
#include <string.h>
#include <stdbool.h>

sepcupart *pptr; 
sepcusys *sptr;
sepcumol *mptr;

float maxcutoff = 2.5;

int iterationnumber = 0;
int resetmomentumfreq = -1;
int neighblistcheckfreq = 2;

int ensemble = 0; // 0: nve, 1: nvt 

bool init = false, initmol = false;


void load_xyz(const char file[]){
	
	pptr = sep_cuda_load_xyz(file);
	sptr = sep_cuda_sys_setup(pptr);
	
	init = true;
}

void load_top(const char file[]){
	
	mptr = sep_cuda_init_mol();
	
	sep_cuda_read_bonds(pptr, mptr, file, 'v');
	sep_cuda_read_angles(pptr, mptr, file, 'v');
	sep_cuda_read_dihedrals(pptr, mptr, file, 'v');
	
	initmol = true;
}

void free_memory(void){
	
	if ( initmol ){
		sep_cuda_free_bonds(mptr);
		sep_cuda_free_angles(mptr);
		sep_cuda_free_dihedrals(mptr);
		
		initmol = false;
	}
	
	if ( init ) {
		sep_cuda_free_memory(pptr);

		init = false;
	}
}

void reset_iteration(void){
	
	sep_cuda_reset_iteration(pptr);

//	if ( iterationnumber%neighbupdatefreq == 0 ) 	
	sep_cuda_update_neighblist(pptr, maxcutoff);

}

void reset_momentum(int freq){
	resetmomentumfreq = freq;
}


void force_lj(const char *types, float *ljparams){
	
	sep_cuda_force_lj(pptr, types, ljparams);
	
}

void force_coulomb(float cf){
	
	sep_cuda_force_sf(pptr, cf);
	
}

void force_bond(int type, float lbond, float ks){

	sep_cuda_force_harmonic(pptr, mptr, type, ks, lbond);

}

void force_angle(int type, float angle, float kangle){

	sep_cuda_force_angle(pptr, mptr, type, kangle, angle);

}

void force_torsion(int type, float *params){
	
	sep_cuda_force_dihedral(pptr, mptr, type, params);

}

void integrate_leapfrog(void){
	
	sep_cuda_integrate_leapfrog(pptr);

	if ( iterationnumber%neighblistcheckfreq==0 )
		sep_cuda_check_neighblist(pptr, sptr->skin);  

	if ( resetmomentumfreq >= 0 && iterationnumber%resetmomentumfreq==0 )
		sep_cuda_reset_momentum(pptr);
	
	iterationnumber ++;
}

void save_xyz(const char filename[]){
	
	sep_cuda_save_xyz(pptr, filename);
	
}

void thermostat_nh(const char *type, float temp0, float mass){
	
	sep_cuda_thermostat_nh(type[0], pptr, temp0, mass);

}

void set_exlusion_molecule(const char rule[]){

	sep_cuda_set_exclusion(pptr, rule);

}

void set_timestep(float dt){

	sptr->dt = dt;

}

void set_maximum_cf(float cf){
	
	maxcutoff = cf;
}

unsigned get_npart(void){

	return sptr->npart;

}

void get_pressure(double *presspointer){
	
	double normalpress, shearpress[3];
	sep_cuda_get_pressure(&normalpress, shearpress, pptr);
	
	presspointer[0] = normalpress;
	for ( int k=1; k<4; k++ ) presspointer[k]=shearpress[k-1];
	
}

void get_energies(double *energypointer){
	
	sep_cuda_get_energies(pptr);
	
	energypointer[0] = sptr->ekin;
	energypointer[1] = sptr->epot;

}

void get_positions(double *positions){

	sep_cuda_copy(pptr, 'x', 'h'); 
	for ( unsigned n=0; n<sptr->npart; n++ ) {
		positions[n] = pptr->hx[n].x;
		positions[n+1] = pptr->hx[n].y;
		positions[n+2] = pptr->hx[n].z;
	}
}


void get_masses(double *masses){

	sep_cuda_copy(pptr, 'x', 'h'); 
	for ( unsigned n=0; n<sptr->npart; n++ ) masses[n] = pptr->hx[n].w;

}

void get_charges(double *charges){

	sep_cuda_copy(pptr, 'v', 'h'); 
	for ( unsigned n=0; n<sptr->npart; n++ ) charges[n] = pptr->hv[n].w;

}

void get_types(char *types){

	for ( unsigned n=0; n<sptr->npart; n++ ) types[n] = (char)pptr->hf[n].w;

}

void get_velocities(double *velocities){

	sep_cuda_copy(pptr, 'v', 'h'); 
	for ( unsigned n=0; n<sptr->npart; n++ ) {
		velocities[n] = pptr->hv[n].x;
		velocities[n+1] = pptr->hv[n].y;
		velocities[n+2] = pptr->hv[n].z;
	}
}


