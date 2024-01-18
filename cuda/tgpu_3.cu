// Minimal and benchmark prg
// At the moment manual neighbour list update 

#include "sepcuda.h"

int main(void){
	
	sepcupart *ptr = sep_cuda_load_xyz("start_slitpore.xyz");
	sepcusys *sptr = sep_cuda_sys_setup(ptr);
	
	sep_cuda_load_lattice_positions(ptr, "start_slitpore.xyz");
	
	int n=0; int nloops = 100000; 
	while ( n<nloops ){
		
		sep_cuda_reset_iteration(ptr);

		if ( n%10==0 )	sep_cuda_update_neighblist(ptr, 2.5);
		
		sep_cuda_force_lj(ptr);
		sep_cuda_force_lattice(ptr, 'w', 1000.0);
		
		sep_cuda_thermostat_nh(ptr, 1.0, 0.1);
		sep_cuda_integrate_leapfrog(ptr);
	
		if ( n%2==0 ) sep_cuda_check_neighblist(ptr, sptr->skin);
	
		n++;
	}
	
	sep_cuda_save_xyz(ptr, "test.xyz");
	
	sep_cuda_free_memory(ptr);
	
	return 0;
}
