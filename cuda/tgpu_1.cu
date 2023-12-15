// Sanity check prg
/*
 *  Ref. seplib. Dens 0.7513, T=2.0, <etot> = -0.82, <press>=4.31 eta_0 = 1.6
 *               Radial distribution function found in rdf_single.dat
 * 
 *  Tested 23 April            
 */
#include "sepcuda.h"


int main(int argc, char **argv){
	
	if ( argc != 2 ) {
		fprintf(stderr, "Provide ensemble option\n");
		exit(EXIT_FAILURE);
	}
	
	char ensemble[10]="nve";
	if ( atoi(argv[1])==1 ) ensemble[2]='t';
	
	printf("Ensemble is %s\n", ensemble);
			
	sepcupart *ptr = sep_cuda_load_xyz("test.xyz");
	sepcusys *sptr = sep_cuda_sys_setup(ptr);

	sepcugh *ghptr = sep_cuda_sample_gh_init(sptr, 100, 5, 10*sptr->dt);
	
	float dump[3];
	float temp0 = 2.0; 	char filestr[100];
	int n = 0; int nloops = 100000; int counter = 0;
	while ( n<nloops ){

		//if ( sep_cuda_logrem(n, 2) ){
		if ( n%1000==0 ){
			sprintf(filestr, "molsim-%05d.xyz", counter);
			sep_cuda_save_xyz(ptr, filestr);
			
			sprintf(filestr, "crossings-%05d.dat", counter);
			sep_cuda_save_crossings(ptr, filestr, n*sptr->dt);
			
			counter ++;
		}
		
		sep_cuda_reset_iteration(ptr);

		sep_cuda_update_neighblist(ptr, 2.5); 
	
		sep_cuda_force_lj(ptr);
			
		if ( atoi(argv[1])==1 )	sep_cuda_thermostat_nh('A', ptr, temp0, 0.1);	

		sep_cuda_integrate_leapfrog(ptr);
		
		if ( n%100==0 ) sep_cuda_reset_momentum(ptr);
			
		if ( n%2==0 ) sep_cuda_check_neighblist(ptr, sptr->skin);
		
		if ( n%10 ==0 ){
			sep_cuda_sample_gh(ghptr, ptr, sptr);
		}
				
		
		if ( n%1000 == 0 ){
			double normalpress, shearpress[3];
			sep_cuda_get_pressure(&normalpress, shearpress, ptr);
			sep_cuda_get_energies(ptr);
			
			printf("%f %f %f %f %f %f %f %f %f \n", 
				   sptr->ekin, sptr->epot, sptr->etot, sptr->temp, normalpress, 
				   shearpress[0], shearpress[1], shearpress[2], sep_cuda_eval_momentum(dump, ptr));
		}
		
		n++;
	}

	sep_cuda_save_xyz(ptr, "test.xyz");
	
	sep_cuda_sample_gh_free(ghptr);
	sep_cuda_free_memory(ptr);
	
	return 0;
}
