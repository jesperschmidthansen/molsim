// Molecules + electro-statics
// Check against Wu et al. SPC/Fw water

#include "sepcuda.h"

int main(void){

	sepcupart *aptr = sep_cuda_load_xyz("start_water.xyz");
	sepcusys *sptr = sep_cuda_sys_setup(aptr);	
	sepcumol *mptr = sep_cuda_init_mol();

	sep_cuda_read_bonds(aptr, mptr, "start_water.top", 'q');
	sep_cuda_read_angles(mptr, "start_water.top", 'q');
	
	sep_cuda_set_exclusion(aptr, "molecule");
	
	float ljparam[4]={1.0,1.0,2.5,1.0};
	
	sptr->dt = 0.0005;

	int sintdip = 1000;int sintstrs = 10; int sintmisc = 20;
	sep_cuda_set_molforcecalc_on(sptr, sintstrs);

	sepcusampler_stress* stresscorr = sep_cuda_sample_stress_init(sptr, 100, 5, sptr->dt*sintstrs);
	sepcusampler_dipole* polcorr = sep_cuda_sample_dipole_init(sptr, 100, 10, sptr->dt*sintdip);
	
	FILE *fout = fopen("test.dat", "w");

	int nloops = 10000; int counter = 0; char filestr[100];
	for ( int n=0; n<nloops; n++ ){
	
		sep_cuda_reset_iteration(aptr);
		
		sep_cuda_update_neighblist(aptr, 2.9); 
	
		sep_cuda_force_lj(aptr, "OO", ljparam);
		sep_cuda_force_sf(aptr, 2.9);
		
		sep_cuda_force_harmonic(aptr, mptr, 0, 68000, 0.316);
		sep_cuda_force_angle(aptr, mptr, 0, 490 , 1.97);
		
		sep_cuda_thermostat_nh(aptr, 3.86, .1);
		sep_cuda_integrate_leapfrog(aptr);

		if ( n%2==0 ) sep_cuda_check_neighblist(aptr, sptr->skin);

		if ( n%sintdip==0 )
		 	sep_cuda_sample_dipole(polcorr, aptr, sptr, mptr);
		if ( n%sintstrs==0 )
			sep_cuda_sample_stress(stresscorr, aptr, sptr, mptr);
			
		if ( n%10000==0 ){
			sprintf(filestr, "molsim-%05d.xyz", counter);
			sep_cuda_save_xyz(aptr, filestr);
			counter ++;
		}

		if ( n%sintmisc==0 ){
			double  P[9]; float dump[3];
			sep_cuda_mol_calc_cmprop(aptr, mptr);
			sep_cuda_mol_calc_dipoles(aptr, mptr); 
			sep_cuda_mol_calc_molpress(P, aptr, mptr);
			double mu = sep_cuda_mol_calc_avdipole(mptr);
	
			sep_cuda_get_energies(aptr);
			
			printf("%f %f %f %f %f ", 
				   sptr->ekin, sptr->epot, sptr->etot, sptr->temp, sep_cuda_eval_momentum(dump, aptr));
		
			//printf("%d %f ", n, mu);
			//for ( int k=0; k<9; k++ ) printf("%f ", P[k]);
			printf("\n");
			fflush(stdout);
		}

	}


	sep_cuda_save_xyz(aptr, "test.xyz");
		
	sep_cuda_sample_stress_free(stresscorr); 
	sep_cuda_sample_dipole_free(polcorr);
	
	sep_cuda_free_memory(aptr);
	
	sep_cuda_free_bonds(mptr);
	sep_cuda_free_angles(mptr);

	return 0;
}
