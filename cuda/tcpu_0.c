
/** @example prg1.c
 *
 * Standard NVT Lennard-Jones simulation (large system)
 *
 * Tests for Nose-Hoover thermostat and linked list + neighbour list.
 * Tests for velocity and stress autocorrelations and mean square displacement
 * samplers.  
 */

#include "sep.h"

int main(int argc, char **argv){
  sepatom *atoms;
  sepsys sys;
  sepret ret;
  double t, dt, rcut, 
    dens, lbox[3], 
    etot, sump;
  int natoms, n, nloops;


	if ( argc != 2 ) {
		fprintf(stderr, "Please provide filename\n");
		exit(EXIT_FAILURE);
	}
	
  // Setting parameter values 
  rcut = 2.5;
  nloops = 100000;
  dt = 0.005;

  // Allocating memory 
  atoms = sep_init_xyz(lbox, &natoms, argv[1], 'v');
   
  // Setting up the system
  // Setting system 
  sys = sep_sys_setup(lbox[0], lbox[1], lbox[2], 2.5, dt, natoms, SEP_LLIST_NEIGHBLIST);

  sep_set_omp(4, &sys);
  
  // Main loop 
  t=0.0; n = 0;
  while ( n<nloops ){
  
    // Reset return values 
    sep_reset_retval(&ret);

    // Reset force 
    sep_reset_force(atoms, &sys);

    // Evaluate forces acting on between part. Particle lable is 'A' as default 
    sep_force_pairs(atoms, "AA", rcut, sep_lj_shift, &sys, &ret, SEP_ALL);
      
    // Integrate particles forward in time
    sep_leapfrog(atoms, &sys, &ret);
    
    t += dt; n++;
  }
  
  sep_save_xyz(atoms, "A", "test.xyz", "w", &sys);
  // Freeing memory 
  sep_close(atoms, natoms);
  sep_free_sys(&sys);

  return 0;
} 
