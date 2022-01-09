/* 
* sepinit.h - This file is a part of the sep-library 
*
* Copyright (C) 2008 Jesper Schmidt Hansen 
* 
* License: GPL - see COPYING for copying conditions.
* There is ABSOLUTELY NO WARRANTY, not even for MERCHANTIBILITY or
* FITNESS FOR A PARTICULAR PURPOSE.
*
* Contact: schmidt@zigzak.net
*/


#ifndef __SEPINIT_H__
#define __SEPINIT_H__

#include <stdio.h>
#include "sepdef.h"
#include "sepstrct.h"
#include "sepmisc.h"
#include "sepmol.h"

/**
 * Allocates memory for the seplib particle structure
 * @param npart Number of particles in the system
 * @param nneighb Maximum number of neighbours per particle. Recommendation: Use SEP_NUM_NEIGHB 
 * @return Pointer to array of the seplib particle structure
 */
seppart *sep_init(size_t npart, size_t nneighb);

/**
 * Frees the memory allocated for the seplib particle structures
 * @param ptr Pointer the array of seplib particle structure
 * @param npart Number of particles in the system
 */
void sep_close(seppart *ptr, size_t npart);

/**
 * Allocates memory for the seplib particle structure. Initializes the particle 
 * positions, velocities, mass, and charge from the input file (xyz-format). Number of 
 * neighbours is set to SEP_NUM_NEIGHB. 
 * @param lbox Pointer to double array of length 3. Stores the sim box lengths.
 * @param npart Pointer to int; stores the total number in system.
 * @param file String specifying the file
 * @param verbose 'q' for silent, 'v' for verbose mode
 * @return Pointer to array of seppart structures
 */
seppart *sep_init_xyz(double *lbox, int *npart, const char *file, char verbose);

/**
 * Sets the initial particle positions on a primitive cubic lattice. Must be called 
 * after initialization of the particle structure array and system structure. 
 * @param ptr Pointer to the seplib particle structure
 * @param sys The seplib system structure 
 */
void sep_set_lattice(seppart *ptr, sepsys sys);

/**
 * Sets the particle velocities corresponding to a kinetic temperature.
 * @param ptr Pointer to the seplib particle structure
 * @param temp The corresponding kinetic temperature
 * @param sys The seplib system structure
 */
void sep_set_vel(seppart*ptr, double temp, sepsys sys);

/**
 * Initialization of the seplib system structure. 
 * @param lengthx Simulation box length x-direction  
 * @param lengthy Simulation box length y-direction
 * @param lengthz Simulation box length z-direction
 * @param maxsyscf The maximum cutoff for the simulation. (Maximum of all potential functions.) 
 * @param dt Integrator step size
 * @param npart Number of particles in the system
 * @param update Specifying the force calculation algorithm. Valid options are: (i) SEP_BRUTE - Brute force npart-squared method. Can be used for very small systems. (ii) SEP_NEIGHBLIST  - Neighbour list method. Can be used for small systems or long ranged interactions. Is depreciated and will likely be removed in future releases. (iii) SEP_LLIST_NEIGBLIST - Linked-list + neighbourlist. Large systems and recommended. Can only be used if maxsyscf + skin < 3*lengthbox. Error message issued if this condition is not met.
 * @return The seplib system structure
 */      
sepsys sep_sys_setup(double lengthx, double lengthy, double lengthz, 
		     double maxsyscf, double dt, size_t npart, size_t update);

/**
 * Frees the memory allocated to the seplib system structure
 * @param ptr Pointer to the system structure 
 */
void sep_free_sys(sepsys *ptr);


#ifndef DOXYGEN_SKIP

//void sep_randomize_lattice(seppart *ptr, sepsys sys);
void sep_set_vel_seed(seppart*ptr, double temp, unsigned int seed,
                      sepsys sys);
void sep_set_vel_type(seppart*ptr, char type, double temp, unsigned int seed,
		      sepsys sys);

#endif
#endif
