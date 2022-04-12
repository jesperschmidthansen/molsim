/* 
* sepmisc.h - This file is a part of the sep-library 
*
* Copyright (C) 2008 Jesper Schmidt Hansen 
* 
* License: GPL - see COPYING for copying conditions.
* There is ABSOLUTELY NO WARRANTY, not even for MERCHANTIBILITY or
* FITNESS FOR A PARTICULAR PURPOSE.
*
* Contact: schmidt@zigzak.net
*/


#ifndef __SEPMISC_H__
#define __SEPMISC_H__


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdarg.h>
#include <string.h>
#include <stdint.h>
#include <float.h>
#include <omp.h>

#include "sepstrct.h"
#include "separray.h"
#include "sepdef.h"
#include "seputil.h"
#include "sepprfrc.h"
#include "sepret.h"
#include "sepmol.h"

/**
 * Flushes standard output
 */
#define SEP_FLUSH fflush(stdout)

/**
 * Macro for defining timer clock_t
 */ 
#define SEP_TICTOC clock_t time_then

/**
 * Set timer
 */
#define SEP_TIC time_then = clock()

/**
 * Get time now
 */
#define SEP_TOC ((int)((clock()-time_then)/CLOCKS_PER_SEC))


/**
 * Generates pseudo random uniform number
 */
#define sep_rand() ( rand()/(RAND_MAX+1.0) )

/**
 * Prints something to the screen and flushes 
 * Usefull for debugging
 */
#define sep_here(x) { printf("%d\n", x); SEP_FLUSH; }

/**
 * Squares input x
 */
#define sep_Sq(x) ( (x)*(x) )

/**
 * Evaluates the absolut value
 */ 
#define sep_Abs(x) ( (x) > 0.0  ? (x) : -(x) )

/**
 * Evaluates the mirror image in periodic bc
 */
#define sep_Wrap( x, y )                          \
 {                                                \
 if ( x > 0.5*y ) x -= y;                         \
 else if  ( x < -0.5*y ) x += y;                  \
 }

/**
 * Perform periodic bc
 */
#define sep_Periodic( x, y )                 \
 {                                           \
 if ( x > y ) x -= y;                        \
 else if  ( x < 0 ) x += y;                  \
 }


/**
 * Formatted error function
 * @param str String printed to stderr before exiting program
 */
void sep_error(char *str, ...);


/**
 * Formatted warning function
 * @param str String printed to stderr. Program continous.
 */
void sep_warning(char *str, ...);

/** 
 * Standard Lennard-Jones potential
 * @param r2 distance between i and j squared
 * @param opt 'f' for force evaluation 'u' for potential
 */
double sep_lj(double r2, char opt);

/** 
 * Shifted Lennard-Jones potential at distance 2.5
 * @param r2 distance between i and j squared
 * @param opt 'f' for force evaluation 'u' for potential
 */
double sep_lj_shift(double r2, char opt);

/**
 * Weeks-Chandler-Andersen potential
 * @param r2 distance between i and j squared
 * @param opt 'f' for force evaluation 'u' for potential
 */
double sep_wca(double r2, char opt);

/** Spring potential due to virtual lattice cite 
 *  at position x0
 * @param r2 distance between i and j squared
 * @param opt 'f' for force evaluation 'u' for potential
 */
double sep_spring_x0(double r2, char opt);

/**
 * Save particle configuration in binary format
 * @param ptr Pointer to seppart structure
 * @param npart Number of particules
 * @param file File name
 */
void sep_save(seppart *ptr, size_t npart, const char *file);


/**
 * Load particle configuration (binary format). Memory allocated.
 * @param ptr Pointer to seppart structure
 * @param nneighb Maximum number of neighbours per particle (use eg SEP_NUM_NEIGHB)
 * @param npart Number of particules
 * @param file File name
 */
void sep_load(seppart *ptr, int nneighb, size_t npart, const char *file);

/**
 * Save particle configuration in xyz format
 * @param ptr Pointer to particle structure
 * @param partnames Particle types that should be saved
 * @param file File name
 * @param Writing mode ("w" for write "a" for append)
 * @param sys sepsys structure
 */
void sep_save_xyz(seppart *ptr, const char *partnames, 
		  const char *file, char *mode, sepsys sys);

/** 
 * Relax/thermostat temperature of specific particle type
 * @param Pointer to particle structure
 * @param type Particle type
 * @param Td Desired/target temperature
 * @param tau Relaxation time (set to eg 0.1)
 * @param sys Pointer to sepsys structure
 */
void sep_relax_temp(seppart *ptr, char type, double Td, double tau, sepsys *sys);


/**
 * Berendsen barostat using atomic pressure; particle positions and box simply 
 * scaled in relaxational fashion and in accordance to desired pressure. Anisotropic 
 * compression: only z-direction is changed. See sep_berendsen_iso for isotropic compression.
 * IMPORTANT NOTE: Should not be used together with sep_compress_box() 
 * @param ptr Pointer to seplib particle structure
 * @param Pd Desired pressure
 * @param beta Strength of relaxation mechanism
 * @param ret Pointer to seplib return structure 
 * @param sys Pointer to seplib system structure
 */
void sep_berendsen(sepatom *ptr, double Pd, double beta, sepret *ret, sepsys *sys);


/**
 * Berendsen barostat using atomic pressure; particle positions and box simply 
 * scaled in relaxational fashion and in accordance to desired pressure.
 * IMPORTANT NOTE: Should not be used together with sep_compress_box() 
 * @param ptr Pointer to seplib particle structure
 * @param Pd Desired pressure
 * @param beta Strength of relaxation mechanism
 * @param ret Pointer to seplib return structure 
 * @param sys Pointer to seplib system structure
 */
void sep_berendsen_iso(sepatom *ptr, double Pd, double beta, sepret *ret, sepsys *sys);

/** 
 * Berendsen barostat using molecular pressure; molecule center of mass  
 * positions and box simply scaled in relaxational fashion 
 * and in accordance to desired pressure
 * IMPORTANT NOTE: Should not be used together with sep_compress_box() 
 * @param ptr Pointer to seplib particle structure
 * @param mol Pointer to seplib molecule structure
 * @param Pd Desired pressure
 * @param beta Strength of relaxation mechanism
 * @param ret Pointer to seplib return structure 
 * @param sys Pointer to seplib system structure
 */
void sep_berendsen_mol(sepatom *ptr, sepmol *mol, double Pd, 
		       double beta, sepret *ret, sepsys *sys);

/**
 * Tests if a particle is in some interval min to max in direction dir. 
 * Returns 1 on success.
 * @param ptr Pointer to particle structure
 * @param min Lower bound in interval 
 * @param max Upper bound in interval
 * @param i Particle index
 * @param dir Direction (0,1,2 for x, y, z)
 */
int sep_is_here(seppart *ptr, double min, double max, int i, int dir);

/**
 * Evaluates the distance and vector displacement between two particles
 * Return the distance. 
 * @param r Pointer to double array of length 3 (the displacement vector).
 * @param i Index of first particle 
 * @param j Index of second particle
 * @param sys Pointer to seplib system structure   
 */
double sep_dist_ij(double *r, seppart *ptr, int i, int j, sepsys *sys);

/**
 * Resets particle forces to zero
 * @param ptr Pointer to particle structure
 * @param sys Pointer to system structure 
 */
void sep_reset_force(seppart *ptr, sepsys *sys);


/**
 * Resets molecular force arrays 
 * @param sys Pointer to seplib system structure
 */
void sep_reset_force_mol(sepsys *sys);

/**
 * Calculates the force acting on particles from virtual site
 * x0
 * @param ptr Pointer to particle structure
 * @param type Particle type/label 
 * @param fun Function specfying the potential function
 * @param sys Pointer to system structure 
 */
void sep_force_x0(seppart *ptr, char type, double (*fun)(double, char), sepsys *sys);

/**
 * Performs second order reaction of type A+B->C+D+heat. The heat is 
 * implemented as momentum conserving thermal kinetic energy. Returns the 
 * number of reaction occurences.
 * @param ptr Pointer to particle structure
 * @param rstr Reaction string eg "ABCD"
 * @param cr Critical reaction radius (particles only react within this radius)
 * @param Pcr Reaction probability ([0;1])
 * @param Q Specifying the heat of reaction
 * @param sys Pointer to system structure
 */
unsigned int sep_reaction_2_order(seppart *ptr,  const char *rstr, 
				  double cr, double Pcr, double Q, 
				  sepsys *sys);

/**
 * Counts the number of particles of a specific type
 * @param ptr Pointer to particle structure
 * @param type Particle type to count
 * @param npart Number of particles in the system (i.e. sys.npart) 
 */
int sep_count_type(seppart *ptr,  char type, int npart);

/**
 * Sets a specific particle type to random particles
 * @param ptr Pointer to particle structure
 * @param type Type to be set eg. 'B'
 * @param numb Number of particles to label (must be smaller than npart)
 * @param npart Number of particles in the system (ie. sys.npart)
 */
void sep_set_type(seppart *ptr,  char type, int numb, int npart);

/**
 * Sets the seplib omp parallisation flag to 'on' and specifies the number
 * of threads/cores to execute in parallel.
 * @param nthreads Number of threads/cores
 * @param sys Pointer to seplib system structure
 */
void sep_set_omp(unsigned nthreads, sepsys *sys);

/**
 * Sets the additional skin/radius in the neighbourlist. Can be varied
 * for optimization. Deault value is 0.25
 * @param sys Pointer to seplib system structure
 * @param value Skin value
 */
void sep_set_skin(sepsys *sys, double value);

/**
 * Sets the particle charge. 
 * @param ptr Pointer to seplib particle structure
 * @param type The particle type/label to give charge
 * @param z Charge value
 * @param sys The seplib system structure 
 */
void sep_set_charge(seppart *ptr, char type, double z, sepsys sys);


/**
 * Sets the particle mass. 
 * @param ptr Pointer to seplib particle structure
 * @param type The particle type/label to give charge
 * @param m Mass value
 * @param sys The seplib system structure 
 */
void sep_set_mass(seppart *ptr, char type, double m, sepsys sys);

/**
 * Sets momentum of particles to zero
 * @param ptr Pointer to seplib particle structure
 * @param type Type/label of particle that is reset
 * @param sys Pointer to seplib system structure
 */
void sep_reset_momentum(seppart *ptr, const char type, sepsys *sys);

/**
 * Sets/specifies the number of degrees of freedom. This is, for example,
 * used in the Nose-Hoover .thermostat. Default is 3*npart-3  
 * @param ndof Degree of freedom
 * @param sys Pointer to seplib system structure
 */
void sep_set_ndof(size_t ndof, sepsys *sys);

/**
 * Evaluates the particle positions using the boundary crossings, hence,
 * the true positions and not the wrapped once. Results are stored in the 
 * xtrue array member of the seplib particle structure.
 * @param ptr Pointer to seplib particle structure
 * @param sys Pointer to the seplib system structure
 */
void sep_eval_xtrue(seppart *ptr, sepsys *sys);

/**
 * Sets the virtual particle positions using the current positions. (x0=x)
 * @param ptr Pointer to seplib particle strucutre
 * @param npart Number of system particles (sys.npart)
 */
void sep_set_x0(seppart *ptr, int npart);

/**
 * Compress/reduce the simulation box size if number density is smaller than
 * a specific value. All box directions are reduced the same. Should be called multiple
 * times in order to reach desired density.
 * @param ptr Pointer to seplib particle structure
 * @param rhoD Desired number density
 * @param xi Reduction parameter
 * @param sys Pointer to seplib system structure
 */
void sep_compress_box(sepatom *ptr, double rhoD, double xi, sepsys *sys);


/**
 * Compress/reduce on side of the simulation box if number density is smaller than
 * a specific value. Other directions are left as is. Should be called multiple
 * times in order to reach desired density.
 * @param ptr Pointer to seplib particle structure
 * @param rhoD Desired number density
 * @param xi Reduction parameter
 * @param dir Direction specifier (0,1,2) for (x,y,z)-directions
 * @param sys Pointer to seplib system structure
 */
void sep_compress_box_dir(sepatom *ptr, double rhoD, double xi, 
			  int dir, sepsys *sys);


/**
 * Pseudo random number from a uniform distribution in interval ]0;1[
 * SLOW  
 */
double sep_rand32(void);

/**
 * Pseudo random number from a normal distribution with zero mean and 
 * unit variance.
 * 
 */
double sep_randn(void);


#ifndef DOXYGEN_SKIP
void sep_scale_vel_index(seppart *ptr, double Td, int Id, int ndim);
void sep_scale_vel(seppart *ptr, char type, double Td, int npart, 
		   int ndim);
void sep_constraint_temp(seppart *ptr, char type, int npart, int ndim);

int sep_reaction_1_order(seppart *ptr, const char *reaction, 
			 double crp, int npart);
int sep_reaction_1_order_binomial(seppart *ptr, const char *reaction, 
			 double crp, int npart);

void sep_get_xp(seppart*ptr, double dt, int npart, int ndim);

int sep_nsubbox(double cf, double delta, double lbox);
double sep_box_length(double dens, int npart, int ndim);

double sep_eval_mom(seppart *ptr, int npart);
double sep_eval_mom_type(seppart *ptr, char type, int dir, int npart);
double sep_eval_momvec(seppart *ptr, int npart);
void sep_set_density(seppart *ptr, char typef, char typet, 
                     double densf, double denst, int npart);

void sep_rand32_init(int *r250_index, int *r521_index, 
                     unsigned long *r250_buffer, 
                     unsigned long *r521_buffer);

void sep_hs_update(sepatom *ptr, int ip, int jp, double tc, int flag,  
		   double lbox, size_t natoms, size_t ndim);
int sep_hscoll(int *ip, int *jp, double *tc,  seppart *ptr, int *list, 
	       sepsys *sys);
int sep_hs_coll(int *ip, int *jp, double *tc,  seppart *ptr, sepsys *sys);
void sep_set_ldiff(sepatom *ptr, char type, double ldiff, sepsys sys);
void sep_set_xn(seppart *ptr, int npart);
double sep_ran0(long *idum);
double sep_ran3(long *idum);

#endif




#endif




