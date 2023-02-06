/* 
* sepprfrc.h - This file is a part of the sep-library 
*
* Copyright (C) 2008 Jesper Schmidt Hansen 
* 
* License: GPL - see COPYING for copying conditions.
* There is ABSOLUTELY NO WARRANTY, not even for MERCHANTIBILITY or
* FITNESS FOR A PARTICULAR PURPOSE.
*
* Contact: schmidt@zigzak.net
*/

#ifndef __SEPPRFRC_H__
#define __SEPPRFRC_H__

#include <stdlib.h>
#include <string.h>
#include "sepdef.h"
#include "sepstrct.h"
#include "separray.h"
#include "sepmisc.h"

#ifdef OMP
#include <omp.h>
#endif


/**
 * Calculates the non-bonding pair force acting on the particles. 
 * If the SEP_BRUTE option is set in the system settings 
 * a brute force method is applied. If 
 * SEP_NEIGHB option is set a neighbour list force method is applied. If 
 * SEP_LLIST_NEIGHB is set linked list+neighbour list method is applied.
 * See sep_sys_setup(). 
 * @param ptr Pointer to the seplib particle structure
 * @param types String specifying the two types interacting eg "WF"
 * @param cf Interaction cutoff. Must be less than or equal the system
 * maximum cutoff
 * @param fun Function pointer specifying the interaction. 
 * First argument is the distance between the particles. Second argument must
 * be either 'f' for force evaluation 'u' for potential. See eg sep_lj_shift().
 * @param sys Pointer to the seplib system structure
 * @param retval Pointer to the seplib return structure
 * @param opt Option specifying exclusions. SEP_ALL: No exclusion. 
 * SEP_EXCL_BONDED: Exclude interactions between bonded particles. 
 * SEP_EXCL_SAME_MOL: Exclude interactions between 
 * particles in same molecule    
 */
int sep_force_pairs(seppart *ptr, const char *types, double cf,
		     double (*fun)(double, char), sepsys *sys, 
		     sepret *retval, const unsigned opt);

/**
 * Calculates the non-bonding forces using the Lennard-Jones interaction.  
 * If the SEP_BRUTE option is set in the system settings 
 * a brute force method is applied. If 
 * SEP_NEIGHB option is set a neighbour list force method is applied. If 
 * SEP_LLIST_NEIGHB is set linked list+neighbour list method is applied.
 * See sep_sys_setup(). A bit faster than sep_force_pairs()
 * @param ptr Pointer to the seplib particle structure
 * @param types String specifying the two types interacting eg "WF"
 * @param param Pointer to array specifying the LJ parameters: {cut-off, sigma, epsilon, 1/r^6 fac}
 * @param sys Pointer to the seplib system structure
 * @param retval Pointer to the seplib return structure
 * @param opt Option specifying exclusions. SEP_ALL: No exclusion. 
 * SEP_EXCL_BONDED: Exclude interactions between bonded particles. 
 * SEP_EXCL_SAME_MOL: Exclude interactions between particles in same molecule    
 */ 
void sep_force_lj(seppart *ptr, const char *types, 
		  const double *param, sepsys *sys, 
		  sepret *retval, const unsigned opt);


/**
 * Calculates the non-bonding pair force acting between dpd particles. The
 * standard linear and purely repulsive force is used.  
 * @param ptr Pointer to the seplib particle structure
 * @param types String specifying the two types interacting eg "WF"
 * @param cf Interaction cutoff. Must be less than or equal the system
 * maximum cutoff
 * @param aij Force repulsion parameter
 * @param temp_desired System temperature
 * @param sigma Force interaction parameter ( see Warren and Groot )
 * @param sys Pointer to the seplib return structure
 * @param opt Option specifying exclusions. SEP_ALL: No exclusion. 
 * SEP_EXCL_BONDED: Exclude interactions between bonded particles. 
 * SEP_EXCL_SAME_MOL: Exclude interactions between particles in same molecule     
 */
void sep_force_dpd(seppart *ptr, const char *types, 
		   const double cf, const double aij, 
		   const double temp_desired, 
		   const double sigma, sepsys *sys, sepret *retval,
		   const unsigned opt);


#ifndef DOXYGEN_SKIP
void sep_dpdforce_neighb(seppart *ptr, const char *types, 
			 const double cf, const double aij, 
			 const double temp_desired, 
			 const double sigma, sepsys *sys, sepret *retval,
			 const unsigned opt);



void sep_force_pair_brute(seppart *ptr, const char *types, double cf,
			  double (*fun)(double, char), sepsys *sys, 
			  sepret *retval, const int opt);

void sep_force_pair_neighb(seppart *ptr, const char *types, double cf,
			   double (*fun)(double, char), sepsys *sys, 
			   sepret *retval, bool parallel);


void sep_make_neighblist(seppart *ptr, sepsys *sys, const unsigned opt);


void sep_neighb_excl_same_mol(seppart *ptr, sepsys *sys);
void sep_neighb(seppart *ptr, sepsys *sys);
void sep_neighb_nonbonded(seppart *ptr, sepsys *sys);

void sep_make_neighblist_from_llist(seppart *ptr, int nneighb, 
				    int *list, sepsys *sys);
void sep_make_neighblist_from_llist_nonbonded(seppart *ptr, int nneighb, 
				   int *list, sepsys *sys);
void sep_make_neighblist_from_llist_excl_same_mol(seppart *ptr, int nneighb, 
				       int *list, sepsys *sys);


int *sep_allocate_celllist(sepsys *sys);
void sep_make_celllist(seppart *ptr, int *list, sepsys *sys);
unsigned int sep_bonded_direct(seppart *ptr, int j1, int j2);


// Lennard-Jones specific systems - a bit faster
void sep_lj_pair_brute(seppart *ptr, const char *types, 
		       const double *p, sepsys *sys, 
		       sepret *retval, const int opt);

void sep_lj_pair_neighb(seppart *ptr, const char *types,
			const double *param, sepsys *sys, 
			sepret *retval, bool parallel);


// Piotr's force  

void sep_add_interaction(sepsys *sys, char types[], double cf, 
			 double (*fun)(double,char));

void _sep_force_pairs(seppart *ptr, sepsys *sys, 
		      sepret *retval, unsigned n, 
		      const unsigned opt);

void _sep_force_pair_brute(seppart *ptr, sepsys *sys, 
			   sepret *retval, const int opt);


void _sep_force_pair_neighb(seppart *ptr,  sepsys *sys, 
			    sepret *retval) ;

void _sep_force_pair_neighb_omp(seppart *ptr, sepsys *sys,  sepret *retval);








void sep_dpdforce_brute(seppart *ptr, const char *types, 
			 const double cf, const double aij, 
			 const double temp_desired, 
			 const double sigma, sepsys *sys, sepret *retval);
#endif
#endif
