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

#ifndef __SEPCOULOMB_H__
#define __SEPCOULOMB_H__

#include <stdlib.h>
#include <string.h>
#include "sepdef.h"
#include "sepstrct.h"
#include "separray.h"
#include "sepmisc.h"
#include "sepprfrc.h"

#ifdef OMP
#include <omp.h>
#endif

/**
 * Calculates the Coulomb force acting on the particles using the 
 * shifted force method. Particle charges are given in the z member 
 * of the seppart structure. See The Journal of Physical Chemistry B, 
 * 116, 5738 (2012).  
 * @param  ptr Pointer to seplib particle structure
 * @param cf Potential function cut-off
 * @param sys Pointer to seplib system structure 
 * @param retval Pointer to seplib return structure
 * @param opt Force method calculation specifier. Valid options: SEP_ALL, SEP_EXCL_BONDED,SEP_EXCL_SAME_MOL  
 */
void sep_coulomb_sf(seppart *ptr, double cf, sepsys *sys, sepret *retval, const unsigned opt);

/**
 * Calculates the Coulomb force acting on the particles using the 
 * Wolf  method. Particle charges are given in the z member 
 * of the seppart structure. 
 * @param  ptr Pointer to seplib particle structure
 * @param alpha 
 * @param cf Potential function cut-off
 * @param sys Pointer to seplib system structure 
 * @param retval Pointer to seplib return structure
 * @param opt Force method calculation specifier. Valid options: SEP_ALL, SEP_EXCL_BONDED,SEP_EXCL_SAME_MOL 
 */
void sep_coulomb_wolf(seppart *ptr, double alpha, double cf, sepsys *sys,
		      sepret *retval, const unsigned opt);


#ifndef DOXYGEN_SKIP

void sep_coulomb_sf_brute(seppart *ptr,  double cf, sepsys *sys,sepret *retval, const int opt);
void sep_coulomb_sf_neighb(seppart *ptr, double cf, sepsys *sys, sepret *retval);
void sep_coulomb_sf_neighb_omp(seppart *ptr, double cf, sepsys *sys, sepret *retval);

void sep_coulomb_wolf_brute(seppart *ptr, double alpha, double rcf, sepsys *sys,
			    sepret *ret, unsigned opt);
void sep_coulomb_wolf_neighb(seppart *ptr, double alpha, double rcf, sep3D *sys,
			     sepret *retval); 

void sep_ewald_direct(sepatom *ptr, int nrep, const sepsys sys);

#endif

#endif
