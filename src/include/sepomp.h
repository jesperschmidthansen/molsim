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

#ifndef __SEPOMP_H__
#define __SEPOMP_H__

#include "sepdef.h"
#include "sepstrct.h"
#include "separray.h"
#include "sepmisc.h"
#include "sepmol.h"


/**
 * Generates the particle neighbourlists. Must be called before usage of 
 * sep_omp_pairs() or sep_omp_pairs_lj() in connected with Openmp Model II.
 * @param ptr Pointer to the seplib particle structure
 * @param cf Cutoff - must be smaller than or equal to maximum cutoff
 * @param sys Pointer to the seplib system structure
 * @param opt Exclusion option - valid macros are SEP_ALL, SEP_NEIGHB_EXCL_BONDED,
 * or SEP_NEIGHB_EXCL_SAME_MOL
 */
void sep_omp_make_neighb(seppart *ptr,double cf,
			 sepsys *sys, const unsigned opt);

/**
 * Calculates the forces acting between the particles using the Openmp Model II. 
 * The forces are 
 * added to ftot. Neighbourlist must be updated using sep_omp_make_neighb()
 * @param ftot Double pointer to double (an array) Forces are added to this array
 * @param ptr Pointer to the seplib particle structure
 * @param types String specifying the particles that interact
 * @param cf Cutoff - must be smaller than or equal to the maximum cutodd
 * @param fun Pointer to the potential interaction function
 * @param sys Pointer to the seplib system structure   
 */
void sep_omp_pairs(double **ftot, const seppart *ptr,
		   const char *types, double cf,
		   double (*fun)(double, char), const sepsys *sys);

/**
 * Calculates the forces acting between the particles using a
 * Lennard-Jones potential using the Openmp Model II. 
 * The forces are added to ftot. Neighbourlist 
 * must be updated using sep_omp_make_neighb()  
 * @param ftot Double pointer to double (an array) Forces are added to this array
 * @param ptr Pointer to the seplib particle structure 
 * @param types String specifying the particles that interact 
 * @param param Array specifying the LJ parameters: {cf, sigma, epsilon}. cf
 * must be smaller than or equal to the maximum cutoff 
 * @param sys Pointer to the seplib system structure
 */
void sep_omp_pairs_lj(double **ftot, const seppart *ptr,
		      const char *types, const double *param, const sepsys *sys);


/**
 * Calculates the electrostatic forces acting between the particles using the
 * shifted Coulomb potential with the Openmp Model II. The forces are added to ftot.
 * Neighbourlist must be updated using sep_omp_make_neighb()    
 * @param ftot Double pointer to double (an array) Forces are added to this array
 * @param ptr Pointer to the seplib particle structure 
 * @param cf Cutoff - must be smaller than or equal to maximum cutoff
 * @param sys Pointer to the seplib system structure
 */
void sep_omp_coulomb(double **ftot, seppart *ptr, double cf, sepsys *sys);

/**
 * Calculates the forces acting between particles shearing a bond with 
 * the Openmp Model II . The bond potential is the standard harmonic bond. 
 * @param ftot Double pointer to double (an array) Forces are added to this array
 * @param aptr Pointer to the seplib particle structure
 * @param type Bond type specifier
 * @param lbond Zero force bond length
 * @param ks Bond spring constant
 * @param sys Pointer to the seplib system structure  
 */
void sep_omp_bond(double **ftot, seppart *aptr, int type, 
		  const double lbond, const double ks, sepsys *sys);

/**
 * Calculates the forces acting on particles defining an angle using the 
 * Openmp Model II. The potential is the cosine squared potential.
 * @param ftot Double pointer to double (an array) Forces are added to this array
 * @param ptr  Pointer to the seplib particle structure
 * @param type Angle type specifier
 * @param angle0 Zero force angle
 * @param k Force amplitude constant
 * @param sys Pointer to the seplib system structure
 */
void sep_omp_angle(double **ftot, seppart *ptr, int type, 
		   const double angle0, const double k, sepsys *sys);

/**
 * Calculates the forces acting on particles defining a diheadral using the 
 * Openmp Model II. The potential is the Ryckaert-Bellemann potential.
 * @param ftot Double pointer to double (an array) Forces are added to this array
 * @param ptr  Pointer to the seplib particle structure
 * @param type Diheadral type specifier
 * @param g Potential parameters
 * @param sys Pointer to the seplib system structure
 */
void sep_omp_torsion(double **ftot, seppart *ptr, int type, 
		     const double g[6],  sepsys *sys);


/**
 * Calculates the DPD forces acting on particles using the Openmp
 * Model II. No exlusion 
 * @param ftot Double pointer to double (an array) Forces are added to this array
 * @param ptr Pointer to the seplib particle structure 
 * @param types Particle type specifier (String)
 * @param cf Cutoff - must be smaller than or equal to maximum cutoff
 * @param aij Repulsion parameter
 * @param temp_desired Desired temperature
 * @param sigma Force interaction parameter ( see Warren and Groot )
 * @param sys Pointer to the seplib system structure
 */
void sep_omp_dpd_pairs(double **f_tot, seppart *ptr, const char *types, 
		       const double cf, const double aij, 
		       const double temp_desired, 
		       const double sigma, sepsys *sys);

#endif
