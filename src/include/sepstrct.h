/* 
* sepstrct.h - This file is a part of the sep-library 
*
* Copyright (C) 2008 Jesper Schmidt Hansen 
* 
* License: GPL - see COPYING for copying conditions.
* There is ABSOLUTELY NO WARRANTY, not even for MERCHANTIBILITY or
* FITNESS FOR A PARTICULAR PURPOSE.
*
*/


#ifndef __SEPSTRCT_H__
#define __SEPSTRCT_H__

#include "sepdef.h"
#include <stdbool.h>

/**
 * The basic particle structure - containing all relevant 
 * information about the particle such as position, velocity, mass etc.
 */
typedef struct {
  double x[3];          /**< Position */
  double v[3];          /**< Velocity */
  double f[3];          /**< Force */
  double a[3];          /**< Acceleration */  
  double m;             /**< Mass */
  char type;            /**< Particle type */
  double z;             /**< Point charge */

  int *neighb;          /**< Neighbour list */

  int cross_neighb[3];  /**< Boundary crossings - for auto neighb update */
  int crossings[3];     /**< Boundary crossings (not reset after neighb list update)*/

  int molindex;         /**< Molecule index (what molecule does this atom belong to) */  
  int bond[SEP_BOND];   /**< Index of bonded particles */
  
  double sigma;         /**< Hardsphere diameter */
  int *collid;          /**< List of collision partners */ 
  double *colltime;     /**< List of collisions times */ 

  double ldiff;          /**< Langevin diff */

  double xtrue[3];       /**< The true position without boundary wrapping */    

  double x0[3];         /**< Initial positions */    
  double xp[3];         /**< Previous positions */
  double xn[3];         /**< Position at last neighborlist update */

  double px[3];         /**< Predicted postion */
  double pv[3];         /**< Predicted velocity */
  double pa[3];         /**< Predicted acceleration */
						      			 
} seppart;

typedef seppart sepatom;


/** 
 * System intra molecular structures - contains basic information about bonds, angles etc,
 * which is used in force calculations and sampling
 */ 
typedef struct {

  unsigned num_mols; /**< Number of molecules in the system */
  unsigned max_nuau; /**< Maximum particles (united atomic units) per molecule */
	
  int flag_bonds, flag_angles, flag_dihedrals; /**< Interaction flags */

  unsigned num_bonds;    /**< Total number of bonds */
  unsigned *blist;       /**< Bond list: (the two bonded part. indicies + bond type)*num_bonds */
  unsigned num_btypes;   /**< Number of bond types */
   
  unsigned num_angles;   /**< Total number of angles */
  unsigned *alist;       /**< Angle list: (the three atom indicies + angle type)*num_angles */
  unsigned num_atypes;   /**< Number of angle types */
  
  unsigned num_dihedrals; /**< Total number of dihedrals */
  unsigned *dlist;        /**< Diheadral list: (the four atom indicies + dihedral type)*num_dihedrals */
  unsigned num_dtypes;    /**< Number of diehdral types */

  double *blengths;       /**< Array of bond lengths */
  double *angles;         /**< Array of angles */
  double *dihedrals;      /**< Array of the dihedral angles */

  unsigned flag_Fij;      /**< Flag for molecular force calculations */  
  float ***Fij;           /**< Mol force-force array (nmol x nmol x 3 ) - for pressure tensor */
  float ***Fiajb;         /**< Force-force array (npart x nmol x 3 ) - for couple tensor */
} sepmolinfo;



/**
 * Defines the system structure. Contains information about system dimensions, time, etc
 */
typedef struct {
  // System definition
  long int npart;       /**< Number of particles in system */
  double length[3];     /**< Simulation box dimensions */
  double volume;        /**< Simulation box volume */ 
  
  // Integrator specifics
  int intgr_type;       /**< Integrator type */
  double dt;            /**< Integrator step */
  double tnow;          /**< Time now */ 
  unsigned ndof;        /**< Number of degrees of freedom */

  // Force stuff
  double cf;               /**< Maximum particle cut-off */  
  double lsubbox[3];       /**< Length of subboxes */ 
  int nsubbox[3];          /**< Number of subboxes */
  double skin;             /**< Skin...! */
  unsigned neighb_update;  /**<  Neighbourlist update */
  unsigned neighb_flag;    /**< Flag for neighblist update */
  unsigned nupdate_neighb; /**< Number of nieghbourlist updates */  
  
  bool omp_flag;           /**< Whether we do parallel sims or not */
  unsigned int nthreads;   /**< Number of threads */
  
  // Molecular info - intra mol. forces
  sepmolinfo *molptr;      /**< In order to access the mol. information */ 

} sep3D;

typedef sep3D sepsys;

/**
 * Defines the molecular structure. Information about the mass, centre-of-mass velocity etc  
 */
typedef struct {
  double m;              /**< Mass */ 
  double x[3], xtrue[3]; /**< Molecular cm - wrap + no wrap */
  double v[3];           /**< Molecular cm velocity */ 

  unsigned nuau;         /**< Number of particles (united atomic units) in molecule */   
  int *index;            /**< Indicies of constituting particles */  

  double ete[3];         /**< End-to-end vector */
  double re2;            /**< End-to-end distance squared */
  double rg;             /**< Ratio of gyration */
  double S[3];           /**< Orientation order parameter */

  double s[3];           /**< Intrinsic angular momentum */
  double inertia[3][3];  /**< Moment of inertia tensor */
  double w[3];           /**< Instrinsic angular velocity */
  int method_w;          /**< Method for determining w (intertia may be rank deficient) */
  double pel[3];         /**< Electrical dipole moment */
  
  char type;             /**< Molecules type */

  unsigned nbonds;       /**< For SHAKE */ 
  double *blength;       /**< For SHAKE */
  int shake_flag;        /**< For SHAKE */

} sepmol;


/**
 * Structure for return values. 
 */ 
typedef struct {
  double etot; /**< System total energy */
  double ekin; /**< System kinetic energy */
  double epot; /**< System potential energy */ 
  double ecoul; /**< System electrostatic energy */
  double sumv2; /** Sum of velocities squared */   

  // Atomic pressure tensor
  double P[3][3];     /**< Atomic pressure tensor */
  double kin_P[3][3]; /**< Kinetic part of atomic pressure tensor */
  double pot_P[3][3]; /**< Potential part of atomic pressure tensor */               
  double p;           /**< Normal pressure  */
  
  // Molecular pressure tensor
  double P_mol[3][3];     /**< Molecular pressure tensor */
  double kin_P_mol[3][3]; /**< Kinetic part of molecular pressure tensor */
  double pot_P_mol[3][3]; /**< Potential part of molecular pressure tensor */
  double p_mol;           /**< Molecular normal pressure */
  
  // Pressure contributions DPD
  double pot_P_conservative[3][3];
  double pot_P_random[3][3];
  double pot_P_dissipative[3][3];
  double pot_P_bond[3][3];
  
  // Couple tensor
  double pot_T_mol[3][3];
  double kin_T_mol[3][3];
  double T_mol[3][3];
  double t_mol;

} sepret;

#endif
