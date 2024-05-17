/* 
* sepsampler.h - This file is a part of the sep-library 
*
* Copyright (C) 2008 Jesper Schmidt Hansen 
* 
* License: GPL - see COPYING for copying conditions.
* There is ABSOLUTELY NO WARRANTY, not even for MERCHANTIBILITY or
* FITNESS FOR A PARTICULAR PURPOSE.
*
* Contact: schmidt@zigzak.net
*/


#ifndef __SEPSAMPLER_H__
#define __SEPSAMPLER_H__


#include "sep.h"
#include "complex.h"
#include <stdbool.h>

/*********************************
 *
 * The different samplers
 *
 ************************************/

// sacf (atomic stress acf -
// contributions from dihedrals and angle not included!!)
typedef struct {
  unsigned i, lvec, nsample, isample;
  double dt, dtsample;

  double *sacf;  
  double **stress;
 
} sepsacf;


// Atomic radial distribution function
typedef struct {
  int i, lvec, nsample, isample, ntypes, ncomb;
  int **hist; // lvec x npairs

  char types[256]; // Lazy...
  
} sepradial;

// Atomic msd and self-intermediate scattering
typedef struct {
  int i, lvec, nsample, isample, nk, npart, ntype;
  double dt, dtsample;
  char type;
  
  double *msd, *msdsq, *k, *time;
  complex double **Fs;
  double **prev_pos, **pos0;
  int **crossings;
  bool logflag;
  int logcounter;
  
} sepmsd;

// Molecular msd and self-intermediate scattering
typedef struct {
  int i, lvec, nsample, isample, nk, nmol, ntype;
  double dt, dtsample;
  char type;
  
  double *msd, *k;
  complex double **Fs;
  double **prev_pos, **pos0;
  int **crossings;

} sepmmsd;

// msacf (molecular stress acf)
typedef struct {
  unsigned i, lvec, nsample, isample;
  double dt, dtsample;

  double **sacf;  
  double **sstress;
  double **astress;
  
} sepmsacf;

// mcacf (couple-stress acf - always molecular)
typedef struct {
  unsigned i, lvec, nsample, isample;
  double dt, dtsample;

  double **cacf;  
  double **scouple;
  double **acouple;
  
} sepmcacf;


// vacf (atomic vel. acf)
typedef struct {
  unsigned i, lvec, nsample, isample;
  double dt, dtsample;

  double *vacf;  
  double **vels;

} sepvacf;

// mvacf (single molecular velocity acf)
typedef struct {
  unsigned i, lvec, nsample, isample;
  double dt, dtsample;

  double *vacf;  
  double **vels;

} sepmvacf;

// mavacf (single molecular angular velocity acf)
typedef struct {
  unsigned i, lvec, nsample, isample;
  double dt, dtsample;

  double *avacf;  
  double **avels;

} sepmavacf;

// gh (atomic generalized hydrodynamics sampler)
typedef struct {
  unsigned i, lvec, nsample, isample, nwave, kdir, tdir, ncalls;
  double dt, dtsample, *k;
  double avekin; 
  
  complex double **fk_tv, **fkm_tv, **c_tv;      // Transverse momentum
  complex double **fk_lv, **fkm_lv, **c_lv;      // Longitudinal momentum
  complex double **fk_rho, **fkm_rho, **c_scatt; // Density
  complex double **fk_e, **fkm_e, **c_e;         // Kin. energy

  // Cross couplings
  complex double **c_re, **c_er, **c_rj, **c_jr, **c_ej, **c_je;
  
  // Template for more options
  complex double **fk_X, **fkm_X, **c_X;

} sepgh;

// mgh (molecular generalized hydrodynamics)
typedef struct {
  unsigned i, lvec, nsample, isample, nwave, ncalls;
  bool safe;
  double dt, dtsample, *k, avekin;


  complex double **fk_tv, **fkm_tv, **c_tv;      // Transverse mometum
  complex double **fk_lv, **fkm_lv, **c_lv;      // Longit. momentum
  complex double **fk_rho, **fkm_rho, **c_scatt; // Density
  complex double **fk_e, **fkm_e, **c_e;         // Transl. kin. energy

  complex double **fk_tav, **fkm_tav, **c_tav;   // Trans. ang. momentum
  complex double **fk_lav, **fkm_lav, **c_lav;   // Longit. ang. momentum
  complex double **fk_vav, **fkm_vav, **c_vav;   // Mom-ang.mol 

  complex double **fk_dip, **fkm_dip, **c_dip;   // Dipole

  // Template for more options
  complex double **fk_X, **fkm_X, **c_X;
  
} sepmgh;


// hprof (atomic flow hydrodynamic profiler - dens, velocity, temperature)
typedef struct {
  unsigned  lvec, nsample, isample;
  char type;
  double dir, dirvel;
  
  double *svel, *momc, *dens, *temp;
  
} sepprofs;

// hprof (molecular version)
typedef struct {
  unsigned  lvec, nsample, isample;
  char type;
  int dir, dirvel, diramom;
  
  double **amom, *momc, *dens, *temp, *numb;
  double ***inertia; 

} sepmprofs;

// scatt (coherent partial scattering function)
typedef struct{
	bool init;
	int index, nsample, isample, lvec, nwave, ncomb;

	double tspan;

	complex double **scatt_p_1, **scatt_m_2;
	complex double **scatt;

	char intropt, types[2]; 
	char outfile[256];
} sepscatt;



/***************************
 *
 *   Main structure 
 *
 ***************************/ 

typedef struct {

  sepmol *molptr;

  // Atomic samplers
  int flag_sacf; sepsacf *sacf;
  int flag_vacf; sepvacf *vacf;
  int flag_gh; sepgh *gh;
  int flag_profs; sepprofs *profs;
  int flag_radial; sepradial *radial;
  int flag_msd; sepmsd *msd;
  int flag_scatt; sepscatt **scatt;

  // Molecular samplers
  int flag_msacf; sepmsacf *msacf;
  int flag_mcacf; sepmcacf *mcacf;
  int flag_mgh; sepmgh *mgh;
  int flag_mprofs; sepmprofs *mprofs;
  int flag_mvacf; sepmvacf *mvacf;
  int flag_mavacf; sepmavacf *mavacf;
  int flag_mmsd; sepmmsd *mmsd;

  long unsigned msd_counter;
} sepsampler;



/**
 * Initialize sampler 
 * @return The seplib sampler structure 
 */ 
sepsampler sep_init_sampler(void);

/**
 * Add an action to the sampler. Output is written to file in local directory.
 * @param sptr Pointer to the seplib sampler structure 
 * @param sampler Sampler specification. Supported samplers: \n \n
 * (1) "sacf" (stress auto-correlation). Additional argument(s): (i) Sample function time duration [type double]\n
 * (2) "msacf" (molecular sacf). Additional argument(s): (i) Sample function time duration [type double]\n 
 * (3) "vacf" (vel. auto-correlation). Additional argument(s): (i) Sample function time duration [type double]\n
 * (4) "mvacf" (mol. vacf). Additional argument(s): (i) Sample function time duration [type double]\n
 * (5) "gh" (hydrodynamic correlation functions). Additional argument(s): 
 * (i) Sample function time duration [type double] (ii) Number of wavevectors [type int] \n
 * (6) "mgh" (mol. hydrodynamic correlation functions). Additional argument(s): 
 * (i) (i) Sample function time duration [type double] (ii) Number of wavevectors [type int] 
 * (iii) Safe mode [type int] \n
 * (7) "msd" (mean square displacement, scattering function). Additional argument(s): 
 * (i) Sample function time duration [type double] (ii) Number of wavevectors [type int] (iii) 
 * Particle type [type char] \n
 * (8) "mmsd" (mean square displacement, scattering function). Additional argument(s): 
 * (i) Sample function time duration [type double] (ii) Number of wavevectors [type int] (iii) 
 * molecule type [type char] \n
 * (9) "profs" (hydrodynamic profiler). Additional argument(s): 
 * (i) Particle type [type char] (ii) Iteration steps between sampling [type int] \n
 * (10) "mprofs" (molecular hydrodynamic profiler). Additional argument(s): 
 * (i) Molecule type [type char] (ii) Iteration steps between sampling [type int] \n
 * (11) "scatt" (All atoms partial coherent scattering function). Additional arguments
 * (i) atom types (ii) Number of wavevector, (iii) sample time span
 */
void sep_add_sampler(sepsampler *sptr, const char *sampler, sepsys sys, int lvec, ...);

/**
 * Close the sampler and free memory
 * @param Pointer to the seplib sampler structure
 */
void sep_close_sampler(sepsampler *ptr);

/**
 * Evaluates the properties specified in the sep_add_sampler() function
 * @param pptr Pointer to the seplib particle structure 
 * @param sptr Pointer to the seplib system structure
 * @param ret Pointer to the seplib return structure
 * @param sys Seplib system structure
 * @param n Simulation iteration step
 */
void sep_sample(seppart *pptr, sepsampler *sptr, sepret *ret, sepsys sys, unsigned n);

/**
 * Setup/prepares the molecular samplers. Must be called prior to the sep_add_sampler() if molecular
 * samplers are used
 * @param sptr Pointer to the seplib sampler structure
 * @param mols pointer to the seplib molecular structure  
 */
void sep_add_mol_sampler(sepsampler *sptr, sepmol *mols);

#ifndef DOXYGEN_SKIP
// sacf (atomic stress acf)
sepsacf *sep_sacf_init(int lvec, double tsample, double dt);
void sep_sacf_sample(sepsacf *sptr, sepret *ret, sepsys sys);
void sep_sacf_close(sepsacf *ptr);


// radial (atomic radial distribution)
sepradial *sep_radial_init(int lvec, int sampleinterval, char *types);
void sep_radial_sample(sepradial *sptr, seppart *atom, sepsys sys);
void sep_radial_close(sepradial *ptr);

// Mean square displacement and inchorent scatt
// (atomic radial distribution)
sepmsd *sep_msd_init(int lvec, double tsample, int nk, char type, sepsys sys);
void sep_msd_sample(seppart *atom, sepmsd *sptr, sepsys sys);
void sep_msd_close(sepmsd *ptr);
void sep_msd_crossing(seppart *atom, sepmsd *sptr, sepsys sys);

// msacf (molecular stress acf)
sepmsacf *sep_msacf_init(int lvec, double tsample, double dt);
void sep_msacf_sample(sepmsacf *sptr, sepatom *atoms, sepmol *mols, 
		      sepret *ret, sepsys sys);
void sep_msacf_close(sepmsacf *ptr);

// vacf (atomic vel. acf)
sepvacf *sep_vacf_init(int lvec, double tsample, sepsys sys);
void sep_vacf_sample(seppart *ptr, sepvacf *vptr, sepsys sys);
void sep_vacf_close(sepvacf *ptr);

// mvacf (single molecular velocity acf)
sepmvacf *sep_mvacf_init(int lvec, double tsample, sepsys sys);
void sep_mvacf_sample(seppart *ptr, sepmvacf *vptr, sepmol *mols, sepsys sys);
void sep_mvacf_close(sepmvacf *ptr);

// gh (atomic generalized hydrodynamics)
sepgh *sep_gh_init(int lvec, double tsample, double dt, int nwave, double Ly);
void sep_gh_close(sepgh *ptr);
void sep_gh_sampler(sepatom *ptr, sepgh *gh, sepsys sys);

// mgh (molecular generalized hydrodynamics)
sepmgh *sep_mgh_init(int lvec, double tsample, double dt, int nwave, 
		     double Ly, int safe);
void sep_mgh_close(sepmgh *ptr);
void sep_mgh_sampler(sepatom *ptr, sepmol *mols, sepmgh *mgh, sepsys sys);

// hprof (atomic flow hydrodynamic profiler - dens, velocity, temperature)
sepprofs *sep_profs_init(char type, int lvec, int isample);
void sep_profs_close(sepprofs *ptr);
void sep_profs_sampler(seppart *pptr, sepprofs *ptr, sepsys sys);

// hprof (molecular version)
sepmprofs *sep_mprofs_init(char type, int lvec, int isample, int dir, int dirvel, 
			   int diramom);
void sep_mprofs_close(sepmprofs *ptr);
void sep_mprofs_sampler(seppart *pptr, sepmol *mols, sepmprofs *ptr, sepsys sys); 

// mcacf (couple-stress acf - always molecular)
sepmcacf *sep_mcacf_init(int lvec, double tsample, double dt);
void sep_mcacf_sample(sepmcacf *sptr, sepatom *atoms, sepmol *mols, 
		      sepret *ret, sepsys sys);
void sep_mcacf_close(sepmcacf *ptr);

// mavacf (single molecular angular velocity acf)
sepmavacf *sep_mavacf_init(int lvec, double tsample, sepsys sys);
void sep_mavacf_sample(seppart *ptr, sepmavacf *vptr, sepmol *mols, sepsys sys);
void sep_mavacf_close(sepmavacf *ptr);


// Molecular mean square displacement and inchorent scatt
// (atomic radial distribution)
sepmmsd *sep_mol_msd_init(int lvec, double tsample, int nk, char type, sepsys sys);
void sep_mol_msd_sample(seppart *atom, sepmol *mol, sepmmsd *sptr, sepsys sys);
void sep_mol_msd_close(sepmmsd *ptr);

// Partial coherent scattering function
sepscatt** sep_init_scatt(const char types[], int lvec, int nwave, double tspan, double dt);
sepscatt* sep_init_sscatt(const char types[], int lvec, int nwave, double tspan, const char file[]);
void sep_sscatt_sample(sepscatt *scatt, sepatom *atoms, char opt, sepsys *sys);
void sep_scatt_sample(sepscatt** scattarray, sepatom *atoms, char opt, sepsys *sys);
void sep_free_sscatt(sepscatt *scatt);
void sep_scatt_close(sepscatt** scattarray);
void sep_set_inter_scatt(sepsampler *sampler);

#endif

#endif

