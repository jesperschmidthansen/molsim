
#ifndef __SEPMOL_H__
#define __SEPMOL_H__

#include "sepmisc.h"
#include "sepstrct.h"
#include <stdbool.h>


/**
 * Reads the molecular topology file and sets up the infrastructure for
 * simulations of molecular systems.
 * @param aptr Pointer to the seplib particle structure
 * @param file File string 
 * @param sys Pointer to seplib system structure
 * @param opt Verbose (v) or quiet (q) option
 */
void sep_read_topology_file(sepatom *aptr, const char *file, sepsys *sys,  char opt); 


/**
 * Calculates the forces acting on the particles due to harmonic bonds of specific
 * type
 * @param aptr Pointer to the seplib particle structure
 * @param type Bond type specifier
 * @patam lbond Zero force bond length
 * @param ks Spring constant
 * @param sys Pointer to seplib system structure 
 * @param ret Pointer to seplib return structure
 */
void sep_stretch_harmonic(sepatom *aptr, int type, const double lbond,
			  const double ks, sepsys *sys, sepret *ret);

/**
 * Calculates the force acting on the particles due to finite extension non-linear
 * elastic (FENE) bonds. 
 * @param aptr Pointer to seplib particle structure
 * @param type Bond type specifier
 * @param R0 Zero force distance
 * @param K Force amplitude
 * @param sys Pointer to seplib system structure
 * @param ret Pointer to seplib return structure
 */
void sep_FENE(sepatom *aptr, int type, const double R0, const double K,
	      sepsys *sys, sepret *ret);

/**
 * Calculates the forces acting on the particles due to cosine squared 
 * angle potential. V = k/2*[cos(angle0)-cos(angle)]^2
 * @param ptr Pointer to seplib particle structure
 * @param type Angle type specifier
 * @param angle0 Zero force angle
 * @param k Force amplitude  
 * @param sys Pointer to the seplib system structure 
 * @param ret Pointer to the seplib return structure 
 */
void sep_angle_cossq(sepatom *ptr, int type, const double angle0,
		     const double k, sepsys *sys, sepret *ret);

/**
 * Calculates the forces acting on the particles due to harmonic angle potential
 * V=k/2*(angle0-angle)^2
 * @param ptr Pointer to the seplib particle structure
 * @param type Angle type specifier
 * @param angle0 Zero force angle
 * @param k Force amplitude
 * @param sys Pointer to the seplib system structure
 * @param ret Pointer to the seplib return structure
 */
void sep_angle_harmonic(sepatom *ptr, int type, 
			const double angle0, const double k, 
			sepsys *sys, sepret *ret);

/**
 * Calculates the forces acting on the particles due the Ryckaert-Belleman potential
 * V = g0 + g1*cos(theta) + ... + g5*cos(theta)^5. Implementation from Rapaport.
 * @param ptr Pointer to the seplib particle structure
 * @param type Type specifier
 * @param g Series coefficients
 * @param sys Pointer to the seplib system structure
 * @param ret Pointer to the seplib return structure
 */
void sep_torsion_Ryckaert(sepatom *ptr, int type, const double g[6],
			  sepsys *sys, sepret *ret);

/**
 * Allocates memory for the seplib molecule structure which is useful when
 * evaluating molecular properties during simulation
 * @param atom Pointer to the seplib particle structure
 * @param sys Pointer to the seplib system structure
 * @return Point to the seplib molecule structure
 */ 
sepmol *sep_init_mol(sepatom *atom, sepsys *sys);

/**
 * Frees memory for the seplib molecule structure array.
 * @param ptr Pointer to the seplib molecule structure
 * @param sys Pointer to the seplib system structure
 */
void sep_free_mol(sepmol *ptr, sepsys *sys);


/**
 * Evaluate the molecular centre-of-mass for each molecule. Value stored in the 
 * x molecule structure  member, e.g. mol[10],x[0], mol[10],x[1], mol[10].x[2]
 * @param ptr Pointer to the seplib particle structure
 * @param mol Pointer to the seplib molecule structure
 * @param sys Pointer to the seplib system structure 
 */
void sep_mol_cm(seppart *ptr, sepmol *mol, sepsys *sys);

/**
 * Evaluate the molecular centre-of-mass velocity for each
 * molecule. Value stored in the v molecule structure member,
 * e.g. mol[10].v[0], mol[10].v[1], mol[10].v[2] 
 * @param ptr Pointer to the seplib particle structure 
 * @param mol Pointer to the seplib molecule structure 
 * @param sys Pointer to the seplib system structure
 */
void sep_mol_velcm(seppart *atom, sepmol *mol, sepsys *sys);


/**
 * Evaluates the distance between atoms with indicies 'a' and 'b' in the molecule.
 * Values are stored in ete molecule structure member
 * @param  atoms Pointer to the seplib particle structure
 * @param mols Pointer the seplib molecule structure
 * @param type Molecule type
 * @param a First atom index  
 * @param b Second atom index
 * @param sys Seplib system structure
 */
void sep_mol_ete(seppart *ptr, sepmol *mols, char type,
		 unsigned a, unsigned b, sepsys sys);

/**
 * Evaluates the molecular moment of inertia and angular momentum. Values stored in the 
 * inertia and s molecular structure members. If 'safe' is false the angular velocity is 
 * evaluated; in case the inertia tensor is singular the average of the trace of the 
 * inertia tensor is used. 
 * @param atom Pointer to the seplib particle structure
 * @param mol Pointer to the seplib molecular structure
 * @param sys seplib system structure
 * @param safe Boolean option: if set to false safe angular velocity calculuated, 
 * if true angular velocity not calculated
 */
void sep_mol_spin(seppart *atom, sepmol *mol, sepsys *sys, bool safe);

/**
 * Evaluates the molecular dipole moment. Molecules must be neutrally charged.
 * Values are stored in 'pel' molecule structure member.
 * @param atom Pointer to the seplib particle structure
 * @parm mol Pointer to the seplib molecule structure
 * @param sys Pointer to the seplib system structure
 */
void sep_mol_dipoles(seppart *atom, sepmol *mol, sepsys *sys);

/**
 * Evaluates the molecular pressure tensor. RECOMMENDATION: Call sep_mol_pressure_tensor 
 * instead.
 * @param atoms Pointer to the seplib particle structure
 * @param mols Pointer to the seplib molecule structure
 * @param ret Pointer to the seplib return structure
 * @param sys Pointer to the seplib system structure
 */
void sep_eval_mol_pressure_tensor(seppart *atoms, sepmol *mols,
				  sepret *ret, sepsys *sys);

/**
 * Evaluates the molecular couple tensor. 
 * @param atoms Pointer to the seplib particle structure
 * @param mols Pointer to the seplib molecule structure
 * @param ret Pointer to the seplib return structure
 * @param sys Pointer to the seplib system structure
 */
void sep_eval_mol_couple_tensor(sepatom *atoms, sepmol *mols,
				sepret *ret, sepsys *sys);

/**
 * Evaluates the average bond length for molecules of a specific type. Useful for
 * quick testing
 * @param type Molecule type
 * @param sys Pointer to the seplib molecule structure
 * @return The average bond length
 */
double sep_average_bondlengths(int type, sepsys *sys);




#ifndef DOXYGEN_SKIP


void sep_save_mol_config(sepatom *atoms, sepmol *mol, 
			 double time, const char *file, sepsys sys);


double sep_calc_mol_temp(sepatom *atoms, sepmol *mols, sepsys sys);

void sep_mol_eval_xtrue(seppart *ptr, sepmol *mol, sepsys sys);

int sep_count_mol_type(sepmol *mols, const char type, sepsys sys);


void sep_write_molconf(seppart *atoms, sepmol *mols, 
		       const char *file, sepsys sys);
void sep_mol_write_config(seppart *atoms, sepmol *mols, sepsys sys);


FILE *sep_set_file_pointer(FILE *fptr, const char *section);

void sep_read_bonds_top(sepatom *aptr, sepmolinfo *ptr, const char *file,
			int npart, char opt);
void sep_read_angles_top(sepatom *aptr, sepmolinfo *ptr, const char *file, int npart, char opt);
void sep_read_dihedrals_top(sepatom *aptr, sepmolinfo *ptr, const char *file, int npart, char opt);

void sep_free_bonds(sepmolinfo *ptr);
void sep_free_angles(sepmolinfo *ptr);
void sep_free_dihedrals(sepmolinfo *ptr);

#endif

#endif
