#
# dihedrals is a class in the molsim package 
# It constains information about dihedral angles between atoms
# 
# Examples: See package examples/ folder
# 
# User-relevant properties
# - matrix: pidx, coeffs 
# - vectors: dtypes 
# - scalar: ndihedrals
#


classdef dihedrals < handle

	properties (Access=public)
		pidx;   # Array with particle indices
 		dtypes; # Dihedral types 
		coeffs; # Coefficients in Ryckaert-Belleman model
		ndihedrals; 
	end

	methods 
		
		## Usage: a = dihedrals(ndihedrals)
		## 
		## Returns an instance of the dihedrals class object   
		function this = dihedrals(ndihedrals)
			this.pidx = zeros(ndihedrals, 4);
			this.dtypes = zeros(ndihedrals, 1);
			this.coeffs = zeros(ndihedrals, 6);
			this.ndihedrals = ndihedrals;
		end

		## Usage: epot = ryckbell(atoms, dihedral type)
		## 
		## Calculates the force acting on the atoms in a dihedral angle of specified type. 
		## The angle potential is the Ryckaert-Bellmann potential
		## Returns the total potential energy. See Rapaport "The Art of Molecular Dynamics Simulation"
		## Cambridge University Press. 	
		##
		## Example:
		## >> epot = sim.dihedrals.ryckbell(sim.atoms, 0)
		function epot = ryckbell(this, atoms, dtype)
	
			epot = ms_dihedrals_ryckaert(atoms.f, atoms.r, this.ndihedrals, dtype, this.dtypes, this.coeffs, 
												this.pidx, atoms.lbox, atoms.natoms);
	 							
		end

	end	


end

