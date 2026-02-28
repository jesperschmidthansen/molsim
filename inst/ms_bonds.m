#
# bonds is a class in the molsim package 
# It constains information about bonds between atoms
# 
# Examples: See package examples/ folder
# 
# User-relevant properties
# - matrix: pidx
# - vectors: springs, l0, btypes 
# - scalar: nbonds
#
classdef ms_bonds < handle

	properties (Access=public)
		pidx;   # Array with particle indices
 		btypes;  # Bond type 
		springs, l0; # Spring constant and zero force lengtih
		nbonds; 
	end

	methods 

		## Usage: b = bonds(nbonds)
		## 
		## Returns an instance of the bonds class object   
		function this = ms_bonds(nbonds)
			this.pidx = zeros(nbonds, 2);
			this.btypes = zeros(nbonds, 1);
			this.springs = zeros(nbonds, 1);
			this.l0 = zeros(nbonds, 1);
			this.nbonds = nbonds;
		end

		## Usage: epot = harmonic(atoms, bond type)
		##
		## Calculates the force acting on the atoms in a bond of specified type. 
		## The bond potential is a simple harmonic potential
		## Returns the total potential energy 	
		##
		## Example:
		## >> epot = sim.bonds.harmonic(sim.atoms, 0)
		function epot = harmonic(this, atoms, btype)
			epot = ms_harmonic(atoms.f, atoms.r, this.nbonds, btype, 
								this.btypes, this.springs, this.l0, this.pidx, atoms.lbox, atoms.natoms);	
		end
	
	end	


end

