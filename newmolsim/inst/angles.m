#
# angles is a class in the molsim package 
# It constains information about angles between atoms
# 
# Examples: See package examples/ folder
# 
# User-relevant properties
# - matrix: pidx
# - vectors: springs, a0, atypes 
# - scalar: nangles
#
classdef angles < handle

	properties (Access=public)
		pidx;        # Array with particle indices
 		atypes;      # Angle types 
		springs, a0; # Spring constant and zero force lengtih
		nangles; 
	end

	methods 

		## Usage: a = angles(nangles)
		## 
		## Returns an instance of the angles class object   
		function this = angles(nangles)
			this.pidx = zeros(nangles, 3);
			this.atypes = zeros(nangles, 1);
			this.springs = zeros(nangles, 1);
			this.a0 = zeros(nangles, 1);
			this.nangles = nangles;
		end

		## Usage: epot = harmonic(atoms, angle type)
		##
		## Calculates the force acting on the atoms in an angle of specified type. 
		## The angle potential is a simple harmonic potential
		## Returns the total potential energy 	
		##
		## Example:
		## >> epot = sim.angles.harmonic(sim.atoms, 0)
		function epot = harmonic(this, atoms, atype)
		
			epot = ms_angle_harmonic(atoms.f, atoms.r, this.nangles, atype, this.atypes, this.springs,
										this.a0, this.pidx, atoms.lbox, atoms.natoms);
												
		end

	end	


end

