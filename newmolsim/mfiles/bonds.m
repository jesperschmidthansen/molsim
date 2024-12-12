
classdef bonds < handle

	properties (Access=public)
		pidx;   # Array with particle indices
 		btypes;  # Bond type 
		springs, l0; # Spring constant and zero force lengtih
		nbonds; 
	end

	methods 

		function this = bonds(nbonds)
			this.pidx = zeros(nbonds, 2);
			this.btypes = zeros(nbonds, 1);
			this.springs = zeros(nbonds, 1);
			this.l0 = zeros(nbonds, 1);
			this.nbonds = nbonds;
		end

		function epot = harmonic(this, atoms, btype)
			epot = ms_harmonic(atoms.f, atoms.r, this.nbonds, btype, 
								this.btypes, this.springs, this.l0, this.pidx, atoms.lbox, atoms.natoms);	
		end
	
	end	
end

