
classdef angles < handle

	properties (Access=public)
		pidx;        # Array with particle indices
 		atypes;      # Angle types 
		springs, a0; # Spring constant and zero force lengtih
		nangles; 
	end

	methods 

		function this = angles(nangles)
			this.pidx = zeros(nangles, 3);
			this.atypes = zeros(nangles, 1);
			this.springs = zeros(nangles, 1);
			this.a0 = zeros(nangles, 1);
			this.nangles = nangles;
		end

		function epot = harmonic(this, atoms, atype)
		
			epot = ms_angle_harmonic(atoms.f, atoms.r, this.nangles, atype, this.atypes, this.springs,
										this.a0, this.pidx, atoms.lbox, atoms.natoms);
												
		end

	end	


end

