

classdef dihedrals < handle

	properties (Access=public)
		pidx;   # Array with particle indices
 		dtypes; # Dihedral types 
		coeffs; # Coefficients in Ryckaert-Belleman model
		ndihedrals; 
	end

	methods 

		function this = dihedrals(ndihedrals)
			this.pidx = zeros(ndihedrals, 4);
			this.dtypes = zeros(ndihedrals, 1);
			this.coeffs = zeros(ndihedrals, 6);
			this.ndihedrals = ndihedrals;
		end

		function epot = ryckbell(this, atoms, dtype)
	
			epot = ms_dihedrals_ryckaert(atoms.f, atoms.r, this.ndihedrals, dtype, this.dtypes, this.coeffs, 
												this.pidx, atoms.lbox, atoms.natoms);
	 							
		end

	end	


end

