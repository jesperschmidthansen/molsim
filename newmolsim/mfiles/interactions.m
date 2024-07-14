
classdef interactions < handle


	methods

		function this = interactions()
			this;
		end

		function epot = lj(this, atoms, ptypes, ljparams, brute=false)
			
			if brute
				epot = lj(atoms.f, atoms.r, atoms.lbox, atoms.natoms);
			else	
				epot = lj(atoms.f, ptypes, ljparams, atoms.r, atoms.t, atoms.nblist, atoms.lbox, atoms.natoms);
			end

		end	

	endmethods 


end
