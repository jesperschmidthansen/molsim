
classdef interactions < handle


	methods

		function this = interactions()
			this;
		end

		function epot = lj(this, atoms)

			atoms.f = zeros(atoms.natoms, 3);
			epot = lj(atoms.f, atoms.r, atoms.lbox, atoms.natoms);

		end	

	end 


end
