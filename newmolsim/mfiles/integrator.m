
classdef integrator < handle

	properties (Access=public)
		dt;
	end

	methods

		function this = integrator(dt=0.005)
			
			this.dt = dt;

		end

		function ekin = step(this, atoms)
							
			ekin = leapfrog(atoms.v, atoms.r, atoms.f, atoms.lbox, atoms.natoms, this.dt);					
	
		end

	end

end
