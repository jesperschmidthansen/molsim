
classdef integrator < handle

	properties (Access=public)
		dt;
	end

	methods

		function this = integrator(dt=0.005)
			
			this.dt = dt;

		end

		function ekin = step(this, atoms, interactions)
							
			ekin = ms_leapfrog(atoms.v, atoms.r, atoms.f, atoms.m, atoms.boxcross, atoms.lbox, atoms.natoms, this.dt);					
			interactions.first_call = true;	
		end

	end

end
