
classdef integrator < handle

	properties (Access=public)
		# time step
		dt;
		# step index	
		sidx;
	end

	methods

		function this = integrator(dt=0.005)
			
			this.dt = dt;
			this.sidx = 0;
		end

		function [ekin Pkin]= lf(this, atoms, interactions)
			
			[ekin Pkin] = ms_leapfrog(atoms.v, atoms.r, atoms.f, atoms.m, atoms.boxcross, atoms.lbox, atoms.natoms, this.dt);					
			
			interactions.first_call = true;	
			
			this.sidx++;
		
		end

		function ekin = dpd(this, atoms, interactions, lambda)

			ekin = ms_verlet_dpd(atoms.r, atoms.v, atoms.f, atoms.pa, atoms.pv, atoms.m, ... 
								atoms.boxcross, atoms.lbox, atoms.natoms, this.dt, this.sidx, lambda);		
			
			interactions.first_call = true;	
			
			this.sidx++;
		end		
			
	end

end
