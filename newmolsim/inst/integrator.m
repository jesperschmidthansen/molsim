#
# integrator is a class in the molsim-package
# It contains relevant properties and methods for integration of the 
# equation of motion (Newton's second law basically)
#
# User-relevant properties
# - Scalar dt: Integrator time step
#
classdef integrator < handle

	properties (Access=public)
		# time step
		dt;
		# step index	
		sidx;
	end

	methods
		
		## Usage: intgr = integrator();
		##        intgr = integrator(time step);  
		##
		## Returns instant of integrator	
		function this = integrator(dt=0.005)
			this.dt = dt;
			this.sidx = 0;
		end

		## Usage: [ekin Pkin] = lf(atoms, interactions);
		## Leap-frog integrator
		##
		## Returns kinetic energy and kinetic part of the pressure 
		function [ekin Pkin]= lf(this, atoms, interactions)
			
			[ekin Pkin] = ms_leapfrog(atoms.v, atoms.r, atoms.f, atoms.m, atoms.bxcrs, 
										atoms.lbox, atoms.natoms, this.dt);					
			
			interactions.first_call = true;	
			this.sidx++;
		
		end

		## Usage: ekin = dpd(atoms, interactions);
		## Verlet dpd integrator - see Groot & Warren, JCP, 107:4423 (1997)
		##
		## Returns kinetic energy  
		function ekin = dpd(this, atoms, interactions, lambda)

			ekin = ms_verlet_dpd(atoms.r, atoms.v, atoms.f, atoms.pa, atoms.pv, atoms.m, ... 
								atoms.bxcrs, atoms.lbox, atoms.natoms, this.dt, this.sidx, lambda);		
			
			interactions.first_call = true;	
			this.sidx++;
		end		
			
	end

end
