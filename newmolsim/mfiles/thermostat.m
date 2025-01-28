# 
# thermostat is a class in the molsim-package. 
# It contains all the relevat properties and methods for thermostating. 
# 
# Examples: See package examples/ folder
#
classdef thermostat < handle

	properties 
		t, ntypes;
		temperature;
		tauQ; 
		xi;
	end

	methods 
		
		## Usage: therm = thermostat(atoms, temperature, atom type, thermostat relaxation parameter) 
		##        
		## Returns an instant of thermostat class object 
		function this = thermostat(atoms, temperature=1.0, ptype='A', tauQ=10.0)
			
						
			this.temperature = temperature;
			this.t = ptype;
			this.tauQ = tauQ;
			this.xi = 0.01;
			this.ntypes = length(find( atoms.t == this.t ));

			if ( length(this.ntypes) < 1 )
				error("Number of particles of the given type is zero");
				return;
			end	

		end

		function relaxtemp(this, atoms)
			ms_relaxtemp(atoms.v, this.t, atoms.natoms, atoms.t, atoms.m, this.tauQ, this.temperature);
		end

		## Usage: nosehoover(atoms, dt)
		##
		## Add force to the atoms of specified type (see constructor) in accordance with the
		## Nose-Hoover scheme. dt is the integrator time step (default is 0.005).
		function nosehoover(this, atoms, dt=5e-3)

			ekin = ms_nosehoover(atoms.f, this.xi, atoms.v, this.t, atoms.natoms, atoms.t, atoms.m);
					
			this.xi = this.xi + dt/this.tauQ*(2/3*ekin - this.ntypes*this.temperature);
		end		

	end

end
