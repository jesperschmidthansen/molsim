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
		# Characterstic relaxation time 
		tauQ; 
		# Reservior variable for Nose-Hoover
		xi;
	end

	methods 
		
		## Usage: therm = thermostat(atoms, temperature, atom type, thermostat relaxation parameter) 
		##        
		## Returns an instant of thermostat class object 
		function this = thermostat(tauQ=10.0)
			this.tauQ = tauQ;
			this.xi = 0.01;
			this.ntypes = -1;
		end
		
		function settype(this, atoms, atype)
			this.t = atype;
			this.ntypes = length(find(atoms.t == atype));
			if ( this.ntypes < 1 )
				error("Number of particles of the given type is zero");
				return;
			end	
		end

		function relaxtemp(this, atoms)
			if this.ntypes == -1
				ms_relaxtemp_all(atoms.v, atoms.natoms, atoms.m, this.tauQ, this.temperature);
			else
				ms_relaxtemp(atoms.v, this.t, atoms.natoms, atoms.t, atoms.m, this.tauQ, this.temperature);
			end
		end

		function relaxttemp(this, atoms, type, temperature, tauQ=0.01)
			ms_relaxtemp(atoms.v, type, atoms.natoms, atoms.t, atoms.m, tauQ, temperature);
		end

		## Usage: nosehoover(atoms, dt)
		##
		## Add force to the atoms of specified type (see constructor) in accordance with the
		## Nose-Hoover scheme. dt is the integrator time step.
		function nosehoover(this, atoms, dt)

			if this.ntypes == -1
				ekin = ms_nosehoover_all(atoms.f, this.xi, atoms.v, atoms.natoms, atoms.m);
				this.xi = this.xi + dt/this.tauQ*(2/3*ekin - atoms.natoms*this.temperature);
			else
				ekin = ms_nosehoover(atoms.f, this.xi, atoms.v, this.t, atoms.natoms, atoms.t, atoms.m);
				this.xi = this.xi + dt/this.tauQ*(2/3*ekin - this.ntypes*this.temperature);
			end
		
		end		

	end

end
