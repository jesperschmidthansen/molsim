
classdef thermostat < handle

	properties 
		ptype, ntypes;
		temperature;
		tauQ; 
		xi;
	end

	methods 
		function this = thermostat(atoms, temperature, ptype='A', tauQ=10.0)
			
			this.temperature = temperature;
			this.ptype = ptype;
			this.tauQ = tauQ;
			this.xi = 0.01;
			this.ntypes = length(find( atoms.t == this.ptype ));

			if ( length(this.ntypes) < 1 )
				error("Number of particles of the given type is zero");
				return;
			end	

		end

		function relaxtemp(this, atoms)
			ms_relaxtemp(atoms.v, this.ptype, atoms.natoms, atoms.t, atoms.m, this.tauQ, this.temperature);
		end

		function nosehoover(this, atoms, dt=0.005)

			ekin = ms_nosehoover(atoms.f, this.xi, atoms.v, this.ptype, atoms.natoms, atoms.t, atoms.m);
					
			this.xi = this.xi + dt/this.tauQ*(0.66667*ekin - this.ntypes*this.temperature);
		end		
	end
end
