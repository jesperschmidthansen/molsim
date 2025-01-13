
classdef thermostat < handle

	properties 
		t, ntypes;
		temperature;
		tauQ; 
		xi;
	end

	methods 
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

		function nosehoover(this, atoms, dt=0.005)

			ekin = ms_nosehoover(atoms.f, this.xi, atoms.v, this.t, atoms.natoms, atoms.t, atoms.m);
					
			this.xi = this.xi + dt/this.tauQ*(0.66667*ekin - this.ntypes*this.temperature);
		end		
	end
end
