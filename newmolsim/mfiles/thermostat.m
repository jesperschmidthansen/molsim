
classdef thermostat < handle

	properties 
		ptype;
		temperature;
		tau; 
	end

	methods 
		function this = thermostat(atoms, temperature, ptype='A', tau=0.01)
			
			this.temperature = temperature;
			this.ptype = ptype;
			this.tau = tau;

			ntypes = find( atoms.t == this.ptype );
			if ( length(ntypes) < 1 )
				error("Number of particles of the given type is zero");
				return;
			end	

		end

		function relaxtemp(this, atoms)
			relaxtemp(atoms.v, this.ptype, atoms.natoms, atoms.t, atoms.m, this.tau, this.temperature);
		end

	end
end
