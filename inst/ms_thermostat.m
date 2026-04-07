# 
# thermostat is a class in the molsim-package. 
# It contains all the relevat properties and methods for thermostating. 
# 
# Examples: See package examples/ folder
#
classdef ms_thermostat < handle

	properties
		 
		types; ntypes; 
		temperature;

		# Characterstic relaxation time for NH thermostat and relax thermostat
		tauQ; 

		# Reservior variable for Nose-Hoover
		xi;
	end

	methods 
		
		## Usage: thermostat = ms_thermostat(atom types, thermostat relaxation parameters) 
		##        thermostat = ms_thermostat(temperature)
		## 
		## Returns an instant of thermostat class object 
		function this = ms_thermostat(atomtypesOrTemperature, temperatureOrTauQ, temperatures)

			if nargin==2 
				this.ntypes = -1; this.xi = 0.01; 
				this.temperature = atomtypesOrTemperature;
				this.tauQ = temperatureOrTauQ;
			elseif nargin==3
				atypes = atomtypesOrTemperature;
				this.ntypes = length(atypes);			

				if this.ntypes != length(temperatures) || this.ntypes != length(temperatureOrTauQ)
					error("Number of atom types and relaxation times must be the same");
				end
			
				this.types = atypes;
				this.temperature = temperatures;	
				this.tauQ = temperatureOrTauQ;
				this.xi = 0.05*rand(1, this.ntypes);
			else
				error("ms_thermostat accepts one or three arguments");
			end

		end
		
		## Usage: relaxttemp(atoms, atom type, temperature, relaxation time)
		function relaxtemp(this, atoms)
			if this.ntypes == -1
				ms_relaxtemp_all(atoms.v, atoms.natoms, atoms.m, this.tauQ, this.temperature);
			else 
				for n=1:this.ntypes
					ms_relaxtemp(atoms.v, this.types(n), atoms.natoms, atoms.t, atoms.m, this.tauQ(n), this.temperature(n));
				end
			end
		end

		## Usage: nosehoover(atoms, dt)
		##
		## Add force to the atoms of specified type (see constructor) in accordance with the
		## Nose-Hoover scheme. dt is the integrator time step.
		function nosehoover(this, atoms, dt)

			if this.ntypes == -1 # We thermostate all atoms
				ekin = ms_nosehoover_all(atoms.f, this.xi, atoms.v, atoms.natoms, atoms.m);
				this.xi = this.xi + dt/this.tauQ*(2/3*ekin - atoms.natoms*this.temperature);
			else # We thermostat types WARNING: Not momentum conserving
				for n=1:this.ntypes
					numthistype = length( find(atoms.t==this.types(n)) );
					ekin = ms_nosehoover(atoms.f, this.xi(n), atoms.v, this.types(n), atoms.natoms, atoms.t, atoms.m);
					this.xi(n) = this.xi(n) + dt/this.tauQ(n)*(2./3.*ekin - numthistype*this.temperature(n));
				end
			end

		end		

	end

end
