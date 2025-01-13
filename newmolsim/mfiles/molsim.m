
classdef molsim < handle
	
	properties
		# Classes
		atoms;
		
		integrator;
		thermostat;

		pairforce;
		bonds;
		angles;
		dihedrals;
		
		# Simulation system 	
		natoms;
		lbox;
		temperature;

		# Misc. properties
		nthreads;
	end

	methods
		
		function this = molsim(fnameOrSize, boxLengths, temperature)
			
			if nargin==1
				this.atoms = atoms(fnameOrSize); 
			elseif nargin==3
				this.atoms = atoms(fnameOrSize, boxLengths, temperature);
			end

			this.natoms = this.atoms.natoms; 
			this.lbox = this.atoms.lbox;
			this.nthreads = 4;

			this.integrator = integrator();
			this.pairforce = prforce(); 
			
			this.thermostat = thermostat(this.atoms);
			if nargin == 3
				this.thermostat.temperature = temperature;
			end

		end

		function set_nthreads(this, nthreads)
			ms_setomp(nthreads);
			this.nthreads = nthreads;
		end
	end

end

