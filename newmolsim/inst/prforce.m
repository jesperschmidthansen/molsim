# 
# prforce is a class in the molsim-package. 
# It contains methods to calculate the forces between atoms. 
# 
# Examples: See package examples/ folder
#
# User-relevant class properties
# - Scalars: max_cutoff, skin
#
# Examples: See package examples/ folder
#
classdef prforce < handle

	properties (Access=public) 
		max_cutoff, skin;
		neighb_updates;
		first_call; 
		first_call_simulation;
	end

	methods
		
		## Usage: prf = prforce();
		##        prf = prforce(maximum cut off, skin);
		##
		## Returns an instance of the prforce class object.  
		function this = prforce(max_cutoff=2.5, skin=0.5)
			this.first_call = true;
			this.first_call_simulation = true; 
			this.max_cutoff = max_cutoff; this.skin = skin;
			this.neighb_updates = 0;
		end

	
		function iteration_start(this, atoms, cutoff)
	
			if ( this.first_call_simulation && cutoff > this.max_cutoff )
				error("There is a pair force with too large cutoff; change max_cutoff in prforce class");
			end
	
			if this.first_call # first_call set to true in integrator
				atoms.f = zeros(atoms.natoms, 3); 
				dr2 = ms_calcdr2(atoms.r, atoms.r0, atoms.lbox, atoms.natoms);
				if this.first_call_simulation || dr2 > this.skin*this.skin 
					ms_neighblist(atoms.nblist, atoms.r, atoms.r0, atoms.lbox, this.max_cutoff, 
									this.skin, atoms.natoms, atoms.exclude);
					this.neighb_updates ++;
				end
				this.first_call = false;
				this.first_call_simulation = false;
			end
		end

		## Usage: [epot Pconf] = lj(atoms, atom types, parameters)
		##
		## Calculates the forces acting between atoms using the Lennard-Jones pair interaction
		## potential. Returns the total potentual  energy and configurational pressure contribution
		## from the interactions. Interactions can be exluded, use 'help molsim.setexclusion' for
		## more details.
		##
		## atoms types are given by a string, e.g. "AA"
		## parameters is a 4-vector [cut-off, sigma, epsilon, aw] 
		## 
		## Example: 
		## >> epot = sim.prforce.lj(sim.atoms, "AB", [2.5, 1.0, 1.0, 1.0]);
		function [epot Pconf] = lj(this, atoms, ptypes, ljparams)
			this.iteration_start(atoms, ljparams(1));
			[epot Pconf] = ms_lj(atoms.f, ptypes, ljparams, atoms.r, atoms.t, atoms.nblist, atoms.lbox, atoms.natoms);
		end	

		## Usage: [epot Pconf] = sf(atoms, cut-off)
		##
		## Calculating the Coulomb forces acting between charges. The method is the shifted-force
		## method. Should not be used for confined systems, and the interaction potential must be 
		## sufficiently large.See The Journ. of Phys. Chem. vol. 166, p 5738 (2012).
		## 
		## Example:
		## >> sim.prforce.max_cutoff = 3.5;
		## >> sim.prforce.sf(sim.atoms, 3.5) 
		function [epot Pconf] = sf(this, atoms, cutoff)
			this.iteration_start(atoms, cutoff);
			[epot Pconf] = ms_sf(atoms.f, atoms.r, atoms.q, atoms.nblist, atoms.lbox, atoms.natoms, cutoff); 	
		end

		function epot = dpd(this, atoms, ptypes, params, temperature)
			this.iteration_start(atoms, params(1));
			epot = ms_dpd(atoms.f, atoms.r, atoms.v, ptypes, params, temperature, atoms.lbox, atoms.t, atoms.nblist, atoms.natoms);
		end
	
	end 


end
