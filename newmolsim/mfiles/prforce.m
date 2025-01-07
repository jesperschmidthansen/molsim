
classdef prforce < handle

	properties 
		first_call, first_call_simulation;
		max_cutoff, skin;
		neighb_updates;
		nthreads;
	end

	methods

		function this = prforce(max_cutoff=2.5, skin=0.5, nthreads=4)
			this.first_call = true;
			this.first_call_simulation = true; 
			this.max_cutoff = max_cutoff; this.skin = skin;
			this.neighb_updates = 0;
			this.nthreads = nthreads;
		
			ms_setomp(this.nthreads);
		end

		function setomp(this, ntreads)
			this.nthreads = nthreads;
			ms_setomp(this.nthreads);
		end

		function iteration_start(this, atoms, cutoff)
	
			if ( this.first_call_simulation && cutoff > this.max_cutoff )
				error("There is a pair force with too large cutoff; change max_cutoff in prforce class");
			end
	
			if this.first_call
				atoms.f = zeros(atoms.natoms, 3);
				dr2 = ms_calcdr2(atoms.r, atoms.r0, atoms.boxcross, atoms.lbox, atoms.natoms);
				if this.first_call_simulation || dr2 > this.skin*this.skin 
					ms_neighblist(atoms.nblist, atoms.r, atoms.r0, atoms.lbox, this.max_cutoff, 
										this.skin, atoms.natoms, atoms.exclude);
					this.neighb_updates ++;
				end
				this.first_call = false;
				this.first_call_simulation = false;
			end
		end
	
		function [epot Pconf] = lj(this, atoms, ptypes, ljparams)
					
			this.iteration_start(atoms, ljparams(1));

			[epot Pconf] = ms_lj(atoms.f, ptypes, ljparams, atoms.r, atoms.t, atoms.nblist, atoms.lbox, atoms.natoms);

		end	

		function [epot Pconf] = sf(this, atoms, cutoff)

			this.iteration_start(atoms, cutoff);
			
			[epot Pconf] = ms_sf(atoms.f, atoms.r, atoms.q, atoms.nblist, atoms.lbox, atoms.natoms, cutoff); 	
		
		end

	end 


end
