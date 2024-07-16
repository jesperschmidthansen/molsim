
classdef interactions < handle

	properties 
		first_call, first_call_simulation;
		max_cutoff, skin;
		neighb_updates;
	endproperties

	methods

		function this = interactions(max_cutoff=2.5, skin=0.5)
			this.first_call = true;
			this.first_call_simulation = true; 
			this.max_cutoff = max_cutoff; this.skin = skin;
			this.neighb_updates = 0;
		end

		function iteration_start(this, atoms)
	
			if this.first_call
				atoms.f = zeros(atoms.natoms, 3);
				dr2 = calcdr2(atoms.r, atoms.r0, atoms.boxcross, atoms.lbox, atoms.natoms);
				if this.first_call_simulation || dr2 > this.skin*this.skin 
					neighblist(atoms.nblist, atoms.r, atoms.r0, atoms.lbox, this.max_cutoff, this.skin, atoms.natoms);
					this.neighb_updates ++;
				end
				this.first_call = false;
				this.first_call_simulation = false;
			end
		end
	
		function epot = lj(this, atoms, ptypes, ljparams)
					
			this.iteration_start(atoms);

			epot = lj(atoms.f, ptypes, ljparams, atoms.r, atoms.t, atoms.nblist, atoms.lbox, atoms.natoms);

		end	

	endmethods 


end
