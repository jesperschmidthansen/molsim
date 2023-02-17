##
## -- [npart, types, pos, vel, mass, charge, lbox]  = molsim_readxyz(filename, option)
## -- [npart, types, pos] = molsim_readxyz(filename)
## 
## Reads xyz-format configuration file
##
## Inputs: 
##
## filename - the xyz-file (string)
## option  - specific xyz format file option (default 'simple')
##           *if option is 'simple' the following format is assumed
##                type x-position y-position z-position
##
##           *if option is 'sep' the seplib xyz format is assumed.
##
## Outputs:
##
## npart, types, pos - numb. particles (scalar), particle types (vector), position (npart x 3 matrix) 
## 
## vel, mass, charge lbox - velocities (npart x 3 matrix), masses (vector), charges (vector), box dimensions (vector)  
##

function [npart, types, pos, vel, mass, charge, lbox] = molsim_readxyz(filename, opt='simple')

	if ( nargin < 1 || nargin > 2  )
		print_usage("Number of input arguments is wrong - see the help text");
	end
	
	fxyz = fopen(filename, 'r');

	npart = fscanf(fxyz, "%d\n", "C");
	
	switch (opt)
		case 'sep'
			[lbox(1), lbox(2), lbox(3)] = fscanf(fxyz, "%f %f %f\n", "C");
			
			pos = zeros(npart,3); vel=zeros(npart,3); mass=zeros(npart,1); charge=zeros(npart,1);
			for n=1:npart
				[types(n), pos(n,1), pos(n,2), pos(n,3), vel(n,1) vel(n,2) vel(n,3), mass(n), charge(n)]= ...
				fscanf(fxyz, "%c%f%f%f%f%f%f%f%f\n", "C");
			end

		case 'simple'
		
			lbox = vel = mass = charge = 0;
			dummy = fgets(fxyz, 1024);
			
			pos = zeros(npart,3);
			for n=1:npart
				[types(n), pos(n,1), pos(n,2), pos(n,3)] = fscanf(fxyz, "%c%f%f%f\n", "C");
			end
			
		otherwise
		
			usage("Note valid option - see the help text");
	end
	
	fclose(fxyz);
end
