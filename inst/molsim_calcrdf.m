
##
## Usage 
##  [radius, radial] = molsim_calcrdf(types, npoints)
##
##  types defaults to "AA", npoints to 100
##

function [radius, radial] = molsim_calcrdf(types="AA", npoints=100)
	
	if ( !exist("molsim_rdf") || !exist("molsim_readxyz") )
		error("molsim_calcrdf dependent on molsim_rdf and molsim_readxyz");
	endif
	
	counter = 0;
	while (1)
		filename = sprintf("molsim-%05d.xyz", counter);
		
		if ( !exist(filename) ) 
			break; 
		endif
		
		[npart, parttypes, pos, vel, mass, charge, lbox] = molsim_readxyz(filename, 'sep'); 
	
		[rdf, r] = molsim_rdf(types, pos, parttypes, lbox, npoints);
	
		if ( counter==0 )
			rdfsum = rdf; radius =  r; 
		else
			rdfsum = rdf + rdfsum;
		endif
		
		counter++;
	endwhile
	
	if ( counter == 0 )
		warning("No files found");
		radial = radius = 0;
	else 
		radial = rdfsum./counter; 
	endif
	
endfunction
