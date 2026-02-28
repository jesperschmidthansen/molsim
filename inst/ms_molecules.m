# 
# molecules contains information about the molecules  
# 
# In future this will enable more complicated systems, eg. with 
# molecule mixture. Intra-molecular interactions should be moved here
#
classdef ms_molecules < handle
	
	properties (Access=public)
		
		# Molecule info
		nmols; 
		nuau;
		atom_idxs; 

	end

	methods
		
		## Usage: molinfo = molecules();
		## 
		## Returns an empty instance of molecules class object
		function this = ms_molecules()
			this;	
		end
	
	end	
end
