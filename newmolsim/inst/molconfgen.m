# Usage: molconfgen(xyz file, top file, no. of molecules, density, seed)
# 
# Writes system configuration file and topology files from single molecules configurations
# and topologies. At the moment only one molecule type is supported.
#
# Writes files "conf.xyz", "bond.top", "angles.top" and "dihedrals.top"
#
# If not specified seed has value 42
#
# Example:  
# >> molconfgen("water.xyz", "water.top", 2000, 0.01);
function molconfgen(xyzfile, topfile, nmols, dens, seed=42)

	ms_molconfig(xyzfile, topfile, nmols, dens, seed);

	# Clear the function as this uses static variables 
	# which remain persistent between calls
	clear -f ms_molconfig 		

end

