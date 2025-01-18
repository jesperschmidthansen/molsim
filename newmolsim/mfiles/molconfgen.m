
function molconfgen(xyzfile, topfile, nmols, dens, seed=42)

	ms_molconfig(xyzfile, topfile, nmols, dens, seed);

	# We must clear the function as this uses static variables 
	# which remain persistent between calls
	clear -f ms_molconfig 		

end

