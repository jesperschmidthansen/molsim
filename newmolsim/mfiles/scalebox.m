
function scalebox(atoms, dens0, prefac=0.999)

	densnow = atoms.natoms/atoms.volume();

	if densnow < dens0
		atoms.lbox = prefac*atoms.lbox
	else 
		atoms.lbox = atoms.lbox./prefac
	end	

end
