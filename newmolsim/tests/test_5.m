
function [dist _angles dihedrals] = test_5()

	addpath("../mfiles/"); addpath("../mex/");

	niter = 1e4;

	sim = molsim("mol.xyz"); 
	sim.integrator.dt = 1.0e-3;

	# Molecular definitions: Each mol is four atoms long. 
	# Linear with three bonds, two angles, one dihedral
	# All this will later be automated using input files
	if rem(sim.natoms,4) != 0
		error("test_5 failed - input file error");
	end
	nmols = int32(sim.natoms/4);

	# Bonds 	
	nbonds = nmols*3; sim.bonds = bonds(nbonds);
	pidx = 0;
	for n=1:nmols
		for m=1:3
			pidx++;
			sim.bonds.pidx(3*(n-1)+m, :) = [pidx, pidx+1]; 
		end
		pidx = pidx + 1;
	end
	sim.bonds.springs = 500*ones(nbonds,1); sim.bonds.l0 = 1.1*ones(nbonds,1); 

	# Angles
	nangles = nmols*2; sim.angles = angles(nangles);
	pidx = 0;
	for n=1:nmols
		for m=1:2
			pidx++;
			sim.angles.pidx(2*(n-1)+m, :) = [pidx, pidx+1, pidx+2];
		end
		pidx = pidx + 2;
	end
	sim.angles.springs = 100*ones(nangles,1); sim.angles.a0 = 2.1*ones(nangles, 1);

	# Dihedrals
	ndihedrals = nmols; sim.dihedrals = dihedrals(ndihedrals);
	pidx = 0;
	for n=1:nmols
		pidx++;
		sim.dihedrals.pidx(n, :) = [pidx, pidx+1, pidx+2, pidx+3];
		pidx = pidx + 3;
	end
	sim.dihedrals.coeffs = zeros(ndihedrals,6); sim.dihedrals.coeffs(:, 2) = -30.0;

	# Exclude pair forces in molecule
	for n=1:4:sim.natoms
		sim.atoms.exclude(n, 1:3) = [n+1, n+2, n+3]; 
		sim.atoms.exclude(n+1, 1:3) = [n, n+2, n+3];
		sim.atoms.exclude(n+2, 1:3) = [n, n+1, n+3];
		sim.atoms.exclude(n+3, 1:3) = [n, n+1, n+2];
	end	

	ekin = zeros(1, niter); epot = zeros(1, niter); 
	dist = zeros(1, niter); _angles = zeros(1, niter); dihedrals = zeros(1, niter);
	
	sim.atoms.resetmom();

	for n=1:niter
		epot(n) = sim.pairforce.lj(sim.atoms, "AA", [2.5, 1.0, 1.0, 1.0]);   
		epot(n) = epot(n) + sim.bonds.harmonic(sim.atoms, 0);
		epot(n) = epot(n) + sim.angles.harmonic(sim.atoms, 0);
		epot(n) = epot(n) + sim.dihedrals.ryckbell(sim.atoms, 0);

		ekin(n) = sim.integrator.step(sim.atoms, sim.pairforce);
		
		dist(n) = sim.atoms.getdist(1,2);
		_angles(n) = sim.atoms.getangle(1,2,3);
		dihedrals(n) = sim.atoms.getdihedral(1,2,3,4);
	end
	
	ekin = ekin./sim.natoms; epot = epot./sim.natoms; etot = epot + ekin;
	spnb = niter/sim.pairforce.neighb_updates; 
	mom = sim.atoms.getmom();

	printf("test_5 output:\n");
	printf("Etot: %1.3e +/- %1.3e \n", mean(etot), std(etot));
	printf("Ekin %1.3e  Momentum: %e %e %e\n", mean(ekin), mom(1), mom(2), mom(3));
	printf("Steps per build  %.3f\n", spnb);
	printf("Av. bond length %.3f, eq. bond length set to %.3f \n", mean(dist), sim.bonds.l0(1));
	printf("Av. angle %.3f, eq. angle set to %.3f \n", mean(_angles), sim.angles.a0(1));
	printf("Av. dihedral %.3f\n", mean(dihedrals));

end
