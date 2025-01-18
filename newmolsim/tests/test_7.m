

function [dist _angles dihedrals] = test_7()

	niter = 1e4;
	dt = 1e-3;

	addpath("../mfiles/"); addpath("../mex/");

	molconfgen("butane.xyz", "butane.top", 500, 0.1);

	sim = molsim();
	sim.setconf("conf.xyz");

	sim.integrator.dt = dt;
	sim.thermostat.temperature = 3.0;

	sim.setbonds("bonds.top");
	nbonds = sim.bonds.nbonds;
	sim.bonds.springs = 33e3*ones(nbonds,1); sim.bonds.l0 = 0.4*ones(nbonds,1); 

	sim.setangles("angles.top");
	nangles = sim.angles.nangles;	
	sim.angles.springs = 866*ones(nangles,1); sim.angles.a0 = 1.9*ones(nangles, 1);

	sim.setdihedrals("dihedrals.top");
	ndihedrals = sim.dihedrals.ndihedrals;
	for n=1:ndihedrals
		sim.dihedrals.coeffs(n,:) = [15.5000,  20.3050, -21.9170, -5.1150,  43.8340, -52.6070]; 
	end
	
	printf("No. bonds %d, no. angles %d, no. dihedrals %d\n", nbonds, nangles, ndihedrals);

	for n=1:4:sim.natoms
		sim.atoms.exclude(n, 1:3) = [n+1, n+2, n+3]; 
		sim.atoms.exclude(n+1, 1:3) = [n, n+2, n+3];
		sim.atoms.exclude(n+2, 1:3) = [n, n+1, n+3];
		sim.atoms.exclude(n+3, 1:3) = [n, n+1, n+2];
	end	

	ekin = zeros(1, niter); epot = zeros(1, niter); 
	dist = zeros(1, niter); _angles = zeros(1, niter); dihedrals = zeros(1, niter);

	for n=1:niter
		epot(n) = sim.pairforce.lj(sim.atoms, "CC", [2.5, 1.0, 1.0, 1.0]);   
		epot(n) = epot(n) + sim.bonds.harmonic(sim.atoms, 0);
		epot(n) = epot(n) + sim.angles.harmonic(sim.atoms, 0);
		epot(n) = epot(n) + sim.dihedrals.ryckbell(sim.atoms, 0);

		ekin(n) = sim.integrator.step(sim.atoms, sim.pairforce);
		
		if rem(n, 10)==0
			sim.scalebox(1.46, 0.999);
		end

		dist(n) = sim.atoms.getdist(1,2);
		_angles(n) = sim.atoms.getangle(1,2,3);
		dihedrals(n) = sim.atoms.getdihedral(1,2,3,4);
	end
	
	ekin = ekin./sim.natoms; epot = epot./sim.natoms; etot = epot + ekin;

	sim.atoms.save("final.xyz");
end
