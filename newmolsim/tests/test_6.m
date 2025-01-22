
function [ekin, epot] = test_6()

	addpath("../mfiles/"); addpath("../mex/");

	niter = 1e3; cutoff = 3.5;
	
	sim = molsim();
	sim.setconf([14,14,14], [15.0, 15.0, 15.0], 1.0);

	sim.atoms.q = (-1).^[1:sim.natoms];
	sim.pairforce.max_cutoff = cutoff;

	ekin = zeros(1, niter); epot = zeros(1, niter);
	for n=1:niter
		epot(n) = sim.pairforce.lj(sim.atoms, "AA", [2.5, 1.0, 1.0, 1.0]);   
		epot(n) = epot(n) + sim.pairforce.sf(sim.atoms, cutoff);
		ekin(n) = sim.integrator.step(sim.atoms, sim.pairforce);
	end
	
	ekin = ekin./sim.natoms; epot = epot./sim.natoms; etot = epot + ekin;
	mom = sim.atoms.getmom();
	 
	index = [1:n];
	plot(index, ekin, ";ekin;", index, epot, ";epot;", index, epot+ekin, ";etot;")
	print("test_6.pdf", '-dpdf');

	printf("test_6 output:\n");
	printf("Total charge %e\n", sum(sim.atoms.q));
	printf("Etot: %1.3e +/- %1.3e   ", mean(etot), std(etot));
	printf("Momentum: %e %e %e\n", mom(1), mom(2), mom(3));

end
