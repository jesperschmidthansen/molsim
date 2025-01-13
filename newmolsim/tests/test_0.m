
function [ekin, epot] = test_0()

	addpath("../mfiles/"); addpath("../mex/");

	niter = 1e4;
	
	sim = molsim([10,10,10], [10.557, 10.557, 10.557], 1.0);

	ekin = zeros(1, niter); epot = zeros(1, niter);

	sim.atoms.m(1:2:end)=1.32; sim.atoms.resetmom();

	for n=1:niter
		epot(n) = sim.pairforce.lj(sim.atoms, "AA", [2.5, 1.0, 1.0, 1.0]);   
		ekin(n) = sim.integrator.step(sim.atoms, sim.pairforce);
	end
	
	ekin = ekin./sim.natoms; epot = epot./sim.natoms; etot = epot + ekin;
	spnb = niter/sim.pairforce.neighb_updates; 
	mom = sim.atoms.getmom();
	 
	index = [1:n];
	plot(index, ekin, ";ekin;", index, epot, ";epot;", index, epot+ekin, ";etot;")
	print("test_0.pdf", '-dpdf');

	printf("test_0 output:\n");
	printf("Etot: %1.3e +/- %1.3e   ", mean(etot), std(etot));
	printf("Momentum: %e %e %e\n", mom(1), mom(2), mom(3));
	printf("Steps per build  %.1f\n", spnb);

end
