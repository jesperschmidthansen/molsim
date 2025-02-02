
function [ekin, epot] = test_7()

	addpath("../mfiles/"); addpath("../mex/");

	niter = 1e4; dt = 0.01; temperature = 1.0;
	
	sim = molsim();
	sim.setconf([10,10,10], [7, 7, 7], temperature);
	sim.pairforce.max_cutoff = 1.0;
	sim.pairforce.skin = 1.0;
	sim.integrator.dt = dt;
	
	# [cf, aij, sigma, dt];
	dpd_params = [1.0, 25.0, 3.0, dt]; 
	dpd_lambda = 0.5;

	ekin = zeros(1, niter); epot = zeros(1, niter);

	for n=1:niter
		epot(n) = sim.pairforce.dpd(sim.atoms, "AA", dpd_params, temperature);   
		ekin(n) = sim.integrator.dpd(sim.atoms, sim.pairforce, dpd_lambda);
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
