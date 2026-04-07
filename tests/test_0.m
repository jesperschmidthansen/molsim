
function [ekin, epot] = test_0()

	addpath("../inst/"); addpath("../src/");

	niter = 1e4;
	
	sim = molsim();
	sim.setconf([10,10,10], [10.557, 10.557, 10.557], 1.0);

	sim.atoms.m(1:2:end)=1.32; sim.atoms.resetmom();

	ekin = zeros(1, niter); epot = zeros(1, niter); P = zeros(3,3);
	for n=1:niter
		[epot(n) Pconf] = sim.pairforce.lj(sim.atoms, "AA", [2.5, 1.0, 1.0, 1.0]);   
		[ekin(n) Pkin] = sim.integrator.lf(sim.atoms, sim.pairforce);

		P = P + Pconf + Pkin;		
	end
	
	P = P./niter;
	
	ekin = ekin./sim.natoms; epot = epot./sim.natoms; etot = epot + ekin;
	spnb = niter/sim.pairforce.neighb_updates; 
	mom = sim.atoms.getmom();
	 
	index = [1:n];
	plot(index, ekin, ";ekin;", index, epot, ";epot;", index, epot+ekin, ";etot;")
	print("test_0.pdf", '-dpdf');

	printf("test_0 output:\n");
	printf("Etot: %1.3e +/- %1.3e  Pressure: %1.2f  \n", mean(etot), std(etot), (P(1,1)+P(2,2)+P(3,3))./3);
	printf("Momentum: %e %e %e\n", mom(1), mom(2), mom(3));
	printf("Steps per build  %.1f\n", spnb);

end
