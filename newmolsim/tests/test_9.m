function [epot, ekin] = test_9()

	niter = 1e4;
	sim = molsim();
	sim.setconf([10,10,10], [15,15,15], 1.0);
	sim.thermostat.temperature = 1.0;
	sim.thermostat.tauQ = 50.0;

	for n=1:100
		sim.addatom(rand.*sim.lbox, 'B', 1.0, 0.);
	end

	ekin = zeros(1, niter); epot = zeros(1, niter);

	epsilon = 0.0;
	for n=1:niter
		epot(n) = sim.lennardjones("AA", [2.5, 1.0, 1.0, 1.0]);   
		epot(n) += sim.lennardjones("AB", [2.5, 1.0, epsilon, 1.0]);   
		epot(n) += sim.lennardjones("BB", [2.5, 1.0, epsilon, 1.0]);   

		sim.nosehoover();
		ekin(n) = sim.leapfrog();

		if epsilon < 1.0
			epsilon += 1e-3;
		end
	end

	ekin = ekin./sim.natoms; epot = epot./sim.natoms; 
	mom = sim.atoms.getmom();
	 
	index = [1:n];
	plot(index, ekin, ";ekin;", index, epot, ";epot;")
	print("test_9.pdf", '-dpdf');

	printf("test_9 output:\n");
	printf("Momentum: %e %e %e\n", mom(1), mom(2), mom(3));

end

