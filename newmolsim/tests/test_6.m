
function [ekin, epot] = test_6()

	addpath("../mfiles/"); addpath("../mex/");

	niter = 1e3; cutoff = 3.5;
	
	p = atoms([14,14,14], [15.0, 15.0, 15.0], 1.0);
	intgr = integrator(); 
	prfrc = prforce(cutoff); 

	p.q = (-1).^[1:p.natoms];

	ekin = zeros(1, niter); epot = zeros(1, niter);
	for n=1:niter
		epot(n) = prfrc.lj(p, "AA", [2.5, 1.0, 1.0, 1.0]);   
		epot(n) = epot(n) + prfrc.sf(p, cutoff);
		ekin(n) = intgr.step(p, prfrc);
	end
	
	ekin = ekin./p.natoms; epot = epot./p.natoms; etot = epot + ekin;
	mom = p.getmom();
	 
	index = [1:n];
	plot(index, ekin, ";ekin;", index, epot, ";epot;", index, epot+ekin, ";etot;")
	print("test_6.pdf", '-dpdf');

	printf("test_6 output:\n");
	printf("Total charge %e\n", sum(p.q));
	printf("Etot: %1.3e +/- %1.3e   ", mean(etot), std(etot));
	printf("Momentum: %e %e %e\n", mom(1), mom(2), mom(3));

end
