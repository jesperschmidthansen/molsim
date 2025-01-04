
function [ekin, epot] = test_0()

	addpath("../mfiles/"); addpath("../mex/");

	niter = 1e4;
	
	p = atoms([10,10,10], [10.557, 10.557, 10.557], 1.0);
	intgr = integrator(); 
	prfrc = prforce(); 

	ekin = zeros(1, niter); epot = zeros(1, niter);

	p.m(1:2:end)=1.32; p.resetmom();

	for n=1:niter
		epot(n) = prfrc.lj(p, "AA", [2.5, 1.0, 1.0, 1.0]);   
		ekin(n) = intgr.step(p, prfrc);
	end
	
	ekin = ekin./p.natoms; epot = epot./p.natoms; etot = epot + ekin;
	spnb = niter/prfrc.neighb_updates; 
	mom = p.getmom();
	 
	index = [1:n];
	plot(index, ekin, ";ekin;", index, epot, ";epot;", index, epot+ekin, ";etot;")
	print("test_0.pdf", '-dpdf');

	printf("test_0 output:\n");
	printf("Etot: %1.3e +/- %1.3e   ", mean(etot), std(etot));
	printf("Momentum: %e %e %e\n", mom(1), mom(2), mom(3));
	printf("Steps per build  %.1f\n", spnb);

end
