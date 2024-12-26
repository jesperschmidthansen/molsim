
function dist = test_5()

	addpath("../mfiles/"); addpath("../mex/");

	niter = 1e4;

	p = atoms("mol.xyz"); 
	intgr = integrator(); 
	prfrc = prforce(); 

	nbonds = p.natoms/2; # Remainder should be checked
	b = bonds(nbonds);
	for n=1:nbonds
		b.pidx(n, :) = [2*n-1, 2*n]; 
	end
	b.springs = 500*ones(nbonds,1); b.l0 = 1.3*ones(nbonds,1); 

	ekin = zeros(1, niter); epot = zeros(1, niter); dist = zeros(1, niter);
	
	p.resetmom();

	for n=1:niter
		epot(n) = prfrc.lj(p, "AA", [2.5, 1.0, 1.0, 1.0]);   
		epot(n) = epot(n) + b.harmonic(p, 0);
		ekin(n) = intgr.step(p, prfrc);
		
		dist(n) = p.getdist(1,2);
	end
	
	ekin = ekin./p.natoms; epot = epot./p.natoms; etot = epot + ekin;
	spnb = niter/prfrc.neighb_updates; 
	mom = p.getmom();

	printf("test_5 output:\n");
	printf("Etot: %1.3e +/- %1.3e   ", mean(etot), std(etot));
	printf("Momentum: %e %e %e\n", mom(1), mom(2), mom(3));
	printf("Steps per build  %.3f\n", spnb);
	printf("Av. bond length %.3f, eq. bond length set to %.3f \n", mean(dist), b.l0(1));

end
