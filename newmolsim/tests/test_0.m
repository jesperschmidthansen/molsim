
function test_0()

	addpath("../mfiles/"); addpath("../mex/");

	niter = 1e4;

	p = atoms("start.xyz"); 
	intgr = integrator(); 
	prfrc = prforce(); prfrc.skin = 0.25;

	ekin = zeros(1, niter); epot = zeros(1, niter);
	p.m(1:2:end)=1.32;

	tic();
	for n=1:niter

		epot(n) = prfrc.lj(p, "AA", [2.5, 1.0, 1.0, 1.0]);   
		ekin(n) = intgr.step(p, prfrc);

	end
	t = toc(); 

	sps = n/t;
	ekin = ekin./p.natoms; epot = epot./p.natoms; etot = epot + ekin;
	spnb = niter/prfrc.neighb_updates;
	 
	index = [1:n];
	plot(index, ekin, ";ekin;", index, epot, ";epot;", index, epot+ekin, ";etot;")
	print("test_0.pdf", '-dpdf');

	printf("test_0 output:\n");
	printf("Etot: %1.3e +/- %1.3e   ", mean(etot), std(etot));
	printf("Steps per seconds: %.0f   Steps per build  %.3f\n", sps, spnb);

endfunction
