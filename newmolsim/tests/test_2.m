
function test_2()

	addpath("../mfiles/"); addpath("../mex/");

	niter = 1e3;
	ekin = zeros(niter,1); epot = zeros(niter,1);	

	for nthreads=1:8 
		sim = molsim([10,10,10], [11, 11, 11], 2.0);
		 
		sim.set_nthreads(nthreads);

		tic();
		for n=1:niter
			sim.pairforce.lj(sim.atoms, "AA", [2.5, 1.0, 1.0, 1.0]);   
			sim.integrator.step(sim.atoms, sim.pairforce);
		end
		t = toc(); 

		sps(nthreads) = n/t; 
	end

	printf("test_2 output:\n");
	for n=1:8
		printf("Num threads: %d -> Steps per second: %.0f  \n", n, sps(n));
	end

	plot(1:8, sps, '-o;sps;');	print("test_2.pdf", '-dpdf');

end
