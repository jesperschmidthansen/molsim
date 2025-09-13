
function test_2()

	addpath("../inst/"); addpath("../src/");

	niter = 1e3;
	ekin = zeros(niter,1); epot = zeros(niter,1);	

	for nthreads=1:8 
		sim = molsim();
		sim.setconf([10,10,10], [11, 11, 11], 2.0);
		 
		sim.setnthreads(nthreads);

		tic();
		for n=1:niter
			sim.pairforce.lj(sim.atoms, "AA", [2.5, 1.0, 1.0, 1.0]);   
			sim.integrator.lf(sim.atoms, sim.pairforce);
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
