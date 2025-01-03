
function test_2()

	addpath("../mfiles/"); addpath("../mex/");

	niter = 1e3;
	ekin = zeros(niter,1); epot = zeros(niter,1);	

	for nthreads=1:8 
		p = atoms("start.xyz"); 
		intgr = integrator(); 
		prfrc = prforce(); 
		 
		ms_setomp(nthreads);

		tic();
		for n=1:niter
			prfrc.lj(p, "AA", [2.5, 1.0, 1.0, 1.0]);   
			intgr.step(p, prfrc);
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
