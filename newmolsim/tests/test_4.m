
function test_4()

	addpath("../mfiles/"); addpath("../mex/");

	niter = 1e3;
	ndims = [10, 13, 15, 20, 25];
	lbox = (ndims.^3./0.8).^(1/3);

	for m=1:length(ndims)
		sim = molsim();
		sim.setconf(ndims(m).*[1, 1, 1], lbox(m)*[1,1,1], 2.0); 

		tic();
		for n=1:niter
			sim.pairforce.lj(sim.atoms, "AA", [2.5, 1.0, 1.0, 1.0]);   
			sim.integrator.lf(sim.atoms, sim.pairforce);
		end
		t = toc(); 

		sps(m) = n/t; 
	end

	natoms = ndims.^3;
	printf("test_4 output:\n");
	for n=1:length(natoms)
		printf("N. atoms: %d -> Steps per second: %.0f  \n", natoms(n), sps(n));
	end

	loglog(natoms, sps, '-o', 'linewidth', 2);
	print("test_4.pdf", '-dpdf');

end
