function test_8()

	addpath("../src/"); addpath("../inst/");

	randn("seed", 42);
	niter = 1e3;
	
	sim = molsim();
	sim.setconf([10,10,10], [10.557, 10.557, 10.557], 1.0);

	ekin = zeros(1, niter); epot = zeros(1, niter);

	sim.atoms.m(1:2:end)=1.32; sim.atoms.resetmom();

	tic();
	for n=1:niter
		epot(n) = sim.pairforce.lj(sim.atoms, "AA", [2.5, 1.0, 1.0, 1.0]);   
		ekin(n) = sim.integrator.lf(sim.atoms, sim.pairforce);
	end
	t_1 = niter/toc();
	etot_1 = mean(epot+ekin);
 
	index = [1:n];
	plot(index, ekin, ";ekin;", index, epot, ";epot;", index, epot+ekin, ";etot;")
	hold on;

	sim.atoms.save("tmp.mat"); clear sim;
	sim = molsim();	sim.setconf("tmp.mat");

	tic();
	for n=1:niter
		epot(n) = sim.pairforce.lj(sim.atoms, "AA", [2.5, 1.0, 1.0, 1.0]);   
		ekin(n) = sim.integrator.lf(sim.atoms, sim.pairforce);
	end
	t_2 =  niter/toc();
	etot_2 = mean(ekin+epot);

	index=[n+1:2*n];	
	plot(index, ekin, ";ekin;", index, epot, ";epot;", index, epot+ekin, ";etot;")
	hold off;

	printf("test_8 output:\n");
	printf("First run: Energy %e. Steps per second %.f\n", etot_1, t_1); 
	printf("Second run: Energy %e. Steps per second %.f\n", etot_2, t_2); 

	print("test_8.pdf", '-dpdf');

end

