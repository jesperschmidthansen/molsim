
function test_1()
	 
	addpath("../mfiles/"); addpath("../mex/");

	T0 = 1.42;
	niter = 2e3;

	p = atoms("start.xyz"); 
	intgr = integrator(); 
	prfrc = prforce();
	thmstat = thermostat(p, T0);

	p.m(1:2:end)=2.0;
	p.setvelocities(T0);

	ekin = zeros(niter,1); momentum = zeros(niter, 3); 
	for n=1:niter
		prfrc.lj(p, "AA", [2.5, 1.0, 1.0, 1.0]);   
		thmstat.relaxtemp(p);
		ekin(n) = intgr.step(p, prfrc);
		momentum(n, :) = p.momentum();
	end

	T = 2/3*mean(ekin(end-100:end))./p.natoms; 

	index = [1:n];
	plot(index, ekin./p.natoms, ";ekin;", index, momentum, ";momentum;");
	print("test_1.pdf", '-dpdf');

	printf("test_1 output:\n");
	printf("Ekin: %.3f +/- %.3f  Norm. mean temperature %f\n", mean(ekin)./p.natoms, std(ekin)./p.natoms, T/T0);

end
