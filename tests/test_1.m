
function test_1()
	 
	addpath("../inst/"); addpath("../src/");

	T0 = 1.1;
	niter = 1e4;
	
	sim = molsim();
	sim.setconf([10,10,10], [11, 11, 11], T0);

	sim.atoms.resetmom();

	sim.atoms.t(501:end)='B';
	sim.setthermostat("nh", T0, 10);
	
	ekin = zeros(niter,1);
	for n=1:niter
		sim.lennardjones("AA", [2.5, 1.0, 1.0, 1.0]);   
		sim.lennardjones("AB", [2.5, 1.0, 1.0, 1.0]);   
		sim.lennardjones("BB", [2.5, 1.0, 1.0, 1.0]);   

		sim.applythermostat();		
		ekin(n) = sim.leapfrog();
	end

	T = 2/3*mean(ekin(end-100:end))./sim.natoms; 
	mom = sim.atoms.getmom();

	plot(2/3*ekin./sim.natoms);
	print("test_1.pdf", '-dpdf');

	printf("test_1 output:\n");
	printf("Ekin: %.3f +/- %.3f  Norm. mean temperature %f\n", mean(ekin)./sim.natoms, std(ekin)./sim.natoms, T/T0);
	printf("Momentum: %1.3e %1.3e %1.3e \n", mom(1), mom(2), mom(3));
	
end
