
function test_3()
	 
	addpath("../inst/"); addpath("../src/");

	T0 = 1.42;
	niter = 1e4;

	sim = molsim();
	sim.setconf([10,10,10], [11, 11, 11], 2.0);

	sim.atoms.t(1:500)='W';	sim.atoms.t(501:end) = 'F';
	sim.thermostat.temperature = T0;
	sim.thermostat.settype(sim.atoms, 'W');
	
	ekin = zeros(niter,1); 
	for n=1:niter
		sim.pairforce.lj(sim.atoms, "FF", [2.5, 1.0, 1.0, 1.0]); 
		sim.pairforce.lj(sim.atoms, "WW", [2.5, 1.0, 1.0, 1.0]); 
  		sim.pairforce.lj(sim.atoms, "FW", [2.5, 1.0, 1.0, 1.0]); 
		
		sim.atoms.tether('W', 300);
	
		sim.thermostat.nosehoover(sim.atoms, sim.integrator.dt);
		ekin(n) = sim.integrator.lf(sim.atoms, sim.pairforce);
	end

	sim.atoms.save("tether.xyz");

	plot3(sim.atoms.r(1:500,1), sim.atoms.r(1:500,2), sim.atoms.r(1:500,3), 'o', 'markerfacecolor', 'red');
	hold on
	plot3(sim.atoms.r(501:end,1), sim.atoms.r(501:end,2), sim.atoms.r(501:end,3), 'o', 'markerfacecolor', 'blue');
	hold off
	view([0,0,0]);
	print("test_3.pdf", '-dpdf');

	printf("Ekin %.2f\n", mean(ekin)/sim.natoms);
end
