
function test_3()
	 
	addpath("../inst/"); addpath("../src/");

	T0 = 1.42;
	niter = 1e4;

	sim = molsim();
	sim.setconf([10,10,10], [11, 11, 11], T0);

	sim.atoms.t(1:500)='W';	sim.atoms.t(501:end) = 'F';
	sim.setthermostat('relax', 'W', T0, 0.01);
	
	ekin = zeros(niter,1); 
	for n=1:niter
		sim.lennardjones("FF", [2.5, 1.0, 1.0, 1.0]); 
		sim.lennardjones("WW", [2.5, 1.0, 1.0, 1.0]); 
  		sim.lennardjones("FW", [2.5, 1.0, 1.0, 1.0]); 
		
		sim.atoms.tether('W', 300);
	
		sim.applythermostat();
		ekin(n) = sim.leapfrog();
	end

	sim.atoms.save("tether.xyz");

	plot3(sim.atoms.r(1:500,1), sim.atoms.r(1:500,2), sim.atoms.r(1:500,3), 'o', 'markerfacecolor', 'red');
	hold on
	plot3(sim.atoms.r(501:end,1), sim.atoms.r(501:end,2), sim.atoms.r(501:end,3), 'o', 'markerfacecolor', 'blue');
	hold off
	view([0,0,0]);
	print("test_3.pdf", '-dpdf');

	printf("Temperature %.2f\n", 2/3*mean(ekin(end-100:end))/sim.natoms);
end
