clear all;

addpath("../inst/"); addpath("../src/");

niter = 1e4;

sim = molsim();
sim.setconf([10,10,10], [10.557, 10.557, 10.557], 1.0);
sim.thermostat.temperature = 1.5;

for n=1:niter
	sim.lennardjones("AA", [2.5, 1.0, 1.0, 1.0]);   
	
	sim.nosehoover();
	sim.leapfrog();
end

