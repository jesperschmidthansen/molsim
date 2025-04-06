clear all;

addpath("../mfiles/"); addpath("../mex/");

niter = 1e4;

sim = molsim();
sim.setconf([10,10,10], [10.557, 10.557, 10.557], 1.0);
sim.thermostat.temperature = 1.5;

for n=1:niter
	sim.pairforce.lj(sim.atoms, "AA", [2.5, 1.0, 1.0, 1.0]);   
	
	sim.thermostat.nosehoover(sim.atoms);
	sim.integrator.lf(sim.atoms, sim.pairforce);
end

