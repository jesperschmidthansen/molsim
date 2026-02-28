clear all;

addpath("../inst/"); addpath("../src/");

niter = 1e4;

sim = molsim();
sim.setconf([10,10,10], [10.557, 10.557, 10.557], 1.0);

ekin = zeros(niter,1); epot = zeros(niter,1);
for n=1:niter
	epot(n) = sim.lennardjones("AA", [2.5, 1.0, 1.0, 1.0]);   
	ekin(n) = sim.leapfrog();
end

plot([1:niter], ekin/sim.natoms, [1:niter], epot/sim.natoms, [1:niter], (epot+ekin)./sim.natoms);
pause(10);

