

clear all;

addpath("../inst/"); addpath("../src/");  addpath("../resources/setup/");


niter = 1e4; dt = 1e-3;
temp0 = 3.0; dens0 = 1.49;

system("rm *top");
molconf("../resources/molconf/butane.xyz", "../resources/molconf/butane.top", [10 10 10], 4.0);

sim = molsim();
sim.setconf("start.xyz");

sim.integrator.dt = dt;
sim.thermostat.temperature= temp0;

sim.settop();

nbonds = sim.bonds.nbonds;
sim.bonds.springs = 33e3*ones(nbonds,1); sim.bonds.l0 = 0.4*ones(nbonds,1); 

nangles = sim.angles.nangles;	
sim.angles.springs = 866*ones(nangles,1); sim.angles.a0 = 1.9*ones(nangles, 1);

ndihedrals = sim.dihedrals.ndihedrals;
for n=1:ndihedrals
	sim.dihedrals.coeffs(n,:) = [15.5000,  20.3050, -21.9170, -5.1150,  43.8340, -52.6070]; 
end

sim.atoms.setexclusions(sim.dihedrals.pidx, "dihedrals");

ekin = zeros(niter,1); epot = zeros(niter,1);

for n=1:niter
	epot(n) = sim.lennardjones("CC", [2.5, 1.0, 1.0, 1.0]);   
	
	epot(n) += sim.harmonicbond(0);
	epot(n) += sim.harmonicangle(0);
	epot(n) += sim.ryckbell(0);

	sim.nosehoover();
	ekin(n) = sim.leapfrog();
	
	if rem(n, 10)==0
		sim.scalebox(dens0);
	end

end
sim.natoms/sim.volume

sim.atoms.save("final.xyz");

plot([1:niter], ekin/sim.natoms, [1:niter], epot/sim.natoms, [1:niter], (epot+ekin)./sim.natoms);
pause(10);

