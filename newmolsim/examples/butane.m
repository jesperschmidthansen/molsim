

clear all;

addpath("../inst/"); addpath("../src/");

niter = 1e4; dt = 1e-3;
temp0 = 3.0; dens0 = 1.49;

molconfgen("../resources/molconf/butane.xyz", "../resources/molconf/butane.top", 500, 0.1);

sim = molsim();
sim.setconf("conf.xyz");

sim.integrator.dt = dt;
sim.thermostat.temperature = temp0;

sim.setbonds("bonds.top");
nbonds = sim.bonds.nbonds;
sim.bonds.springs = 33e3*ones(nbonds,1); sim.bonds.l0 = 0.4*ones(nbonds,1); 

sim.setangles("angles.top");
nangles = sim.angles.nangles;	
sim.angles.springs = 866*ones(nangles,1); sim.angles.a0 = 1.9*ones(nangles, 1);

sim.setdihedrals("dihedrals.top");
ndihedrals = sim.dihedrals.ndihedrals;
for n=1:ndihedrals
	sim.dihedrals.coeffs(n,:) = [15.5000,  20.3050, -21.9170, -5.1150,  43.8340, -52.6070]; 
end

sim.atoms.setexclusions(sim.dihedrals.pidx, "dihedrals");


for n=1:niter
	sim.lennardjones("CC", [2.5, 1.0, 1.0, 1.0]);   
	
	sim.harmonicbond(0);
	sim.harmonicangle(0);
	sim.ryckbell(0);

	sim.nosehoover();
	sim.leapfrog();
	
	if rem(n, 10)==0
		sim.scalebox(dens0);
	end

end

sim.atoms.save("final.xyz");
