clear all;

addpath("../mfiles/"); addpath("../mex/");

niter = 1e4; dt = 5e-4;
temp0 = 3.0;dens0 = 3.16; temp0 = 298.15/78.2; 

cutoff = 2.9;

lbond = 0.316; kspring = 68421; 
angle = 1.97; kangle = 490;

#system("cp ../resources/molconf/water.xyz ./");
#system("cp ../resources/molconf/water.top ./");i

molconfgen("../resources/molconf/water.xyz", "../resources/molconf/water.top", 2000, 0.1);

sim = molsim();
sim.setconf("conf.xyz");

sim.pairforce.max_cutoff = cutoff;
sim.integrator.dt = dt;
sim.thermostat.temperature = temp0;

sim.setbonds("bonds.top");
nbonds = sim.bonds.nbonds;
sim.bonds.springs = kspring*ones(nbonds,1); sim.bonds.l0 = lbond*ones(nbonds,1); 

sim.setangles("angles.top");
nangles = sim.angles.nangles;	
sim.angles.springs = kangle*ones(nangles,1); sim.angles.a0 = angle*ones(nangles, 1);

sim.atoms.setexclusions(sim.angles.pidx, "angles");


for n=1:niter
	sim.pairforce.lj(sim.atoms, "OO", [2.5, 1.0, 1.0, 1.0]);   
	sim.pairforce.sf(sim.atoms, cutoff);
	sim.bonds.harmonic(sim.atoms, 0);
	sim.angles.harmonic(sim.atoms, 0);

	sim.thermostat.nosehoover(sim.atoms);
	sim.integrator.lf(sim.atoms, sim.pairforce);
	
	if rem(n, 10)==0
		sim.scalebox(dens0);
	end

end

sim.atoms.save("final.xyz");
