clear all;

addpath("../inst/"); addpath("../src/"); addpath("../resources/setup/");

niter = 1e4; dt = 5e-4;
dens0 = 3.16; temp0 = 298.15/78.2; 

cutoff = 2.9;

lbond = 0.316; kspring = 68421; 
angle = 1.97; kangle = 490;

molconf("../resources/molconf/water.xyz", "../resources/molconf/water.top", [15, 15, 15]);

sim = molsim();
sim.setconf("configuration.xyz");
sim.settop("topology.mat");

sim.atoms.setvels(temp0);
sim.pairforce.max_cutoff = cutoff;
sim.integrator.dt = dt;
sim.setthermostat("nh", temp0, 10.0);

nbonds = sim.bonds.nbonds;
sim.bonds.springs = kspring*ones(nbonds,1); 
sim.bonds.l0 = lbond*ones(nbonds,1); 

nangles = sim.angles.nangles;	
sim.angles.springs = kangle*ones(nangles,1); 
sim.angles.a0 = angle*ones(nangles, 1);

sim.atoms.setexclusions(sim.angles.pidx, "angles");

for n=1:niter
	sim.lennardjones("OO", [2.5, 1.0, 1.0, 1.0]);   
	sim.sfcoulomb(cutoff);
	sim.harmonicbond();
	sim.harmonicangle();

	sim.applythermostat();
	ekin = sim.leapfrog();
	
	if rem(n, 50)==0
		if sim.natoms/sim.volume < dens0 
			sim.scalebox(dens0, [1:3], 0.999);
		end

		printf("\r It. no. %d, dens. %f, temp. %f ", ... 
			n, sim.natoms/sim.volume, ekin*2/(3*sim.natoms)); 
		fflush(stdout);
	end

end
printf("\n");

sim.atoms.save("water-sys.xyz"); sim.atoms.save("water-sys.mat");

