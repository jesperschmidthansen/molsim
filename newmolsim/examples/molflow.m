

clear all;

addpath("../inst/"); addpath("../src/");

niter = 1e4; dt = 1e-3;
temp0 = 3.0; 
dztot = 0.0; dz = 1e-3;

molconf("diatomic.xyz", "diatomic.top", [5 5 40], 2.0);
molslitpore("start.xyz", [3.0 3.0], [1.0, 1.0]);

sim = molsim();
sim.setconf("molslitpore.xyz");

sim.integrator.dt = dt;

sim.setbonds("bonds.top");
nbonds = sim.bonds.nbonds;
sim.bonds.springs = 200*ones(nbonds,1); sim.bonds.l0 = 1.0*ones(nbonds,1); 

sim.atoms.setexclusions(sim.bonds.pidx, "bonds");

for n=1:niter
	epot(n) = sim.lennardjones("CC", [2^(1/6), 1.0, 1.0, 1.0]);
	epot(n) = sim.lennardjones("Cw", [2^(1/6), 1.0, 1.0, 1.0]);
	epot(n) = sim.lennardjones("CW", [2^(1/6), 1.0, 1.0, 1.0]);

	epot(n) += sim.harmonicbond(0);

	sim.atoms.tether('W', 300);
	sim.atoms.tether('w', 300);

	sim.thermostat.relaxttemp(sim.atoms, 'w', 1.0);
	sim.thermostat.relaxttemp(sim.atoms, 'W', 1.0);

	ekin(n) = sim.leapfrog();

	if dztot < 10.0	
		sim.atoms.mvlattice('w', [0, 0, -dz]);
		dztot += dz;
	end

end

sim.atoms.save("final.xyz");


