clear all;

function [epot, ekin] = mdsim(sim, ewald=false, nloops=1, cutoff=3.0)

	epot = zeros(nloops, 1); ekin = zeros(nloops, 1);
	for n=1:nloops
		epot(n) = sim.lennardjones("AA", [2.0^(1/6), 1.0, 1.0, 1.0]);   
		epot(n) += sim.lennardjones("BB", [2.0^(1/6), 1.0, 1.0, 1.0]);   
		epot(n) += sim.lennardjones("AB", [2.0^(1/6), 1.0, 1.0, 1.0]);   

		if ewald
			epot(n) += sim.ewald([1, cutoff, 1, 10]);
		else	
			epot(n) += sim.sfcoulomb(cutoff);
		end
			
		sim.applythermostat();
		ekin(n) = sim.leapfrog();
	end

end

# Remove if script is run after package installation 
addpath("../inst/"); addpath("../src/");

# State point etc
nx = 10; npart = nx^3;
dens = 0.368; temp = 0.0177; 
lbox = (npart/dens)^(1/3);
nloops = 1e3;
cutoff = 3.0;

# Simulation setup
sim = molsim();

sim.setconf([nx,nx,nx], [lbox, lbox, lbox], temp);

sim.atoms.t(1:2:end) = 'A'; 
sim.atoms.t(2:2:end) = 'B';

sim.atoms.q = (-1).^[1:npart];

sim.setthermostat("relax", 0.2, 0.1);
sim.pairforce.max_cutoff = cutoff;

[epot, ekin] = mdsim(sim, true, nloops, cutoff);

epot = epot./sim.natoms; ekin = ekin./sim.natoms;

# Plot the energies
plot([1:nloops], epot, '-s', [1:nloops], ekin, '-o');

