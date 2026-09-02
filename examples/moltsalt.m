clear all;

function [epot, ekin] = mdstep(sim, ewald=false, cutoff=3.0, nloops = 1)

	for n=1:nloops
		epot = sim.lennardjones("AA", [2.0^(1/6), 1.0, 1.0, 1.0]);   
		epot += sim.lennardjones("BB", [2.0^(1/6), 1.0, 1.0, 1.0]);   
		epot += sim.lennardjones("AB", [2.0^(1/6), 1.0, 1.0, 1.0]);   

		if ewald
			epot += sim.ewald([1, cutoff, 1]);
		else	
			epot += sim.sfcoulomb(3.0);
		end

		sim.applythermostat();

		ekin = sim.leapfrog();

		if rem(n,100)==0
			printf("\r Did %d ", n); fflush(stdout);
		end
	end

	if nloops > 100 
		printf("\n");
	end

end

# Remove if script is run after package installation 
addpath("../inst/"); addpath("../src/");

# State point etc
nx = 10; npart = nx^3;
dens = 0.368; temp = 0.0177; 
lbox = (npart/dens)^(1/3);
nloops = 1e2;
cutoff = 3.0;

# Simulation setup
sim = molsim();

sim.setconf([nx,nx,nx], [lbox, lbox, lbox], temp);

sim.atoms.t(1:2:end) = 'A'; 
sim.atoms.t(2:2:end) = 'B';

sim.atoms.q = (-1).^[1:npart];

sim.setthermostat("relax", 0.2, 0.1);
sim.pairforce.max_cutoff = cutoff;

mdstep(sim, false, cutoff, 5e3);

printf("Did eq.\n"); fflush(stdout);

sim.setautosave(100);
sim.thermostat.temperature = temp;
# Main MD loop
ekin = zeros(nloops,1); epot = zeros(nloops,1);
for n=1:nloops
	[epot(n) ekin(n)] = mdstep(sim, true);
end


# Plot the energies
plot([1:nloops], 2/3*ekin/sim.natoms, [1:nloops], epot/sim.natoms, [1:nloops], (epot+ekin)./sim.natoms);

