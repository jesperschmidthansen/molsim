
# Example of calibration of a molecular system - here the flexible water model with Coulomb SF
# See Wu et al., J.Chem.Phys., 124: 024503 (2006)

clear all;

addpath("../inst/"); addpath("../src/"); addpath("../resources/setup/");

niter = 100e3; dt = 5e-4;
dens0 = 3.16; temp0 = 298.15/78.2; 

cutoff = 2.9;

lbond = 0.316; kspring = 68421; 
angle = 1.97; kangle = 490;

# If system too small after compression the optimization algorithm issues a warning 
# - and the simulation will likely crash
molconf("../resources/molconf/water.xyz", "../resources/molconf/water.top", [15, 15, 15], 2.0);

# Load the generated system files
sim = molsim();
sim.setconf("configuration.xyz");
sim.settop("topology.mat");

# Misc. settings
sim.atoms.setvels(temp0);
sim.pairforce.max_cutoff = cutoff;
sim.integrator.dt = dt;
sim.setthermostat("nh", temp0, 10.0);

# Intra-molecular details
nbonds = sim.bonds.nbonds;
sim.bonds.springs = kspring*ones(nbonds,1); 
sim.bonds.l0 = lbond*ones(nbonds,1); 

nangles = sim.angles.nangles;	
sim.angles.springs = kangle*ones(nangles,1); 
sim.angles.a0 = angle*ones(nangles, 1);

sim.atoms.setexclusions(sim.angles.pidx, "angles");

# Main MD loop
for n=1:niter
	sim.lennardjones("OO", [2.5, 1.0, 1.0, 1.0]);   
	sim.sfcoulomb(cutoff);
	sim.harmonicbond();
	sim.cossqangle();

	ekin = sim.leapfrog();
	sim.atoms.scalevels(temp0);
	
	if rem(n, 50)==0
		# Compress
		if sim.natoms/sim.volume < dens0 
			sim.scalebox(dens0, [1:3], 0.9995);
		end

		printf("\r It. no. %d, dens. %f, temp. %f ", ... 
			n, sim.natoms/sim.volume, ekin*2/(3*sim.natoms)); 
		fflush(stdout);
	end

end
printf("\n");

sim.atoms.save("water-sys.xyz"); sim.atoms.save("water-sys.mat");

