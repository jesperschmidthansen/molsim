
# Simulation of Ryckaert-Bellmann butane model 
# System configuration is load from the file butane_config_system.xyz and is calibrated  

clear all;

# Remove if script is run after package installation 
addpath("../inst/"); addpath("../src/");  

# Variables holding number of MD loops, time step size, simulation temperature and density 
nloops = 1e4; dt = 1e-3;
temp0 = 3.0; dens0 = 1.49;

# Instance of molsim 
sim = molsim();

# Set configureation from file 
sim.setconf("../resources/molconf/butane_config_system.xyz");

# Reads in topology  
sim.settop("../resources/molconf/butane_topology_system.mat");

# Set integrator time step and thermostat tauQ = 10.0 
sim.integrator.dt = dt;
sim.setthermostat("nh", temp0, 10.0);

# Sets the bond spring constant and zero force bond length - only one bond type, namely, 0
nbonds = sim.bonds.nbonds;
sim.bonds.springs = 33e3*ones(nbonds,1); sim.bonds.l0 = 0.4*ones(nbonds,1); 

# Sets the angle spring constant and zero force angle (radians) - only one angle type, namely, 0
nangles = sim.angles.nangles;	
sim.angles.springs = 866*ones(nangles,1); sim.angles.a0 = 1.9*ones(nangles, 1);

# Sets the dihedral coefficients for the Ryckert-Belleman mode - only one dihedral, namely, 0 
ndihedrals = sim.dihedrals.ndihedrals;
for n=1:ndihedrals
	sim.dihedrals.coeffs(n,:) = [15.5000,  20.3050, -21.9170, -5.1150,  43.8340, -52.6070]; 
end

# Exclude pair forces for atoms in same dihedral angle (here this is the same as the molecule) 
sim.atoms.setexclusions(sim.dihedrals.pidx, "dihedrals");

# Main MD loop
for n=1:nloops
	# Calculate the pair forces 	
	sim.lennardjones("CC", [2.5, 1.0, 1.0, 1.0]);   
	
	# Calcuate the intra-molecular forces
	sim.harmonicbond();
	sim.cossqangle();
	sim.ryckbell();

	# Apply Nose-Hoover thermostat
	sim.applythermostat();
	
	# Integrate forward in time
	sim.leapfrog();

	# Here you can call a sampler or collect data
end

