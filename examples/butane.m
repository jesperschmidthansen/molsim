
clear all;

# Remove if script is run after package installation 
addpath("../inst/"); addpath("../src/");  addpath("../resources/setup/");

# Variables holding number of MD loops, time step size, simulation temperature and density 
nloops = 1e4; dt = 1e-3;
temp0 = 3.0; dens0 = 1.49;

# Write configuation and topology files: start.xyz bonds.top, angles.top, dihedrals.top  
# Configuration is 10x10x10 molecules with center-of-mass distance of 4
system("rm -f *top");
molconf("../resources/molconf/butane.xyz", "../resources/molconf/butane.top", [10 10 10], 4.0);

# Instance of molsim 
sim = molsim();

# Set configureation from start.xyz
sim.setconf("start.xyz");

# Set integrator time step and thermostat temperature
sim.integrator.dt = dt;
sim.thermostat.temperature= temp0;

# Reads in topology files 
sim.settop();

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
ekin = zeros(nloops,1); epot = zeros(nloops,1);
for n=1:nloops
	# Calculate the pair forces 	
	epot(n) = sim.lennardjones("CC", [2.5, 1.0, 1.0, 1.0]);   
	
	# Calcuate the intra-molecular forces
	epot(n) += sim.harmonicbond(0);
	epot(n) += sim.harmonicangle(0);
	epot(n) += sim.ryckbell(0);

	# Apply Nose-Hoover thermostat
	sim.nosehoover();

	# Integrate forward in time
	ekin(n) = sim.leapfrog();
	
	# Rescale simulation box to obtain desired density
	if rem(n, 10)==0
		sim.scalebox(dens0);
	end

end

# Save the final configuration
sim.atoms.save("final.xyz");

