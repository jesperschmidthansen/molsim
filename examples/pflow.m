clear all;

addpath("../inst/"); addpath("../src/"); 

# Number of MD loops, cutoff, desired fluis density. A WCA fluid
nloops = 1e5; cutoff = 2.0^(1/6); rhoFluid = 0.7;

# Instance of molsim object and setting cutoff
sim = molsim();
sim.pairforce.maxcutoff = cutoff;

# Set number of particles and simulation box lengths in each direction. 
# Temperature is initially set to 1.0
sim.setconf([10,10,13], [10, 10, 13], 1.0);

# Set atom types: 'F' for fluid and 'W' for wall
sim.atoms.t(1:end) = 'F';
idxW = find( sim.atoms.r(:,3) > 9.0 ); sim.atoms.t(idxW) = 'W';

# Set dummy particles such that density is rhoFluid
idxF = find( sim.atoms.t == 'F'); nFluid = length(idxF); 
nDesired = int64(nFluid*rhoFluid);  nDummy = nFluid-nDesired;  sim.atoms.t(1:nDummy) = 'D';

# Get the final indices for fluids particles. Used later to add force
idxF = find( sim.atoms.t == 'F'); nFluid = length(idxF); 

# Thermostat
sim.setthermostat('relax', 'W', 1.0, 0.01);

# Main MD loop
for n=1:nloops
	# Calculate the forces acting between particles 
	sim.lennardjones("FF", [cutoff, 1.0, 1.0, 1.0]);   
	sim.lennardjones("WF", [cutoff, 1.0, 1.0, 1.0]);   
	sim.lennardjones("WW", [cutoff, 1.0, 1.0, 1.0]);

	# Aply spring force to wall atoms   
	sim.atoms.tether('W', 500);

	# Apply external force
	sim.atoms.f(idxF, 1) += 0.01;	
	
	# Apply thermostat	
	sim.applythermostat();

	# Integrate forward in time using the leap-frog algorithm
	sim.leapfrog();
end

