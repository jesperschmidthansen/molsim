clear all;

# Remove if script is run after package installation 
addpath("../inst/"); addpath("../src/");

# Number of MD loops
nloops = 1e4;

# Instance of molsim object
sim = molsim();

# Set number of particles and simulation box lengths in each direction. 
# Temperature is initially set to 1.0
sim.setconf([10,10,10], [10.557, 10.557, 10.557], 1.0);

# Main MD loop
ekin = zeros(nloops,1); epot = zeros(nloops,1);
for n=1:nloops
	# Calculates the forces acting on particles of type A - default type
	# Lennard-Jones parameters are [cutoff, epsilon, sigma, aW]
	epot(n) = sim.lennardjones("AA", [2.5, 1.0, 1.0, 1.0]);   

	# Integrate forward in time using the leap-frog algorithm
	ekin(n) = sim.leapfrog();
end

# Plot the energies
plot([1:nloops], ekin/sim.natoms, [1:nloops], epot/sim.natoms, [1:nloops], (epot+ekin)./sim.natoms);


