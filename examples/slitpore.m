clear all;

cutoff = 2.5; epsilon = 1.0; sigma = 1.0; aw=1.0;
kspring = 250;
nloops = 1e6;

molsim('load', 'xyz', 'slitpore.xyz');

molsim('set', 'virtualsites');

molsim('sample', 'profiles', 'F', 200, 50);

for n=1:nloops
	% Reset forces etc
	molsim('reset');
	
	molsim('calcforce', 'lj', 'FF', cutoff, sigma, epsilon, aw);
	molsim('calcforce', 'lj', 'WW', cutoff, sigma, epsilon, aw);
	molsim('calcforce', 'lj', 'WF', cutoff, sigma, epsilon, aw);

	molsim('calcforce', 'lattice', 'W', 250);

	molsim('add', 'force', 'F', [0.01, 0.0, 0.0]);

	% Integrate forward in time - use leapfrog alogrithm
	molsim('integrate', 'leapfrog');
	molsim('thermostat', 'relax', 'W', 1.5, 0.1),

	molsim('sample', 'do');	

	if rem(n, 10000)==0
		data = load("profs.dat");
		dens = [data(20:end,3); data(1:19,3)];	temp = [data(20:end,4); data(1:19,4)];
 		vel = [data(20:end,5); data(1:19,5)]; 
		z = linspace(data(1,1), data(end,1), length(data));
		subplot(3,1,1); plot(z, dens); subplot(3,1,2); plot(z, temp); subplot(3,1,3); plot(z, vel);
		pause(0.01);
		printf("\r Fininshed %2.1f %%", n/nloops*100); fflush(stdout);
	end
end

% Free memory allocated
molsim('clear');
