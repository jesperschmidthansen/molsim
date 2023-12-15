# Water
clear

cutoff = 2.5; epsilon = 1.0; sigma = 1.0; aw = 1.0;
timestep = 5.0e-4;

cmolsim('load', 'xyz', 'start_water.xyz');
cmolsim('load', 'top', 'start_water.top');

cmolsim('set', 'timestep', timestep);
cmolsim('set', 'exclusion', 'molecule');
cmolsim('set', 'cutoff', 2.9); 

for n=1:10000
	
	cmolsim('reset');

	cmolsim('calcforce', 'lj', 'OO', cutoff, sigma, epsilon, aw);
	cmolsim('calcforce', 'coulomb', 'sf', 2.9);
	
	cmolsim('calcforce', 'bond', 0, 0.316, 68000);
	cmolsim('calcforce', 'angle', 0, 1.97, 490);
 		
	cmolsim('thermostat', 'nosehoover', 'O', 3.86, 0.1);
	cmolsim('integrate', 'leapfrog');

	if rem(n,1000)==0 
		energies = cmolsim('get', 'energies');
		printf("%d %f %f %f %f \n", n, energies(1), energies(2),  sum(energies), energies(1)*2/3);
		fflush(stdout);
	endif
	
endfor
	
