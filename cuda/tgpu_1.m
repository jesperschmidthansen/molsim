#########################################
##
## Scanity check programme 
##
##########################################

clear

nvt = true;
cutoff = 2.5; epsilon = 1.0; sigma = 1.0; aw = 1.0;

cmolsim('load', 'xyz', 'start_singleAN1000.xyz');
cmolsim('set', 'momresetfrq', 1000);

for n=1:100000
	
	cmolsim('reset');
	cmolsim('calcforce', 'lj', 'AA', cutoff, sigma, epsilon, aw);
	if nvt
		cmolsim('thermostat', 'nosehoover', 'A', 2.0, 0.1);
	endif	
	cmolsim('integrate', 'leapfrog');

	if rem(n,1000)==0 
		pressure = cmolsim('get', 'pressure');
		energies = cmolsim('get', 'energies');
		printf("%d %f %f %f %f %f\n", ... 
			n,  pressure(1), energies(1), energies(2),  sum(energies), energies(1)*2/3);
		fflush(stdout);
	endif
	
end

cmolsim('save', 'test.xyz');

% Free memory allocated
cmolsim('clear');
