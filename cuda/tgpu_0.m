##########################################
##
## Simple benchmark script
## system('export LD_LIBRARY_PATH=.');
##
##########################################

clear

tic();

cutoff = 2.5; epsilon = 1.0; sigma = 1.0; aw = 1.0;

cmolsim('load', 'xyz', 'start_singleAN1000.xyz');

for n=1:100000
	
	cmolsim('reset');
	
	cmolsim('calcforce', 'lj', 'AA', cutoff, sigma, epsilon, aw);

	cmolsim('integrate', 'leapfrog');

end

cmolsim('save', 'test.xyz');

% Free memory allocated
cmolsim('clear');

toc()
