
%%
%% molsim_runbutane: Example of butane simulation
%% 

function _butane(opt)

	temp0 = 4.0;
	tau = 0.01;
	dt  = 0.002;

	bondlength = 0.4;
	bondconstant =  33615;
	bondangle = 1.9;
	angleconstant = 866.0;
	torsionparam = [15.5000,  20.3050, -21.9170, -5.1150,  43.8340, -52.6070];

	molsim('load', 'xyz', 'butane.xyz');
	molsim('load', 'top', 'butane.top');

	molsim('set','timestep', dt);
	molsim('set', 'temperature', temp0);

	if opt==1
		molsim('set', 'exclusion', 'molecule');
	else
		molsim('set', 'exclusion', 'bonded');
		molsim('set', 'omp', 2);
	end	
	
	for n=1:1000

		molsim('reset')

		molsim('calcforce', 'lj', 'CC', 2.5, 1.0, 1.0, 1.0);
		molsim('calcforce', 'bond', 0, bondlength, bondconstant);
		molsim('calcforce', 'angle', 0, bondangle, angleconstant);
		molsim('calcforce', 'torsion', 0, torsionparam);

		molsim('integrate', 'leapfrog')
		molsim('thermostat', 'relax', 'C', temp0, tau);

	end

	lbonds = molsim('get', 'bondlengths');
	angles = molsim('get', 'angles');
	torsions = molsim('get', 'torsions');	

	save butane.mat lbonds angles torsions;

	molsim('clear');
	fprintf('\n');  

end  

