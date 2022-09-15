%%
%% _lj: Simple Lennard-Jones simulation (liquid, super crit. fluids, gasses)
%%
%%

function [epot, ekin, sum_mom]=_lj(nloops, temp0, optThermostat, optOmp)
  cutoff = 2.5; epsilon = 1.0; sigma = 1.0; aw=1.0;
  
  molsim('load', 'xyz', 'start.xyz');
 	 
  molsim('sample', 'vacf', 100, 5.0);
  molsim('sample', 'radial', 100, 50, 'A');

  if ( optOmp )
    molsim('set', 'omp', 2);
  end

  m = 1;
  for n=1:nloops
    molsim('reset');

    molsim('calcforce', 'lj', 'AA', cutoff, sigma, epsilon, aw);
    
    if ( optThermostat )
      molsim('integrate', 'leapfrog');
      molsim('thermostat', 'relax', 'A', temp0, 0.01);
    else
      molsim('integrate', 'leapfrog');
    end

    molsim('sample', 'do');

    if ( rem(n,100) == 0) 
      energies(m,:) = molsim('get', 'energies');

      v = molsim('get', 'velocities');
      sum_mom(m, :) = [sum(v(:,1)), sum(v(:,2)), sum(v(:,3))];

      m=m+1;
    end
    
  end

  ekin = energies(:,1)./1000;
  epot = energies(:,2)./1000;

  molsim('save', 'A', 'final.xyz');
  
  molsim('clear');
  fprintf('\n');
end

