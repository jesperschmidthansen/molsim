clear all;

%%
%% _lj: Simple Lennard-Jones simulation (liquid, super crit. fluids, gasses)
%%
%% Test: 1) Energy conserveration (no thermostat)
%%       2) Momentum conservation
%%       3) Diffusion coefficient (Rowly & Painter)
%%       4) omp
%%
function [epot, ekin, sum_mom, D]=_lj(nloops, dens0, temp0, optThermostat, optOmp)
 cutoff = 2.5; epsilon = 1.0; sigma = 1.0; aw=1.0;

  lbox = (1000/dens0)^(1/3);
  molsim('set', 'lattice', [10 10 10], [lbox lbox lbox]);
  molsim('load', 'xyz', 'start.xyz');

  molsim('set', 'temperature', temp0);

  molsim('sample', 'vacf', 100, 5.0);
  
  if ( optOmp )
    molsim('set', 'omp', 2);
  end

  m = 1;
  for n=1:nloops

    molsim('reset');
    molsim('calcforce', 'lj', 'AA', cutoff, sigma, epsilon, aw);
    molsim('integrate', 'leapfrog');

    if ( optThermostat )
      molsim('thermostat', 'relax', 'A', temp0, 0.01);
    end

    if ( rem(n,100) == 0 )
      energies(m,:) = molsim('get', 'energies');

      v = molsim('get', 'velocities');
      sum_mom(m, :) = [sum(v(:,1)), sum(v(:,2)), sum(v(:,3))];

      m=m+1;
    endif

  end

  ekin = energies(:,1)./1000;
  epot = energies(:,2)./1000;

  data = load("vacf.dat");
  D = trapz(data(:,1), data(:,2));
  
  molsim('clear');
  fprintf('\n');

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Runing the tests %%%%%%%%%%%%%%%%%%%%%%%%%%%


testArray = [true, true];


%%%% Test 1 LJ
if ( testArray(1) )
    
  printf("\n --- Test 1: LJ conservation --- \n \n");
  fflush(stdout);
  
  dens0 = 0.8; temp0=1.2;
  
  [epot, ekin, sum_mom]=_lj(100000, dens0, temp0, false, false);
  
  printf("\n *Result*: ");
  printf("Single CPU; energy %.2e +/- %.2e, momentum %.4e \n\n", ...
	 mean(epot+ekin), std(epot+ekin), mean(sum_mom(:,1)));

  t = linspace(0, length(epot), length(epot));
  plot(t, epot, 'b', t, ekin, 'b', t, epot + ekin, 'b');
  
  if ( abs(mean(sum_mom)) > 1e-10 )
    printf("TEST FAILED! BAILING OUT");
    return
  end
       
  
  [epot, ekin, sum_mom]=_lj(100000, dens0, temp0, false, true);

  printf("\n *Result*: ")
  printf("Multi CPU; energy %.2e +/- %.2e, momentum %.4e \n\n", ...
	 mean(epot+ekin), std(epot+ekin), mean(sum_mom(:,1)));

  hold on
  plot(t, epot, 'b', t, ekin, 'b', t, epot + ekin, 'b');  
  hold off

  
  if ( abs(mean(sum_mom)) > 1e-10 )
    printf("TEST FAILED! BAILING OUT");
    return
  end
  
end

%%%% Test 2 LJ
if ( testArray(2) )
  
  printf("\n --- Test 2: LJ dynamics --- \n \n");
  fflush(stdout);
  
  dens0 = 0.8; temp0=1.2;
  
  [epot, ekin, sum_mom, D]=_lj(1000000, dens0, temp0, true, true);

  printf("\n *Result*: ")
  printf("Diffusion %f   (lit. 0.079 +/- 0.001)", D);
  fflush(stdout);

  data = load("vacf.dat");
  plot(data(:,1), cumtrapz(data(:,1), data(:,2)));
  
end
