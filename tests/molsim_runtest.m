clear all;

%%
%% _lj: Simple Lennard-Jones simulation (liquid, super crit. fluids, gasses)
%%
%%

function [epot, ekin, sum_mom]=_lj(nloops, temp0, optThermostat, optOmp)
  cutoff = 2.5; epsilon = 1.0; sigma = 1.0; aw=1.0;
  
  molsim('load', 'xyz', 'start.xyz');
  molsim('set', 'temperature', temp0);

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
      #molsim('thermostat', 'nosehoover', 'A', temp0, 10.0);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

testArray = [true,  
	     true,
	     true];

log_file = fopen("test.log", "w");


%%%% Test 1 LJ
if ( testArray(1) )
  fprintf(stdout, "\n --- Test 1: LJ conservation --- \n \n");
  fflush(stdout);
  fprintf(log_file, "\n --- Test 1: LJ conservation --- \n \n");
  
  dens0 = 0.8; temp0=1.2;

  lbox = (1000/dens0)^(1/3);
  molsim('set', 'lattice', [10 10 10], [lbox lbox lbox]);
  
  [epot, ekin, sum_mom]=_lj(10000, temp0, false, false);

  fprintf(stdout, "\n  *Result*: "); fprintf(log_file, "\n  *Result*: ");
  printf("Single CPU; energy %.5e +/- %.5e, momentum %.4e \n\n", ...
	 mean(epot+ekin), std(epot+ekin), mean(sum_mom(:,1)));
  fprintf(log_file, "Single CPU; energy %.5e +/- %.5e, momentum %.4e \n\n", ...
	 mean(epot+ekin), std(epot+ekin), mean(sum_mom(:,1)));

  t = linspace(0, length(epot), length(epot));
  plot(t, epot, 'k', t, ekin, 'k', t, epot + ekin, 'b');
  
  [epot, ekin, sum_mom]=_lj(10000, temp0, false, true);

  fprintf(stdout, "\n  *Result*: "); fprintf(log_file, "\n  *Result*: ")
  printf("Multi CPU; energy %.5e +/- %.5e, momentum %.4e \n", ...
	 mean(epot+ekin), std(epot+ekin), mean(sum_mom(:,1)));
  fprintf(log_file, "Multi CPU; energy %.5e +/- %.5e, momentum %.4e \n", ...
	 mean(epot+ekin), std(epot+ekin), mean(sum_mom(:,1)));
  
end

%%%% Test 2 LJ Ref. Morsali et al., Chem. Phys., 310:11 (2005)
if ( testArray(2) )

  fprintf(stdout, "\n --- Test 2: LJ structure --- \n \n"); fflush(stdout);
  fprintf(log_file, "\n--- Test 2: LJ structure --- \n \n");
  
  dens0 = 0.9; temp0=1.5;

  lbox = (1000/dens0)^(1/3);
  molsim('set', 'lattice', [10 10 10], [lbox lbox lbox]);
  
  [epot, ekin, sum_mom]=_lj(10000, temp0, true, false);


  system("cp final.xyz start.xyz");
  [epot, ekin, sum_mom]=_lj(100000, temp0, true, false);
  
  rdf_ref = load("rdf_ref.dat");
  rdf = load("radial.dat"); rdf(:,2) = rdf(:,2)./rdf(end,2);

  a = max(rdf_ref(:,2));
  b = max(rdf(:,2));
  fprintf(stdout, "\n  *Result*: ");  fprintf(log_file, "\n  *Result*: ")
  fprintf(stdout, "Rdf maximum difference: %f %%  \n", (a-b)/b.*100);
  fprintf(log_file, "Rdf maximum difference: %f %%  \n", (a-b)/b.*100);
  
  plot(rdf(:,1), rdf(:,2), 'k-;molsim;', ...
       rdf_ref(:,1), rdf_ref(:,2), 'bo;Mosali et al;');
  xlabel('r'); ylabel('rdf');
  print('rdf.eps', '-deps');
end

if ( testArray(3) )

  fprintf(stdout, "\n --- Test 3: Mem. management --- \n \n"); fflush(stdout);
  fprintf(log_file, "\n --- Test 3: Mem. management --- \n \n");
  
  dens0 = 0.8;
  lbox = (1000/dens0)^(1/3);
  molsim('set', 'lattice', [10 10 10], [lbox lbox lbox]);

  pid=getpid();
  str=sprintf("ps aux | grep %d > ps.out", pid);
  system(str);
  system("awk '{ print $4 }' ps.out > mem.out");
  load mem.out; mem_0 = max(mem);

  for n=1:1000
    molsim('load', 'xyz', 'start.xyz');  
    molsim('clear');
  end
  
  str=sprintf("ps aux | grep %d > ps.out", pid);
  system(str);
  system("awk '{ print $4 }' ps.out > mem.out");
  load mem.out; mem_1 = max(mem);

  
  fprintf(stdout, "\n  *Result*: ");  fprintf(log_file, "\n  *Result*: ");
  fprintf(stdout, "Mem. before test %f after test %f\n", mem_0, mem_1);
  fprintf(log_file, "Mem. before test %f after test %f\n", mem_0, mem_1);
  
end

fclose(log_file);
