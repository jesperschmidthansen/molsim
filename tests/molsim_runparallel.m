%%
%% molsim_runparallel: Example of simple Lennard-Jones simulation
%%
%% usage: execution_time = molsim_runparallel(ndir, density, num threads);
%%
%%        where ndir is the number of particles per direction
%%

function t= molsim_runparallel(nxyz, density, nthreads)

  cutoff = 2.5; epsilon = 1.0; sigma = 1.0; aw=1.0;

  lbox =  (nxyz^3/density)^(1.0/3.0);
  
  molsim('set', 'lattice', [nxyz nxyz nxyz], [lbox lbox lbox]);
  molsim('load', 'xyz', 'start.xyz');

  molsim('set', 'omp', nthreads);

  tic;
  for n=1:1000

    molsim('reset');
    molsim('calcforce', 'lj', 'AA', cutoff, sigma, epsilon, aw);
    molsim('integrate', 'leapfrog');

    if rem(n, 100)==0
      molsim('print');
    end

  end
  t = toc();

  molsim('clear');
  fprintf('\n');
  
end
