
%%
%% molsim_runlj: Simple Lennard-Jones simulation
%%
%% Test: 
%%

function molsim_runlj()

  cutoff = 2.5; epsilon = 1.0; sigma = 1.0; aw=1.0;

  molsim('set', 'lattice', [10 10 10], [12 12 12]);
  molsim('load', 'xyz', 'start.xyz');

  for n=1:10000

    molsim('reset');
    molsim('calcforce', 'lj', 'AA', cutoff, sigma, epsilon, aw);
    molsim('integrate', 'leapfrog');

    if rem(n, 100)==0
      molsim('print');
    end

  end

  molsim('save', 'A', 'start.xyz');
  molsim('clear');
  fprintf('\n');

  
  
end
