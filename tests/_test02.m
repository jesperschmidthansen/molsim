
printf("\n --- Test 2: LJ structure --- \n \n"); fflush(stdout);

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
printf("\n  *Result*: ");  
printf("Rdf maximum difference: %f %%  \n", (a-b)/b.*100);
  
plot(rdf(:,1), rdf(:,2), 'k-;molsim;', ...
     rdf_ref(:,1), rdf_ref(:,2), 'bo;Mosali et al;');
xlabel('r'); ylabel('rdf');
print('rdf.eps', '-deps');

  
