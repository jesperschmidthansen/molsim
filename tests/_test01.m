
%%%% Test 1 LJ

file=fopen("test01.log", "w");
printf(file, "\n --- Test 1: LJ conservation --- \n \n");

dens0 = 0.8; temp0=1.2;
lbox = (1000/dens0)^(1/3);

molsim('set', 'lattice', [10 10 10], [lbox lbox lbox]);
  
[epot, ekin, sum_mom]=_lj(10000, temp0, false, false);

fprintf(file, "\n  *Result*: "); 
fprintf(file, "Single CPU; energy %.5e +/- %.5e, momentum %.4e \n\n", ...
	mean(epot+ekin), std(epot+ekin), mean(sum_mom(:,1)));

t = linspace(0, length(epot), length(epot));
plot(t, epot, 'k', t, ekin, 'k', t, epot + ekin, 'b');

[epot, ekin, sum_mom]=_lj(10000, temp0, false, true);

fprintf(file, "\n  *Result*: ");
fprintf(file, "Multi CPU; energy %.5e +/- %.5e, momentum %.4e \n", ...
	mean(epot+ekin), std(epot+ekin), mean(sum_mom(:,1)));

fclose(file);
