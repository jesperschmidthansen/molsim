
file = fopen("test04.m", "w");
fprintf(file, "\n --- Test 4: Molecular simulations --- \n \n"); fflush(stdout);

_water();

load water.mat;

fprintf(file, "Dipole: %1.2f p/m %1.2f  (2.39 \pm 0.16 Wu et al.) \n", ...
	mean(dipoles)*0.584, std(dipoles));

fprintf(file,"Bond length: %1.2f p/m %1.2f  (0.316) \n", mean(lbonds), std(lbonds));

fprintf(file,"Angels: %1.2f p/m %1.2f (1.97)\n\n", mean(angels), std(angels));

fclose(file);



