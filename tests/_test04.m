

printf("\n --- Test 4: Molecular simulations --- \n \n"); fflush(stdout);

_water();

load water.mat;

printf("Dipole: %1.2f p/m %1.2f  (2.39 \pm 0.16 Wu et al.) \n", ...
       mean(dipoles)*0.584, std(dipoles));

printf("Bond length: %1.2f p/m %1.2f  (0.316) \n", mean(lbonds), std(lbonds));

printf("Angels: %1.2f p/m %1.2f (1.97)\n\n", mean(angels), std(angels));




