
filestr = fopen("test06.log", "w");
for n=0:1
	_butane(n);

	load butane.mat
	fprintf(filestr, "Run id %d: %f %f %f \n", n, mean(angles), mean(lbonds), mean(torsions));
end

fclose(filestr);

