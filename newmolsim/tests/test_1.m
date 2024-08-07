
function test_1()
	 
	addpath("../mfiles/"); addpath("../mex/");

	T0 = 1.1;
	niter = 1e4;

	p = atoms("start.xyz"); 
	intgr = integrator(); 
	prfrc = prforce();
	thmstat = thermostat(p, T0);

	p.setvels(T0);

	counter = 1; P = zeros(3,3);
	for n=1:niter
		[epot, Pconf] = prfrc.lj(p, "AA", [2.5, 1.0, 1.0, 1.0]);   
		thmstat.relaxtemp(p);
		[ekin Pkin] = intgr.step(p, prfrc);

		if rem(n, 10)==0
			Ekin(counter) = ekin;
			P = P + Pconf + Pkin;	
			counter++;
		end
	end

	T = 2/3*mean(Ekin(end-100:end))./p.natoms; 
	P = 10/niter*P;

	index = [1:counter-1]; plot(index, Ekin./p.natoms, ";ekin;");
	print("test_1.pdf", '-dpdf');

	printf("test_1 output:\n");
	printf("Ekin: %.3f +/- %.3f  Norm. mean temperature %f\n", mean(ekin)./p.natoms, std(ekin)./p.natoms, T/T0);
	printf("Av. normal press: %.3f  \n", (P(1,1)+P(2,2)+P(3,3))/3);

end
