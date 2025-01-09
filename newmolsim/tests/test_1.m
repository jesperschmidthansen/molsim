
function test_1()
	 
	addpath("../mfiles/"); addpath("../mex/");

	T0 = 1.1;
	niter = 1e4;
	
	p = atoms([10,10,10], [11, 11, 11], T0);
	intgr = integrator(); 
	prfrc = prforce();
	thmstat = thermostat(p, T0);

	p.setvels(T0);

	counter = 1; P = zeros(3,3);
	for n=1:niter
		[epot, Pconf] = prfrc.lj(p, "AA", [2.5, 1.0, 1.0, 1.0]);   
		thmstat.nosehoover(p);
		[ekin Pkin] = intgr.step(p, prfrc);

		if rem(n, 10)==0 
			Ekin(counter) = ekin;
			P = P + Pconf + Pkin;	
			counter++;
		end

	end

	T = 2/3*mean(Ekin(end-100:end))./p.natoms; 
	P = P/(counter-1);

	index = [1:counter-1]; plot(index, 2/3*Ekin./p.natoms, ";T;");
	print("test_1.pdf", '-dpdf');

	printf("test_1 output:\n");
	printf("Ekin: %.3f +/- %.3f  Norm. mean temperature %f\n", mean(ekin)./p.natoms, std(ekin)./p.natoms, T/T0);
	printf("Av. normal press: %.3f  \n", (P(1,1)+P(2,2)+P(3,3))/3);

end
