
function test_1()
	 
	addpath("../mfiles/"); addpath("../mex/");

	T0 = 1.1;
	niter = 1e6;

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
			Pnow = (Pconf + Pkin)./p.volume();
			P = P + Pnow;	Pxy(counter) = Pnow(1,2); 
			counter++;
		end
	end

	T = 2/3*mean(Ekin(end-100:end))./p.natoms; 
	P = 10/niter*P;
	css = molsim_calccf(Pxy', Pxy', length(Pxy)/200).*p.volume(); 
	tcss = linspace(0, 0.005*10*length(css), length(css));
	eta0 = trapz(tcss, hann(css))./T0;

	subplot(2,1,1);
	index = [1:counter-1]; plot(index, Ekin./p.natoms, ";ekin;");
	subplot(2,1,2); 
	plot(tcss, css, "-o;Stress corr;");
	print("test_1.pdf", '-dpdf');

	printf("test_1 output:\n");
	printf("Ekin: %.3f +/- %.3f  Norm. mean temperature %f\n", mean(ekin)./p.natoms, std(ekin)./p.natoms, T/T0);
	printf("Av. normal press: %.3f  Visc. %.2f\n", 	(P(1,1)+P(2,2)+P(3,3))/3, eta0);
	
	volume = p.volume()
	save shearpressure.mat Pxy T0 volume
end
