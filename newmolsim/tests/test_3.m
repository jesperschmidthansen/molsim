
function test_3()
	 
	addpath("../mfiles/"); addpath("../mex/");

	T0 = 1.42;
	niter = 1e4;

	p = atoms([10,10,10], [11, 11, 11], 2.0);
	p.t(1:500)='W';	p.t(501:end) = 'F';

	intgr = integrator(); 
	prfrc = prforce();
	thmstat = thermostat(p, T0, 'W');
	
	ekin = zeros(niter,1); 
	for n=1:niter
		prfrc.lj(p, "FF", [2.5, 1.0, 1.0, 1.0]); 
		prfrc.lj(p, "WW", [2.5, 1.0, 1.0, 1.0]); 
  		prfrc.lj(p, "FW", [2.5, 1.0, 1.0, 1.0]); 
		
		p.tether('W', 300);
	
		thmstat.relaxtemp(p);
		ekin(n) = intgr.step(p, prfrc);
	end

	p.save("tether.xyz");

	plot3(p.r(1:500,1), p.r(1:500,2), p.r(1:500,3), 'o', 'markerfacecolor', 'red');
	hold on
	plot3(p.r(501:end,1), p.r(501:end,2), p.r(501:end,3), 'o', 'markerfacecolor', 'blue');
	hold off
	view([0,0,0]);
	print("test_3.pdf", '-dpdf');

end
