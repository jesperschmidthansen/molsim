
function test_3()
	 
	addpath("../mfiles/"); addpath("../mex/");

	T0 = 1.42;
	niter = 1e4;

	p = atoms("start.xyz"); 
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
	
end
