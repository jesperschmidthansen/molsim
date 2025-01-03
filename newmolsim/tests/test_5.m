
function [dist _angles dihedrals] = test_5()

	addpath("../mfiles/"); addpath("../mex/");

	niter = 1e4;

	p = atoms("mol.xyz"); 
	intgr = integrator(0.001); 
	prfrc = prforce(); 

	# Molecular definitions: Each mol is four atoms long. 
	# Linear with three bonds, two angles, one dihedral
	# All this will later be automated using input files
	if rem(p.natoms,4) != 0
		error("test_5 failed - input file error");
	end
	nmols = int32(p.natoms/4);

	# Bonds 	
	nbonds = nmols*3; b = bonds(nbonds);
	pidx = 0;
	for n=1:nmols
		for m=1:3
			pidx++;
			b.pidx(3*(n-1)+m, :) = [pidx, pidx+1]; 
		end
		pidx = pidx + 1;
	end
	b.springs = 500*ones(nbonds,1); b.l0 = 1.1*ones(nbonds,1); 

	# Angles
	nangles = nmols*2; a = angles(nangles);
	pidx = 0;
	for n=1:nmols
		for m=1:2
			pidx++;
			a.pidx(2*(n-1)+m, :) = [pidx, pidx+1, pidx+2];
		end
		pidx = pidx + 2;
	end
	a.springs = 100*ones(nangles,1); a.a0 = 2.1*ones(nangles, 1);

	# Dihedrals
	ndihedrals = nmols; d = dihedrals(ndihedrals);
	pidx = 0;
	for n=1:nmols
		pidx++;
		d.pidx(n, :) = [pidx, pidx+1, pidx+2, pidx+3];
		pidx = pidx + 3;
	end
	d.coeffs = zeros(ndihedrals,6); d.coeffs(:, 2) = -30.0;

	# Exclude pair forces in molecule
	for n=1:4:p.natoms
		p.exclude(n,1:3) = [n+1, n+2, n+3]; 
		p.exclude(n+1, 1:3) = [n, n+2, n+3];
		p.exclude(n+2, 1:3) = [n, n+1, n+3];
		p.exclude(n+3, 1:3) = [n, n+1, n+2];
	end	

	ekin = zeros(1, niter); epot = zeros(1, niter); 
	dist = zeros(1, niter); _angles = zeros(1, niter); dihedrals = zeros(1, niter);
	
	p.resetmom();

	for n=1:niter
		epot(n) = prfrc.lj(p, "AA", [2.5, 1.0, 1.0, 1.0]);   
		epot(n) = epot(n) + b.harmonic(p, 0);
		epot(n) = epot(n) + a.harmonic(p, 0);
		epot(n) = epot(n) + d.ryckbell(p, 0);
		ekin(n) = intgr.step(p, prfrc);
		
		dist(n) = p.getdist(1,2);
		_angles(n) = p.getangle(1,2,3);
		dihedrals(n) = p.getdihedral(1,2,3,4);
	end
	
	ekin = ekin./p.natoms; epot = epot./p.natoms; etot = epot + ekin;
	spnb = niter/prfrc.neighb_updates; 
	mom = p.getmom();

	printf("test_5 output:\n");
	printf("Etot: %1.3e +/- %1.3e \n", mean(etot), std(etot));
	printf("Ekin %1.3e  Momentum: %e %e %e\n", mean(ekin), mom(1), mom(2), mom(3));
	printf("Steps per build  %.3f\n", spnb);
	printf("Av. bond length %.3f, eq. bond length set to %.3f \n", mean(dist), b.l0(1));
	printf("Av. angle %.3f, eq. angle set to %.3f \n", mean(_angles), a.a0(1));
	printf("Av. dihedral %.3f\n", mean(dihedrals));

end
