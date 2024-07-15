
clear all;

addpath("../mfiles/"); addpath("../mex/");

niter = 1e4;

p = atoms("start.xyz"); 
intgr = integrator(); 
inter = interactions();

ekin = zeros(1, niter); epot = zeros(1, niter);

tic();
for n=1:niter

	p.f = zeros(p.natoms, 3);
	dr2(n) = calcdr2(p.r, p.r0, p.boxcross, p.lbox, p.natoms);
	if n==1 || dr2(n)>0.5*0.5 
		neighblist(p.nblist, p.r, p.r0, p.lbox, 2.5, 0.5, p.natoms);
	end

	epot(n) = inter.lj(p, "AA", [2.5, 1.0, 1.0, 1.0]);   
	ekin(n) = intgr.step(p);

end

t = toc(); sps = n/t;
printf("%f\n", sps);

ekin = ekin./p.natoms; epot = epot./p.natoms;
index = [1:n];
plot(index, ekin, ";ekin;", index, epot, ";epot;", index, epot+ekin, ";etot;")

