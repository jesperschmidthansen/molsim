clear all;

addpath("../mfiles/"); addpath("../mex/");

p = atoms("start.xyz"); 
intgr = integrator(); 
inter = interactions();

tic();
for n=1:10000

	p.f = zeros(p.natoms, 3);
	if rem(n-1, 10)==0
		neighblist(p.nblist, p.r, p.lbox, 2.5, 0.5, p.natoms);
	end

	epot(n) = lj(p.f, "AA", [2.5, 1.0, 1.0, 1.0], p.r, p.t, p.nblist, p.lbox, p.natoms);   
	ekin(n) = intgr.step(p);

end
t = toc(); sps = n/t;
printf("%f\n", sps);

ekin = ekin./p.natoms; epot = epot./p.natoms;
index = [1:n];
plot(index, ekin, ";ekin;", index, epot, ";epot;", index, epot+ekin, ";etot;")

