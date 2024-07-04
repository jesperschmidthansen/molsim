clear all;

addpath("../mfiles/"); addpath("../mex/");

p = atoms("start.xyz"); 
intgr = integrator(); 
inter = interactions();

for n=1:1000

	epot(n) = inter.lj(p); 
	ekin(n) = intgr.step(p);
end

ekin = ekin./p.natoms; epot = epot./p.natoms;
index = [1:n];
plot(index, ekin, ";ekin;", index, epot, ";epot;", index, epot+ekin, ";etot;")
