
clear all;

addpath("../mfiles/"); addpath("../mex/");

niter = 1e4;

p = atoms("start.xyz"); 
intgr = integrator(); 
prfrc = prforce();

ekin = zeros(1, niter); epot = zeros(1, niter);
p.m(1:2:end)=1.32;

tic();
for n=1:niter

	epot(n) = prfrc.lj(p, "AA", [2.5, 1.0, 1.0, 1.0]);   
	ekin(n) = intgr.step(p, prfrc);

end

t = toc(); sps = n/t;
printf("%f\n", sps);

ekin = ekin./p.natoms; epot = epot./p.natoms;
index = [1:n];
plot(index, ekin, ";ekin;", index, epot, ";epot;", index, epot+ekin, ";etot;")

