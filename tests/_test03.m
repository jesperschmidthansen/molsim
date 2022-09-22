
file=fopen("test03.log", "w");

fprintf(file, "\n --- Test 3: Mem. management --- \n \n"); 
  
dens0 = 0.8;
lbox = (1000/dens0)^(1/3);
molsim('set', 'lattice', [10 10 10], [lbox lbox lbox]);

pid=getpid();
str=sprintf("ps aux | grep %d > ps.out", pid);
system(str);
system("awk '{ print $4 }' ps.out > mem.out");
load mem.out; mem_0 = max(mem);

for n=1:10000
  molsim('load', 'xyz', 'start.xyz');
  x = molsim('get', 'positions');
  t = molsim('get', 'types');
  molsim('clear');
end
  
str=sprintf("ps aux | grep %d > ps.out", pid);
system(str);
system("awk '{ print $4 }' ps.out > mem.out");
load mem.out; mem_1 = max(mem);

fprintf(file, "\n  *Result*: ");  
fprintf(file, "Mem. before test %f after test %f\n", mem_0, mem_1);

for n=1:1000
  _lj(110, 1.0, false, false);
end

str=sprintf("ps aux | grep %d > ps.out", pid);
system(str);
system("awk '{ print $4 }' ps.out > mem.out");
load mem.out; mem_1 = max(mem);
  
fprintf(file, "\n  *Result*: ");  
fprintf(file, "Mem. before test %f after test %f\n", mem_0, mem_1);

fclose(file);
  
