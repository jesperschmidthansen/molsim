
printf("\n --- Test 3: Mem. management --- \n \n"); fflush(stdout);
  
dens0 = 0.8;
lbox = (1000/dens0)^(1/3);
molsim('set', 'lattice', [10 10 10], [lbox lbox lbox]);

pid=getpid();
str=sprintf("ps aux | grep %d > ps.out", pid);
system(str);
system("awk '{ print $4 }' ps.out > mem.out");
load mem.out; mem_0 = max(mem);

for n=1:1000
  molsim('load', 'xyz', 'start.xyz');  
  molsim('clear');
end
  
str=sprintf("ps aux | grep %d > ps.out", pid);
system(str);
system("awk '{ print $4 }' ps.out > mem.out");
load mem.out; mem_1 = max(mem);
  
printf("\n  *Result*: ");  
printf("Mem. before test %f after test %f\n", mem_0, mem_1);

  
