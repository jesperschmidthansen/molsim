
file = fopen("test05.log", "w");
fprintf(file,"\n --- Test 5: Effect of memory clearence --- \n \n"); fflush(stdout);

tic();
_lj(2000, 2.0, false, false);
x=molsim('get', 'positions');
v=molsim('get', 'velocities');
t=molsim('get', 'types');
_lj(2000, 2.0, false, false);
timer_1=toc();

save tmp.mat timer_1;

clear all

tic();
_lj(2000, 2.0, false, false);
x=molsim('get', 'positions');
v=molsim('get', 'velocities');
t=molsim('get', 'types');
clear all
_lj(2000, 2.0, false, false);
timer_2=toc();

load tmp.mat;

fprintf(file,"Timers: %f sec.  %f sec. \n\n", timer_1, timer_2);
fclose(file);



