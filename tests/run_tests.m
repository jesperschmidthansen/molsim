
clear all
system("make -C ../src/ -B");

test_files = glob("test_*.m");
for n=1:length(test_files)
	printf("---------> Running test %d ...\n", n-1); fflush(stdout);
	run(test_files{n,1});
end
