

test_files = glob("test_*.m");
for n=1:length(test_files)
	run(test_files{n,1});
end
