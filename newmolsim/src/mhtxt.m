#
# Help text generator
#
files = glob("*.c");

for n=1:length(files)
	fname = files{n};
	nfname = sprintf("%s.m", fname(1:end-2));	

	fid = fopen(nfname, "w");
	fprintf(fid, "#%s This function is a part of the molsim package\n", fname(1:end-2));
	fprintf(fid, "# It is likely a mex function for internal package use and there is currently no help text\n");
	fclose(fid);
end
