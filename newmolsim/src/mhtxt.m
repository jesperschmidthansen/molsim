#
# Script to generate help txt - rather explaine why there is no help txt
#
files = glob("*.c");

for n=1:length(files)
	fname = files{n};
	nfname = sprintf("%s.m", fname(1:end-2))	

	fid = fopen(nfname, "w");
	fprintf(fid, "# This function is a part of the molsim package\n");
	fprintf(fid, "# It is likely a mex function for internal package use and there is currently no help text\n");
	fclose(fid);
end
