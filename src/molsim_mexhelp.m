
# Usage: molsim_mexhelp()
# 
# Generating generic help texts to mex-files

function molsim_mexhelp()

	cfiles = glob("*.c");

	for n=1:length(cfiles)
		
		filename = cfiles{n};
		mhfile = sprintf("%s.m", filename(1:end-2));
		
		fout = fopen(mhfile, "w");
		if fout == -1
			warning("Couldn't open m-help file - skipping...\n");
			pause(3);
			break;
		end

		fprintf(fout, "# %s is a part of the molsim package\n", filename);
		fprintf(fout, "# It is not inteded that it is used directly by the user, ");
		fprintf(fout, "# but only through the package class methods\n");
		fprintf(fout, "# For usage see the package source https://github.com/jesperschmidthansen/molsim");

		fclose(fout);
	end

end
