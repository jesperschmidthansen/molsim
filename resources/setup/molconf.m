#
# Usage [nmols, atomsidx] = molconf(single mol xyz-file, single-mol top-file, numdim, offset,verbose)
#
# This is a function in the molsim package.
# 
# Generates a system configuration file and topology files from single molecule files. 
# Saves configuration in 'start.xyz' and topology files in 'bonds.top', 'angles.top' and
# 'dihedral.top' These files are overridden if they exist.
#
# Returns:
#  Number of molecules
#  Row-wise array of atom indices per molecule
#  
# Input:
#  single mol xyz-file: Molecule configuration  
#  single-mol top-file: Molecule topology
#  numdim: 3-array specifying number of molecules in x,y, and z directions.  
#  offset: Distance between molecular centre-of-mass. Default 10;
#  verbose: Set to true if some verbose statements should printed. Default false
#
# Examples: water.m and butane.m in example folder
# 
function [nmols, atom_idxs] = molconf(xyzfile, topfile, numdim, offset = 10.0, verbose=false)

	nmols = prod(numdim);
	
	### Positions
	[t, pos, mass, charge] = rpos(xyzfile);
	
	nuau = length(t);	

	# Wrt to cm and shift 
	molmass = sum(mass);
	for k=1:3 
		rcm(k) = sum(mass.*pos(:,k))/mass;
		pos(:,k) = pos(:,k) - rcm(k);
		pos(:,k) = pos(:,k) - min(pos(:,k));
	end

	if verbose
		printf("Found %d atoms in molecule - writing start.xyz file.\n", nuau);
	end

	wpos(t, pos, mass, charge, numdim, offset);

	### Bonds
	binfo = rbond(topfile);
	if binfo == -1 
		warning("No bonds found\n");
	else  
		if verbose
			printf("Found %d bonds in molecule - writing bonds.top file.\n", rows(binfo));
		end
		wbonds(numdim, nuau, binfo);
	end
	
	### Angles 
	ainfo = rangle(topfile);
	if ainfo != -1
		if verbose
			printf("Found %d angles in molecule - writing angles.top file.\n ", rows(ainfo));
		end
		wangles(numdim, nuau, ainfo);
	end

	### Dihedrals
	dinfo = rdihedral(topfile);
	if dinfo != -1
		if verbose
			printf("Found %d dihedrals in molecule - writing dihedrals.top file. \n ", rows(dinfo));	
		end	
		wdihedrals(numdim, nuau, dinfo);
	end
	
	### Mol. details to return 
	atom_idxs = zeros(nmols, nuau);
	for n=1:nmols
		for m=1:nuau
			atom_idxs(n,m) = (n-1)*nuau + m;
		end
	end

	save molinfo.mat nmols atom_idxs;

endfunction

function [t, pos, mass, charge] = rpos(xyzfile)

	### Position configuration	
	fin = fopen(xyzfile, "r");
	if fin == -1
		error("Couldn't open xyz-file");
	end
 
	nuau = fscanf(fin, "%d\n", "C");
	dummy = fgets(fin, 1024);

	for n=1:nuau
		[t(n), pos(n,1), pos(n,2), pos(n,3), mass(n), charge(n)] = ... 
												fscanf(fin, "%c%f%f%f%f%f\n", "C");
	end	

	fclose(fin);	

endfunction


function wpos(t, pos, mass, q, numdim, offset)

		
	nuau = length(t);
	nmol = prod(numdim);

	fout = fopen("start.xyz", "w");
	if fout == -1
		error("Couldn't open output file");
	end

	fprintf(fout, "%d\n%f %f %f\n", ... 
			nuau*nmol, numdim(1)*offset, numdim(2)*offset, numdim(3)*offset);


	idx=1; 
	for n=1:nmol
		for m=1:nuau
			amass(idx) = mass(m);
			idx++;
		end
	end
	avel = randn(nuau*nmol, 3);

	for k=1:3	
		mom = sum(amass'.*avel(:,k)); 
		avel(:,k) = avel(:,k) - mom/sum(amass);
	end
	
	idx=1;
	for n=1:numdim(3)
		for m=1:numdim(2)
			for k=1:numdim(1)
				for i=1:nuau
					xnuau = (k-1)*offset + pos(i,1); 
					ynuau = (m-1)*offset + pos(i,2); 
					znuau = (n-1)*offset + pos(i,3);
				
					fprintf(fout, "%c %f %f %f %f %f %f %f %f\n", ... 
			t(i), xnuau, ynuau, znuau, avel(idx,1), avel(idx,2), avel(idx,3), mass(i), q(i));				
					idx++;
				end
			end
		end
	end	

	# Momentum resetting missing!
	
	fclose(fout);

end



function bondinfo = rbond(topfile)

	fin = fopen(topfile, "r");
	if fin == -1
		error("Couldn't read top file");
	end

	readsection = false;
	nreadline = bidx = 0;
 
	while nreadline < 1e3

		strline = fgetl(fin);

		if strcmp(strline, "[ bonds ]")
			readsection = true;
			# Comment
			fgetl(fin);
				
			while true
				[mid, idx1, idx2, bt, count, err] = fscanf(fin, "%d %d %d %d\n", "C");		
				if count==4
					bidx++;
					bondinfo(bidx, :) = [idx1, idx2, bt]; 			
				else
					break;
				end			
			endwhile
		endif

		if readsection
			break;
		end

		nreadline++;
	endwhile

	if !readsection
		bondinfo = -1;
	end

	fclose(fin);

endfunction


function wbonds(numdim, nuau, bondinfo)

	fout = fopen("bonds.top", "w");
	if fout==-1
		error("Couldn't open bond top output file");
	end

	nmol = prod(numdim); nbonds = rows(bondinfo);
	for n=1:nmol
		for m=1:nbonds
			aidx = bondinfo(m, 1) + nuau*(n-1);
			bidx = bondinfo(m, 2) + nuau*(n-1);
			fprintf(fout, "%d %d %d %d\n", n, aidx, bidx, bondinfo(m,3));  
		end		
	end
	
	fclose(fout);
	
endfunction


function angleinfo = rangle(topfile)

	fin = fopen(topfile, "r");
	if fin == -1
		error("Couldn't read top file");
	end

	readsection = false;
	nreadline = aidx = 0;
 
	while nreadline < 1e3

		strline = fgetl(fin);

		if strcmp(strline, "[ angles ]")
			readsection = true;
			# Comment
			fgetl(fin);
				
			while true
				[mid, idx1, idx2, idx3, at, count, err] = fscanf(fin, "%d %d %d %d %d\n", "C");		
				if count==5
					aidx++;
					angleinfo(aidx, :) = [idx1, idx2, idx3, at]; 			
				else
					break;
				end			
			endwhile
		endif

		if readsection
			break;
		end

		nreadline++;
	endwhile

	if !readsection
		angleinfo = -1;
	end

	fclose(fin);

endfunction

function wangles(numdim, nuau, angleinfo)

	fout = fopen("angles.top", "w");
	if fout == -1
		error("Couldn't open bond top output file");
	end

	nmol = prod(numdim); nangles = rows(angleinfo);
	for n=1:nmol
		for m=1:nangles
			aidx = angleinfo(m, 1) + nuau*(n-1);
			bidx = angleinfo(m, 2) + nuau*(n-1);
			cidx = angleinfo(m, 3) + nuau*(n-1);
			fprintf(fout, "%d %d %d %d %d\n", n, aidx, bidx, cidx, angleinfo(m,4));  
		end		
	end
	
	fclose(fout);
	
endfunction


function dihedralinfo = rdihedral(topfile)

	fin = fopen(topfile, "r");
	if fin == -1
		error("Couldn't read top file");
	end

	readsection = false;
	nreadline = didx = 0;
 
	while nreadline < 1e3

		strline = fgetl(fin);

		if strcmp(strline, "[ dihedrals ]")
			readsection = true;
			# Comment
			fgetl(fin);
				
			while true
				[mid, idx1, idx2, idx3, idx4, dt, count, err] = ...
					fscanf(fin, "%d %d %d %d %d %d\n", "C");		
				if count==6
					didx++;
					dihedralinfo(didx, :) = [idx1, idx2, idx3, idx4, dt]; 			
				else
					break;
				end			
			endwhile
		endif

		if readsection
			break;
		end
		
		nreadline++;
	endwhile

	if !readsection
		dihedralinfo = -1;
	end

	fclose(fin);

endfunction

function wdihedrals(numdim, nuau, dihedralinfo)

	fout = fopen("dihedrals.top", "w");
	if fout == -1
		error("Couldn't open dihedral top output file");
	end

	nmol = prod(numdim); ndihedrals = rows(dihedralinfo);
	for n=1:nmol
		for m=1:ndihedrals
			aidx = dihedralinfo(m, 1) + nuau*(n-1);
			bidx = dihedralinfo(m, 2) + nuau*(n-1);
			cidx = dihedralinfo(m, 3) + nuau*(n-1);
			didx = dihedralinfo(m, 4) + nuau*(n-1);

			fprintf(fout, "%d %d %d %d %d %d\n", ...
					n, aidx, bidx, cidx, didx, dihedralinfo(m,5));  
		end		
	end
	
	fclose(fout);
	
endfunction


