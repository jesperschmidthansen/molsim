#
# Usage molslitpore(xyz file, wall thickness, whall density)
#
# This is a function in the molsim package.
# 

function molslitpore(xyzmol, wthickness, rho)


	### Open and read molecule configuration file
	fin	= fopen(xyzmol, "r");
	if fin==-1
		error("Couldn't open xyz file for molecules");
	end
	
	natomsMols = fscanf(fin, "%d\n", "C");
	[lbox(1), lbox(2), lbox(3)] = fscanf(fin, "%f%f%f\n", "C");

	for n=1:natomsMols 
		[t(n), pos(n,1), pos(n,2), pos(n,3), vel(n,1), vel(n,2), vel(n,3), mass(n), charge(n)] = ... 
			fscanf(fin, "%c%f%f%f%f%f%f%f%f\n", "C");
	end	

	fclose(fin);	
		

	### Build wall 1
	[t_w1 pos_w1 vel_w1 mass_w1 charge_w1] = wall('w', lbox, wthickness(1), rho(1));
	
	natomsWall1 = length(t_w1);

	### Build wall 2 
	[t_w2 pos_w2 vel_w2 mass_w2 charge_w2] = wall('W', lbox, wthickness(2), rho(2));
	
	natomsWall2 = length(t_w2);

	### Write augmented configuration file 

	fout = fopen("molslitpore.xyz", "w");
	if  fout == -1
		error("Couldn't open xyz output file\n");
	end

	natoms = natomsMols + natomsWall1 + natomsWall2;

	fout = fopen("molslitpore.xyz", "w");	
	fprintf(fout, "%d\n", natoms);
	fprintf(fout, "%f %f %f\n", lbox(1), lbox(2), lbox(3)+wthickness(1)+wthickness(2));
	
	# just copying	
	for n=1:natomsMols
		fprintf(fout, "%c %f %f %f %f %f %f %f %f\n", t(n), pos(n,1), pos(n,2), pos(n,3), ...
				vel(n,1), vel(n,2), vel(n,3), mass(n), charge(n));
	end	

	pos_w1(:,3) = pos_w1(:,3)+lbox(3);
	for n=1:natomsWall1
		fprintf(fout, "%c %f %f %f %f %f %f %f %f\n", t_w1(n), ...
		   	pos_w1(n,1), pos_w1(n,2),pos_w1(n,3), vel_w1(n,1), vel_w1(n,2), vel_w1(n,3), ...
			mass_w1(n), charge_w1(n));
	end

	pos_w2(:,3) = pos_w2(:,3)+lbox(3)+wthickness(1);
	for n=1:natomsWall2
		fprintf(fout, "%c %f %f %f %f %f %f %f %f\n", t_w2(n), ...
		   	pos_w2(n,1), pos_w2(n,2),pos_w2(n,3), vel_w2(n,1), vel_w2(n,2), vel_w2(n,3), ...
			mass_w2(n), charge_w2(n));
	end

	fclose(fout);

end


function [t pos vel mass charge] = wall(type, lbox, wallthickness, rho=1.0)

	dx = (1.0/rho)^(1/3);
	nx = int64(lbox(1)/dx); ny = int64(lbox(2)/dx); nz = int64(wallthickness/dx);
	natoms = nx*ny*nz;

	idx = 1;
	for n=1:nz
		for m=1:ny
			for k=1:nx
				pos(idx,1) = (k-1)*dx; 
				pos(idx,2) = (m-1)*dx; 
				pos(idx,3) = (n-1)*dx;
				t(idx) = type;
				idx++;
			end
		end
	end

	vel = randn(rows(pos), 3); mass=ones(1, rows(pos)); charge = zeros(1, rows(pos));

end
