%
% usage; molslitconf(xyzFile, topFile, densFluid, height, numberMol, atype, lbond)
%	
%  xyzFile: Single molecule xyz file
%  topFile: Single molecule topology file
%  densFluid: Desired molecular fluid density
%  height: Desired slitpore height (confining direction)
%  numberMol: Number of molecules 
%  atype: Atomic type (character in molecule) OPTIONAL - defaults to 'A'
%  lbond: Bond length OPTIONAL - defaults to 1.0
%

function molslitconf(xyzFile, topFile, densFluid, height, numberMol, atype, lbond)

		if nargin < 5 || nargin > 7 
				error("Number of arguments not correct")
		end	
				
		if nargin == 5
			lbond = 1.0;
			atype = 'A';
		elseif nargin == 6 
			lbond = 1.0;
		end	
		

		function [npart str] = readheader(fin, opt)

		  npart = fscanf(fin, "%d\n", "C");

		  if strcmp(opt, "string")
		    str = fgets(fin, 1024);
		  else 
		    str = '\0';
		  endif

		endfunction

		function add_wallforce(lbox)
		  
		  x = molsim('get', 'positions');
		  
		  f = 1./x(:,3).^12 - 1./(x(:,3)-lbox).^12;
		  
		  molsim('add', 'force', f, 3);
		  
		endfunction


		function write_config(lbox, boff)
			
			nxy = int32(lbox(1));
			if ! rem(nxy,2)==0
				nxy = nxy - 1;
			endif
		
			str=sprintf("sep_lattice -b -n=%d,%d,4 -l=%f,%f,4.0 -f=wall.xyz", nxy, nxy, lbox(1), lbox(2))

			[status, output] = system(str);

			fin_mol = fopen("molecules.xyz", 'r');
			fin_wall = fopen("wall.xyz", 'r');

			npart_wall = readheader(fin_wall, 'string');
			npart_mol = readheader(fin_mol, 'string');

			fout = fopen("slitpore.xyz", 'w');

			fprintf(fout, "%d\n", npart_wall*2+npart_mol);
			fprintf(fout, "%f %f %f\n", lbox(1), lbox(2), lbox(3)*2);

			offset = 3.0 + boff; zmax = 0.0;
			for n=1:npart_mol
			[t, x, y, z, vx, vy, vz, m, q] = fscanf(fin_mol, "%c %f %f %f %f %f %f %f %f\n", "C");
			fprintf(fout, "%c %f %f %f %f %f %f %f %f\n",t, x, y, z+offset, vx, vy, vz, m, q);

			if z > zmax
			  zmax = z;
			end

			endfor

			offset = zmax + 2*boff + 3.0;

			for n=1:npart_wall
				[t, x, y, z, vx, vy, vz, m, q] = fscanf(fin_wall, "%c %f %f %f %f %f %f %f %f\n", "C");
				fprintf(fout, "w %f %f %f %f %f %f %f %f\n", x, y, z, vx, vy, vz, m, q);
				fprintf(fout, "W %f %f %f %f %f %f %f %f\n", x, y, offset+z, vx, vy, vz, m, q)
			endfor

			fclose(fin_mol);fclose(fin_wall);fclose(fout);

			printf("Wrote configuration in slitpore.xyz\n");
			printf("Topology file is start.top\n");
	  		printf("Wall density is set to %.3f\n", nxy^2*4/(lbox(1)^2*4.0));
			
		endfunction


		function Lbox = compress(atype, lb, numbermol, dens0, lbox_z, xyzfile, topfile)

		  temp0 = 4.0; ks = 1000.0;
		  types = [atype, atype];
			
		  molsim('set', 'molconfig', xyzfile, topfile, numbermol, 0.01, int32(rand*100));
		  molsim('load', 'xyz', 'start.xyz');
		  molsim('load', 'top', 'start.top');

		  molsim('set','timestep', 0.001);
		  molsim('set', 'temperature', temp0);
		  molsim('set', 'exclusion', 'molecule');
		  molsim('set', 'compressionfactor', 0.99995);

		  npart = molsim('get', 'numbpart');

		  Lbox = molsim('get', 'box');
		  lbox_xy = sqrt(npart/(lbox_z*dens0)); 
		  
		  while ( Lbox(1) > lbox_xy || Lbox(3) > lbox_z )

		    molsim('reset')
		    
		    molsim('calcforce', 'lj', types, 2.5, 1.0, 1.0, 1.0);
		    
		    molsim('calcforce', 'bond', 0, lb, ks);
		    
		    add_wallforce(molsim('get', 'box')(3));
		    
		    molsim('integrate', 'leapfrog')
		    molsim('thermostat', 'relax', atype, temp0, 0.01);
		    
		    molsim('compress', lbox_xy, 1);
		    molsim('compress', lbox_xy, 2);
		    molsim('compress', lbox_z, 3);

		    Lbox = molsim('get', 'box');
		  end

		  molsim('save', atype, 'molecules.xyz');
		  molsim('clear');
		  
		end

		lbox = compress(atype, lbond, numberMol, densFluid, height, xyzFile, topFile);
		write_config(lbox, 1.5)

		clear all;
endfunction
