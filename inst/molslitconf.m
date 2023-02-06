%
% usage; molslitconf(xyzFile, topFile, densFluid, height, numberMol, atype, lbond, bcc)
%	
%  xyzFile: Single molecule xyz file
%  topFile: Single molecule topology file
%  densFluid: Desired molecular fluid density
%  height: Desired slitpore height (confining direction)
%  numberMol: Number of molecules 
%  atype: Atomic type (character in molecule) OPTIONAL - defaults to 'A'
%  lbond: Bond length OPTIONAL - defaults to 1.0
%  bcc: set body centered crystal - default set to true 
%

function molslitconf(xyzFile, topFile, densFluid, height, numberMol, atype, lbond, bcc=true)

	if nargin < 5 || nargin > 8 
			error("Number of arguments not correct")
	end	
			
	if nargin == 5
		lbond = 1.0;
		atype = 'A';
	elseif nargin == 6 
		lbond = 1.0;
	end	
		

	function [read_npart read_str] = readheader(read_fin, read_opt)

		read_npart = fscanf(read_fin, "%d\n", "C");

		if strcmp(read_opt, "string")
			read_str = fgets(read_fin, 1024);
		else 
			read_str = '\0';
		endif

	endfunction

	function add_wallforce(add_lbox)
		
		add_x = molsim('get', 'positions');
		
		add_f = 1./add_x(:,3).^12 - 1./(add_x(:,3)-add_lbox).^12;
		
		molsim('add', 'force', add_f, 3);
		
	endfunction


	function write_config(w_lbox, w_boff, w_bcc)
		
		w_nxy = int32(w_lbox(1));

		if w_bcc
			if ! rem(w_nxy,2)==0
				w_nxy = w_nxy - 1;
			endif
	
			w_str=sprintf("sep_lattice -b -n=%d,%d,4 -l=%f,%f,4.0 -f=wall.xyz", ...
																w_nxy, w_nxy, w_lbox(1), w_lbox(2));
		else
			w_str=sprintf("sep_lattice -n=%d,%d,3 -l=%f,%f,3.0 -f=wall.xyz", ...
																w_nxy, w_nxy, w_lbox(1), w_lbox(2));
		endif
		
		[w_status, w_output] = system(w_str);

		w_fin_mol = fopen("molecules.xyz", 'r');
		w_fin_wall = fopen("wall.xyz", 'r');

		w_npart_wall = readheader(w_fin_wall, 'string');
		w_npart_mol = readheader(w_fin_mol, 'string');

		w_fout = fopen("slitpore.xyz", 'w');

		fprintf(w_fout, "%d\n", w_npart_wall*2+w_npart_mol);
		fprintf(w_fout, "%f %f %f\n", w_lbox(1), w_lbox(2), w_lbox(3)*2);

		w_offset = 3.0 + w_boff; w_zmax = 0.0;
		for w_n=1:w_npart_mol
			
			[w_t, w_x, w_y, w_z, w_vx, w_vy, w_vz, w_m, w_q] = ...
			fscanf(w_fin_mol, "%c %f %f %f %f %f %f %f %f\n", "C");
			fprintf(w_fout, "%c %f %f %f %f %f %f %f %f\n",w_t, w_x, w_y, w_z+w_offset, 
					w_vx, w_vy, w_vz, w_m, w_q);

			if w_z > w_zmax
				w_zmax = w_z;
			end

		endfor

		w_offset = w_zmax + 2*w_boff + 3.0;

		for w_n=1:w_npart_wall
			[w_t, w_x, w_y, w_z, w_vx, w_vy, w_vz, w_m, w_q] = ...
			 fscanf(w_fin_wall, "%c %f %f %f %f %f %f %f %f\n", "C");
			fprintf(w_fout,"w %f %f %f %f %f %f %f %f\n", w_x, w_y, w_z, w_vx, w_vy, w_vz, w_m,w_q);
			fprintf(w_fout,"W %f %f %f %f %f %f %f %f\n",w_x, w_y,w_offset+w_z,w_vx,w_vy,w_vz,w_m,w_q);
		endfor

		fclose(w_fin_mol);fclose(w_fin_wall);fclose(w_fout);

		printf("Wrote configuration in slitpore.xyz\n");
		printf("Topology file is start.top\n");
		printf("Wall density is set to %.3f\n", w_nxy^2*4/(w_lbox(1)^2*4.0));
		
	endfunction


	function c_Lbox = compress(c_atype, c_lb, c_numbermol, c_dens0, c_lbox_z, c_xyzfile, c_topfile)

		c_temp0 = 4.0; c_ks = 1000.0;
		c_types = [c_atype, c_atype];
		
		molsim('set', 'molconfig', c_xyzfile, c_topfile, c_numbermol, 0.01, int32(rand*100));
		molsim('load', 'xyz', 'start.xyz');
		molsim('load', 'top', 'start.top');
		printf("1\n"); fflush(stdout);
		
		molsim('set','timestep', 0.001);
		molsim('set', 'temperature', c_temp0);
		molsim('set', 'exclusion', 'molecule');
		molsim('set', 'compressionfactor', 0.99999);

		c_npart = molsim('get', 'numbpart');

		c_Lbox = molsim('get', 'box');
		c_lbox_xy = sqrt(c_npart/(c_lbox_z*c_dens0)); 
		
		while ( c_Lbox(1) > c_lbox_xy || c_Lbox(3) > c_lbox_z )

			molsim('reset')
		
			molsim('calcforce', 'lj', c_types, 2.5, 1.0, 1.0, 1.0);
		
			molsim('calcforce', 'bond', 0, c_lb, c_ks);
		
			add_wallforce(molsim('get', 'box')(3));
		
			molsim('integrate', 'leapfrog')
			molsim('thermostat', 'relax', c_atype, c_temp0, 0.001);
		
			molsim('compress', c_lbox_xy, 1);
			molsim('compress', c_lbox_xy, 2);
			molsim('compress', c_lbox_z, 3);

			Lbox = molsim('get', 'box');

			energies = molsim('get', 'energies');

			%molsim('print');
		end

		molsim('save', c_atype, 'molecules.xyz');
		molsim('clear');
		
	end

	lbox = compress(atype, lbond, numberMol, densFluid, height, xyzFile, topFile);
	write_config(lbox, 1.5, bcc);

	clear all;
endfunction
