
classdef atoms < handle

	properties (Access=public)
		r, v, f;
		m, q, t; 
		nblist, max_nnb, update_nblist;
		boxcross, r0;
		rl;
		natoms,	lbox;
	end

	methods
		
		function this = atoms(filename)

			if ( nargin == 0 )
				this;
			elseif ( nargin==1 )
				if !exist(filename)
					error("Configuration file does not exists");
				end	
				
				format = filename(end-2:end);
	
				if strcmp(format, "xyz")
						
					fptr = fopen(filename, "r");
				
					natoms = fscanf(fptr, "%d\n", "C");  
					[Lx Ly Lz] = fscanf(fptr, "%f %f %f\n", "C");
					
					for n=1:natoms
						[t(n), x(n), y(n), z(n), vx(n), vy(n), vz(n), m(n), q(n)] = ...
													  fscanf(fptr, "%c %f %f %f %f %f %f %f %f\n", "C");
					endfor 	
		
					fclose(fptr);
				
					this.t = t'; this.m = m'; this.q = q';
					this.r = [x', y', z']; this.v = [vx', vy', vz']; this.f = zeros(natoms, 3);
					this.lbox = [Lx, Ly, Lz]; 
					this.natoms = natoms; 
					
					this.boxcross = int32(zeros(natoms, 3));
					this.update_nblist = true;
					this.r0 = [x', y', z']; this.rl = [x', y', z'];
					this.max_nnb = 500; this.nblist = -1*int32(ones(natoms, this.max_nnb)); ## ACHTUNG WITH 500
				elseif strcmp(format, "mat")
					load(filename);					
					this.r=r; this.v=v; this.f=f; this.m=m; this.q=q; this.t=t; 
					this.rl = rl; this.natoms = natoms; this.lbox = lbox;
	
					this.boxcross = int32(zeros(natoms, 3));
					this.update_nblist = true;
					this.r0 = r; #[x', y', z'];
					this.max_nnb = 500; this.nblist = -1*int32(ones(natoms, this.max_nnb)); ## ACHTUNG WITH 500
				else
					error("Format not supported");
				end
      		endif

		endfunction	
	
		function save(this, filename, write_opt="w")
	
			format = filename(end-2:end);

			if ( nargin > 2 )
				error("save takes maximum two arguments");
			end

			if strcmp(format, "xyz")
				fptr = fopen(filename, write_opt);
			
				fprintf(fptr, "%d\n", this.natoms);  
				fprintf(fptr, "%f %f %f\n", this.lbox(1), this.lbox(2), this.lbox(3));	
				for n=1:this.natoms
					fprintf(fptr, "%c %f %f %f %f %f %f %f %f\n", ...
							this.t(n), this.r(n,1), this.r(n,2), this.r(n,3), ... 
							this.v(n,1), this.v(n,2), this.v(n,3), this.m(n), this.q(n));
				endfor 	
	
				fclose(fptr);
			elseif strcmp(format, "mat")		
				r=this.r; v=this.v; f=this.f; m=this.m; q=this.q; t=this.t; 
				rl = this.rl; natoms=this.natoms;	lbox=this.lbox;
				
				save(filename, "r", "v", "f", "m", "q", "t", "rl", "natoms", "lbox");	
			else
				error("Format not supported");
			end
				
		endfunction

		function tether(this, ptype, kspring)
				
			tether(this.f, this.r, this.rl, ptype, this.t, kspring, this.lbox, this.natoms);
			
		endfunction

	endmethods

end

