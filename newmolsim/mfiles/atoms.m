
classdef atoms < handle

	properties (Access=public)
		r, v, f;
		m, q, t; 
		nblist, max_nnb; 
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

				this.max_nnb = 500; this.nblist = -1*int32(ones(natoms, this.max_nnb));


      		endif

		end	
	
		function save(this, filename="config.mat")
			
			r = this.r; v = this.v; f=this.f; 
			m = this.m; q = this.q; t = this.t;

			save(filename, "r", "v", "f", "m", "q", "t");	
			
		end	
	
	end

end
