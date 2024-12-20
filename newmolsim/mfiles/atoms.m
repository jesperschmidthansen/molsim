
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

				this.resetmom();
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
				
			ms_tether(this.f, this.r, this.rl, ptype, this.t, kspring, this.lbox, this.natoms);
			
		endfunction

		function mvlattice(this, ptype, dr)
			
			ms_mvlattice(this.rl, ptype, dr, this.t, this.lbox, this.natoms); 
			
		endfunction
	
		function mom = getmom(this)
			
			for k=1:3; mom(k) = sum(this.v(:,k).*this.m); end
		
		endfunction
	
		function resetmom(this)

			mom = this.getmom()./sum(this.m);
			for k=1:3; this.v(:,k) = this.v(:,k) - mom(k); end

		endfunction

		function setvels(this, temperature)
			
			this.v = randn(this.natoms, 3);
		
			ekin = 0.0;
			for k=1:3
				ekin = ekin + 0.5*(this.m'.*this.v(:,k)')*this.v(:,k);
			end

			temp = 2/3*ekin/this.natoms;			
			scale = sqrt(temperature/temp); 
			this.v = this.v*scale;
			
			this.resetmom();
		endfunction
		
		function vol = volume(this)
			vol = this.lbox(1)*this.lbox(2)*this.lbox(3);
		endfunction

		function lbox = getlbox(this)
			lbox = this.lbox;
		end

		function dist = getdist(this, i, j)
			
			dr = this.r(i,:)-this.r(j,:);
			for k=1:3
				if dr(k) > 0.5*this.lbox(k)
					dr(k) = dr(k) - this.lbox(k);
				elseif dr(k) < -0.5*this.lbox(k)
					dr(k) = dr(k) + this.lbox(k);
 				end
			end

			dist = sqrt(dot(dr,dr));
		end

	endmethods

end

