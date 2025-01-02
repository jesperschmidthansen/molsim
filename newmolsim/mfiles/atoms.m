
classdef atoms < handle

	properties (Access=public)
		# Atom properties
		r, v, f;
		m, q, t; 
		# Neighbourlist etc
		nblist, max_nnb, update_nblist;
		# Simulation box crossing
		boxcross;
		# Initial and last positions
		r0; rl;
		# Pair interaction exclusion list
		exclude, max_exclude;
		# Number of atoms and simulation box length 
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
					end 	
		
					fclose(fptr);
				
					this.t = t'; this.m = m'; this.q = q';
					this.r = [x', y', z']; this.v = [vx', vy', vz']; this.f = zeros(natoms, 3);
					this.lbox = [Lx, Ly, Lz]; 
					this.natoms = natoms; 
					
					this.boxcross = int32(zeros(natoms, 3));
					this.update_nblist = true;
					this.r0 = [x', y', z']; this.rl = [x', y', z'];

					this.max_nnb = 500; this.nblist = -1*int32(ones(natoms, this.max_nnb)); ## ACHTUNG WITH 500
					this.max_exclude = 5; this.exclude = -1*int32(ones(natoms, this.max_exclude));

				elseif strcmp(format, "mat")
					load(filename);					
					this.r=r; this.v=v; this.f=f; this.m=m; this.q=q; this.t=t; 
					this.rl = rl; this.natoms = natoms; this.lbox = lbox;
	
					this.boxcross = int32(zeros(natoms, 3));
					this.update_nblist = true;
					this.r0 = r; #[x', y', z'];
				
					this.max_nnb = 500; this.nblist = -1*int32(ones(natoms, this.max_nnb)); ## ACHTUNG WITH 500
					this.max_exclude = 5; this.exclude = -1*int32(ones(natoms, this.max_exclude));
				else
					error("Format not supported");
				end

				this.resetmom();
      		end

		end	
	
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
				end 	
	
				fclose(fptr);
			elseif strcmp(format, "mat")		
				r = this.r; v = this.v; f = this.f; m = this.m; q = this.q; t = this.t; 
				rl = this.rl; natoms=this.natoms; lbox=this.lbox;
				
				save(filename, "r", "v", "f", "m", "q", "t", "rl", "natoms", "lbox");	
			else
				error("Format not supported");
			end
				
		end

		function tether(this, ptype, kspring)
				
			ms_tether(this.f, this.r, this.rl, ptype, this.t, kspring, this.lbox, this.natoms);
			
		end

		function mvlattice(this, ptype, dr)
			
			ms_mvlattice(this.rl, ptype, dr, this.t, this.lbox, this.natoms); 
			
		end
	
		function mom = getmom(this)
			
			for k=1:3; mom(k) = sum(this.v(:,k).*this.m); end
		
		end
	
		function resetmom(this)

			mom = this.getmom()./sum(this.m);
			for k=1:3; this.v(:,k) = this.v(:,k) - mom(k); end

		end

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
		end
		
		function vol = volume(this)
			vol = this.lbox(1)*this.lbox(2)*this.lbox(3);
		end

		function lbox = getlbox(this)
			lbox = this.lbox;
		end

		function dist = getdist(this, i, j)
			
			dr = this.r(i,:)-this.r(j,:);
			dr = wrap(dr, this.lbox);

			dist = sqrt(dot(dr,dr));
		end

		function _angle = getangle(this, a, b, c)
			
			dr1 = wrap(this.r(b,:) - this.r(a,:), this.lbox);
			dr2 = wrap(this.r(c,:) - this.r(b,:), this.lbox);			

			c11 = dot(dr1, dr1); c12 = dot(dr1, dr2); c22 = dot(dr2, dr2);

			_angle = pi - acos(c12/sqrt(c11*c22)); 

		end

		function dihedral = getdihedral(this, a, b, c, d)
			dr1 = wrap(this.r(b,:) - this.r(a,:), this.lbox);
			dr2 = wrap(this.r(c,:) - this.r(b,:), this.lbox);			
			dr3 = wrap(this.r(d,:) - this.r(c,:), this.lbox);

			c11 = dot(dr1, dr1); c12 = dot(dr1, dr2); c13 = dot(dr1, dr3);
			c22 = dot(dr2, dr2); c23 = dot(dr2, dr3); c33 = dot(dr3, dr3);

			cA = c13*c22 - c12*c23;
			cB1 = c11*c22 - c12*c12;
			cB2 = c22*c33 - c23*c23;
				
			cD = sqrt(cB1*cB2); cc = cA/cD;
			
			dihedral = pi - acos(cc);	
		end
	end

end

