# 
# atoms is a class in the molsim-package. 
# It contains all the relevat properties and methods for atoms. 
# 
# User-relevant class properties
# - Matrices (natoms by 3): r (positions), v (velocities), f (forces), rl (lattice positions)
# - Vectors (natoms elements): m (masses), q (charges) 
# - String (natoms elements): t (atom types)  
#
# Examples: See package examples/ folder
#
classdef atoms < handle

	properties (Access=public)
		# Atom properties
		r, v, f;
		m, q, t; 
		# Simulation box crossing
		boxcross;
		# Initial and last positions
		r0; rl;
		# Number of atoms and simulation box length 
		natoms,	lbox;
		# Neighbourlist etc
		nblist, max_nnb, update_nblist;
		# Pair interaction exclusion list
		exclude, max_exclude;
		# For dpd
		pv, pa;
	end

	methods
	
		## Usage: part = atoms(Configuration file name) 
		##        part = atoms([Nx, Ny, Nz], [Lx, Ly, Lz], temperature)
		##
		## Returns an instant of atoms class object 
		function this = atoms(fnameOrSize, boxLengths, temperature)

			if nargin == 0 

				this; return;

			elseif nargin == 1 

				if !exist(fnameOrSize)
					error("Configuration file does not exists");
				end	
				
				format = fnameOrSize(end-2:end);
	
				if strcmp(format, "xyz")
						
					fptr = fopen(fnameOrSize, "r");
				
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

					this.r0 = [x', y', z']; this.rl = [x', y', z'];
					this.boxcross = int32(zeros(natoms, 3)); 

				elseif strcmp(format, "mat")
					load(fnameOrSize);					
					this.r = r; this.v = v; this.f = f; this.m = m; this.q = q; this.t = t; 
					this.rl = rl; this.r0 = r0; this.boxcross = int32(boxcross);
					this.natoms = natoms; this.lbox = lbox; 
				else
					error("Format not supported");
				end

			elseif nargin == 3
				
				nx = fnameOrSize(1); ny = fnameOrSize(2); nz = fnameOrSize(3);	
				Lx = boxLengths(1); Ly = boxLengths(2); Lz = boxLengths(3); 
					
				dx = Lx/nx; dy = Ly/ny; dz = Lz/nz;
				natoms = nx*ny*nz;
			
				idx = 1;
				for n=1:nz
					for m=1:ny
						for k=1:nx
							x(idx) = (k-1)*dx; y(idx) = (m-1)*dy; z(idx) = (n-1)*dz;
							this.t(idx) = 'A';
							idx++;
						end
					end
				end	
				this.t = char(this.t');
				this.r = [x', y', z'];
				this.m = ones(natoms,1); 
				this.q = zeros(natoms,1);
				this.f = zeros(natoms, 3);
				this.lbox = [Lx, Ly, Lz]; 
				this.natoms = natoms;  
				this.r0 = [x', y', z']; this.rl = [x', y', z'];
				this.setvels(temperature);
				this.boxcross = int32(zeros(natoms, 3)); 
			end
			
			# DPD init
			for n=1:natoms
				for k=1:3	
		   			this.pv(n,k) = this.v(n,k);
				    this.pa(n,k) = 0.0;		
				end
			end

			this.update_nblist = true;
			this.max_nnb = 3000; this.nblist = -1*int32(ones(natoms, this.max_nnb)); 
			this.max_exclude = 20; this.exclude = -1*int32(ones(natoms, this.max_exclude));
	
		end	
	
		## Usage: save(filename);
		##        save(filename, option);
		##
		## Saves the current configuration to file. Supported file formats are .xyz and .mat. For xyz-format 
		## option can be be "w" for writing and "a" for appending
		##
		## Example:
		## >> atoms.save("dump-conf.mat");
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
				rl = this.rl; r0 = this.r0; boxcross = this.boxcross;
				natoms = this.natoms; lbox = this.lbox; 				
				
				save(filename, "r", "v", "f", "m", "q", "t", "rl", "r0", "natoms", "lbox", "boxcross");	
			else
				error("Format not supported");
			end
				
		end

		## Usage: thether(atom type, spring constant);
		##
		## Calculates the tethering force acting on the atoms due to the lattice 
		## coordinates specified by the propery rl. 
		## The atom type and restoring spring constant must be specified.  
		##
		## Example:
		## >> atoms.tether('W', 100.0);
		function tether(this, ptype, kspring)
				
			ms_tether(this.f, this.r, this.rl, ptype, this.t, kspring, this.lbox, this.natoms);
			
		end

		## Usage: mvlattice(atom type, distance)
		## 
		## Moves the lattice coordinates of specified atom type by a distance
		##
		## Example:
		## >> vel_wall = 0.01;
		## >> dx = vel_wall*dt;
		## >> atoms.mvlattice('W', [dx, 0.0, 0.0]);
		function mvlattice(this, ptype, dr)
			
			ms_mvlattice(this.rl, ptype, dr, this.t, this.lbox, this.natoms); 
			
		end

		## Usage: momenta = getmom();
		##
		## Returns the total x,y,z-momenta.
		##
		## Example:
		## >> mom = atoms.getmom();
		## >> printf("%e %e %e\n", mom(1), mom(2), mom(3));	
		function mom = getmom(this)
			
			for k=1:3; mom(k) = sum(this.v(:,k).*this.m); end
		
		end
	
		## Usage: resetmom();
		##
		## Resets the total momentum in all direction to zero
		function resetmom(this)

			mom = this.getmom()./sum(this.m);
			for k=1:3; this.v(:,k) = this.v(:,k) - mom(k); end

		end

		## Usage: setvels(temperature);
		##
		## Sets atom velocity corresponding to a temperature
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

		## Usage: [distance dr] = getdist(i, j)
		##
		## Returns the distance and the position vector between atoms with indices i and j  
		function [dist dr] = getdist(this, i, j)
			
			dr = this.r(i,:)-this.r(j,:);
			dr = wrap(dr, this.lbox);

			dist = sqrt(dot(dr,dr));
		end

		## Usage: angle = getangle(i,j,k)
		##
		## Returns the angle defined by the positions of atoms with indices i,j,k
		function _angle = getangle(this, a, b, c)
			
			dr1 = wrap(this.r(b,:) - this.r(a,:), this.lbox);
			dr2 = wrap(this.r(c,:) - this.r(b,:), this.lbox);			

			c11 = dot(dr1, dr1); c12 = dot(dr1, dr2); c22 = dot(dr2, dr2);

			_angle = pi - acos(c12/sqrt(c11*c22)); 

		end

		## Usage: dihedral  = getdihedral(i,j,k,l)
		##
		## Returns the dihedral angle defined by the positions of atoms with indices i,j,k,l
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

		## Usage: setexclusions(exclusion array, exclusion specifier);
		##
		## Sets exclusions for pair interactions (lj and sf). Valid specifiers are "bonds",
		## "angles", and "dihedrals".
		##
		## Examples  -  From bonds
		## >> sim.setbonds("bonds.top");
		## >> atoms.setexclusions(sim.bonds.pidx, "bonds");
		## or from dihedrals
		## >> sim.setdihedrals("dihedrals.top");
		## >> atoms.setexclusions(sim.dihedrals.pidx, "dihedrals");
		function setexclusions(this, exarray, specifier)
				
			nr = rows(exarray);	
			counter = zeros(this.natoms,1);
	
			switch (specifier)

				case "bonds"
					for n=1:nr
						for m=1:2
							idx(m) = exarray(n,m);
						end
						counter(idx(1))++; this.exclude(idx(1), counter(idx(1))) = idx(2);					

						counter(idx(2))++; this.exclude(idx(2), counter(idx(2))) = idx(1);
					end
	
				case "angles"
					for n=1:nr
						for m=1:3
							idx(m) = exarray(n,m);
						end
						
						counter(idx(1))++; 	this.exclude(idx(1), counter(idx(1)))=idx(2);		
						counter(idx(1))++; 	this.exclude(idx(1), counter(idx(1)))=idx(3);		

						counter(idx(2))++; 	this.exclude(idx(2), counter(idx(2)))=idx(1);		
						counter(idx(2))++; 	this.exclude(idx(2), counter(idx(2)))=idx(3);		
						
						counter(idx(3))++; 	this.exclude(idx(3), counter(idx(3)))=idx(1);		
						counter(idx(3))++; 	this.exclude(idx(3), counter(idx(3)))=idx(2);		
					end

				case "dihedrals"
					for n=1:nr
						for m=1:4
							idx(m) = exarray(n,m);
						end
						
						counter(idx(1))++; 	this.exclude(idx(1), counter(idx(1)))=idx(2);		
						counter(idx(1))++; 	this.exclude(idx(1), counter(idx(1)))=idx(3);		
						counter(idx(1))++; 	this.exclude(idx(1), counter(idx(1)))=idx(4);		

						counter(idx(2))++; 	this.exclude(idx(2), counter(idx(2)))=idx(1);		
						counter(idx(2))++; 	this.exclude(idx(2), counter(idx(2)))=idx(3);		
						counter(idx(2))++; 	this.exclude(idx(2), counter(idx(2)))=idx(4);	

						counter(idx(3))++; 	this.exclude(idx(3), counter(idx(3)))=idx(1);		
						counter(idx(3))++; 	this.exclude(idx(3), counter(idx(3)))=idx(2);
						counter(idx(3))++; 	this.exclude(idx(3), counter(idx(3)))=idx(4);
		
						counter(idx(4))++; 	this.exclude(idx(4), counter(idx(4)))=idx(1);		
						counter(idx(4))++; 	this.exclude(idx(4), counter(idx(4)))=idx(2);
						counter(idx(4))++; 	this.exclude(idx(4), counter(idx(4)))=idx(3);
					end

				otherwise
					error("Not a valid exlusion specifier");	
			end
	

		end
	end

end

