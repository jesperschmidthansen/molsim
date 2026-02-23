# 
# molsim is a (wrapper) class in the molsim-package. 
# It contains all the relevat properties for a simulation. 
# 
# class properties
# - sub-classes: atoms, integrator, pairforce, bonds, angels, dihedrals 
# - vectors: lbox 
# - scalars: natoms, volume,temperature 
#
# Examples: See package examples/ folder
#
classdef molsim < handle
	
	properties (Access=public)
		# Classes
		atoms;
		
		integrator;  
		thermostat;

		pairforce;  
		bonds;
		angles;
		dihedrals;
	
		molecules;
	
		# Simulation system 	
		natoms;
		lbox; volume;
		temperature;

		# Molecule info
		nmols; atom_idxs; 

		# Misc. properties
		nthreads;
	end

	methods
		
		## Usage: sim = molsim();
		## 
		## Returns an empty instance of molsim class object
		function this = molsim()
			this;	
		end

		## Usage: setconf(configuration-file) 
		##        setconf([nx, ny, nz],[Lx, Ly, Lz], temperature);     
		##
		## Sets system configuration from file or from specified dimensions
		##
		## Example:
		## >> sim = molsim();
		## >> sim.setconf("conf.xyz");
		## or
		## >> sim = molsim();
		## >> sim.setconf([10, 12, 10], [11, 14, 12], 2.0);
		function setconf(this, fnameOrSize, boxLengths, temperature)	
	
			if nargin==2
				this.atoms = atoms(fnameOrSize); 
			elseif nargin==4
				this.atoms = atoms(fnameOrSize, boxLengths, temperature);
			else
				error("Invalid call to setconf");
			end

			this.natoms = this.atoms.natoms; 
			this.lbox = this.atoms.lbox;
			this.volume = this.lbox(1)*this.lbox(2)*this.lbox(3);
			this.nthreads = 4;

			this.integrator = integrator(); 
			this.pairforce = prforce(); 
			this.thermostat = thermostat();
			this.molecules = molecules();
		end

		## Usage: setbonds(top-file) 
		##
		## Sets the bonds between atoms from file with bond information
		## 
		## Example:
		## >> sim = molsim();
		## >> sim.setbonds("bonds.top");
		##
		## See molconf 
		function setbonds(this, fname, molinfo=false)
			
			_bonds = load(fname); 
			_nbonds = rows(_bonds);
			
			this.bonds = bonds(_nbonds);
			this.bonds.pidx = _bonds(:,2:3) + 1;
			this.bonds.btypes = _bonds(:,4);			
 
			if molinfo
				load molinfo.mat 
				this.nmols = nmols;
				this.atom_idxs = atom_idxs;
	
				this.molecules.nmols = nmols;
				this.molecules.atom_idxs = atom_idxs;
			end
 
		end

		## Usage: setangles(top-file) 
		##
		## Sets the angles between atoms from file with angle information
		## 
		## Example:
		## >> sim = molsim();
		## >> sim.setangles("angles.top");
		##
		## See molconf 
		function setangles(this, fname)
			
			_angles = load(fname); 
			_nangles = rows(_angles);
			
			this.angles = angles(_nangles);
			this.angles.pidx = _angles(:,2:4) + 1;
			this.angles.atypes = _angles(:,5);			
 
		end

		## Usage: setdihedrals(top-file) 
		##
		## Sets the dihedrals between atoms from file with dihedral information
		## 
		## Example:
		## >> sim = molsim();
		## >> sim.setdihedrals("dihedrals.top");
		##
		## See molconf 
		function setdihedrals(this, fname)
			
			_dihedrals = load(fname); 
			_ndihedrals = rows(_dihedrals);
			
			this.dihedrals = dihedrals(_ndihedrals);
			this.dihedrals.pidx = _dihedrals(:,2:5) + 1;
			this.dihedrals.dtypes = _dihedrals(:,6);			
 
		end

		## Usage: settop()
		##
		## Read topology files generated molconf 
		## 
		## Example:
		## >> system("rm bonds.top angles.top dihedrals.top"); 
		## molconf("water.xyz", "water.top", [10 10 10], 2.0);  
		## sim = molsim(); 
		## sim.setconf("start.xyz"); 
		## sim.settop();
		##
		## See also water.m in examples/
		function settop(this)

			if exist("bonds.top")
				this.setbonds("bonds.top", true); 
			end
			if exist("angles.top")
				this.setangles("angles.top"); 
			end
			if exist("dihedrals.top")
				this.setdihedrals("dihedrals.top");
			end

		end

		## Usage: scalebox(target density)
		##        scalebox(target density, directions, scale factor)
		##
		## Scales the simulation box with scale factor (if not specifed this defaults to 0.999)
		##
		## Example:
		## >> sim=molsim();
		## >> sim.scalebox(0.98, [1:2], 0.9999);
		function scalebox(this, dens0, dirs=[1:3], prefac=0.999)

			densnow = this.natoms/this.volume;

			if densnow < dens0
				this.lbox(dirs) = prefac*this.lbox(dirs);
			else 
				this.lbox(dirs) = this.lbox(dirs)./prefac;
			end	
			
			this.volume = prod(this.lbox);
			this.atoms.lbox = this.lbox;

			# Ensures neighbourlist rebuild
			this.pairforce.first_call_simulation = true;				
		end

		## Usage: setnthreads(number of threads)
		##
		## Sets the number of threads for parallel computing (by default set to 4)
		##
		## Example:
		## >> sim=molsim();
		## >> sim.setnthreads(4);
		function setnthreads(this, nthreads=4)

			this.nthreads = nthreads;
			ms_setomp(nthreads);

		end

		## Usage addatom(position, atom type, mass, charge)
		##
		## Inserts an atom at a particular position. NOTE: position 
		## should be different from any other atom's position 
		##
		## Example
		## >> sim = molsim();
		## >> sim.setconf([10,10,10],[11, 11, 23], 1.0);
		## >> sim.addatom(rand*sim.lbox, 'B', 2.0, -2.);
		function addatom(this, r, atype, mass, charge)
		
			this.natoms = this.natoms + 1;
			this.atoms.natoms = this.natoms;

			this.atoms.r = [this.atoms.r; r];
			this.atoms.rl = [this.atoms.rl; r];
			this.atoms.r0 = [this.atoms.r0; r];
			this.atoms.v = [this.atoms.v; rand(1,3)];
			this.atoms.f = [this.atoms.f; [0,0,0]];
				
			this.atoms.t = [this.atoms.t; atype(1)];
			this.atoms.m = [this.atoms.m; mass];
			this.atoms.q = [this.atoms.q; charge];

			this.atoms.nblist = [this.atoms.nblist; -1*int32(ones(1, this.atoms.max_nnb))];
			this.atoms.exclude = [this.atoms.exclude; -1*int32(ones(1, this.atoms.max_exclude))];
			this.atoms.bxcrs = [this.atoms.bxcrs; [0,0,0]];
		end


		## Usage: lennardjones(atom types, lj params)
		## 
		## Calculates the Lennard-Jones pair force between particles
		## lj parameters is a 4-vector [cut-off, epsilon, sigma, aw]
		##
		## Example
		## >> sim = molsim();
		## >> sim.setconf([10,10,10],[11, 11, 23], 1.0);
		## >> [epot, Pconf] = sim.lennardjones("AA", [2.5,1.0,1.0,0.0]);
		function [epot Ppot] = lennardjones(this, atypes, params)
			[epot, Ppot] = this.pairforce.lj(this.atoms, atypes, params);		
		end

		## Usage: sfcoulomb(cutoff)
		##
		## Calculates the Coulomb interactions using a shifted force method
		## See e.g. Hansen et al., J. Phys. Chem. B, 5738:116 (2012) 
		##
		## Example: see water.m 
		function [epot Ppot] = sfcoulomb(this, cutoff)
			[epot, Ppot] = this.pairforce.sf(this.atoms, cutoff);
		end
 
		## Usage: harmonicbond(bond type)
		##
		## Calculates the bond forces using a harmonic bond potential
		##
		## Example: see butane.m
		function epot = harmonicbond(this, btype)
			epot = this.bonds.harmonic(this.atoms, btype);
		end

		## Usage: harmonicangle(angle type)
		##
		## Calculates the angle forces using a harmonic angle potential
		##
		## Example: see butane.m
		function epot = harmonicangle(this, atype)
			epot = this.angles.harmonic(this.atoms, atype);
		end

		## Usage: ryckbell(dihed. type)
		## 
		## The Ryckaert-Bellemann dihedral model
		##
		## Example: see butane.m
		function epot = ryckbell(this, dtype)
			epot = this.dihedrals.ryckbell(this.atoms, dtype);
		end

		## Usage: leapfrog()
		##
		## Integrate forward one step using the leap-frog algorithm
		##
		## Example: see lj.m
		function [ekin Pkin] = leapfrog(this)
			[ekin, Pkin] = this.integrator.lf(this.atoms, this.pairforce);
		end

		function setthermostat(this, atype, temperature, tauQ=10.0)

			this.thermostat.settype(this.atoms, atype);
			this.thermostat.temperature = temperature;

		end 


		## Usage: nosehoover()
		##
		## Apply the Nose-Hoover termostat for NVT simulations
		## See also documentation for thermostat
		function nosehoover(this)
			this.thermostat.nosehoover(this.atoms, this.integrator.dt);
		end

			
		## Usage: [rmol masses] = getmolpos()
		##
		## Get molecular centre of masses and molecular masses
		##
		## Example:
		## >> molconf("butane.xyz", "butane.top", [5 5 20], 3.0);
		## >> sim = molsim();
		## >> sim.setconf("start.xyz"); sim.settop();
		## >> sim.getmolpos()
		function [rmols mmols] = getmolpos(this)
			nuau = columns(this.atom_idxs);
			[rmols mmols] = ms_calcmolpos_oct(this.atoms.r, this.atoms.m, nuau*this.nmols, ...
												this.atom_idxs, nuau, this.lbox); 


		end
		
		## Usage: rmol = getmolvel()
		##
		## Get molecular centre of velocity
		##
		## Example:
		## >> molconf("butane.xyz", "butane.top", [5 5 20], 3.0);
		## >> sim = molsim();
		## >> sim.setconf("start.xyz"); sim.settop();
		## >> sim.getmolvel()
		function vmols = getmolvel(this)
			nuau = columns(this.atom_idxs);
			vmols = ms_calcmolvel(this.atoms.r, this.atoms.m, nuau*this.nmols, this.atom_idxs, nuau); 
		end
		
	end

end

