
classdef molsim < handle
	
	properties
		# Classes
		atoms;
		
		integrator;
		thermostat;

		pairforce;
		bonds;
		angles;
		dihedrals;
		
		# Simulation system 	
		natoms;
		lbox; volume;
		temperature;

		# Misc. properties
		nthreads;
	end

	methods
		
		function this = molsim()
			this;	
		end

		function setconf(this, fnameOrSize, boxLengths, temperature)	
	
			if nargin==2
				this.atoms = atoms(fnameOrSize); 
			elseif nargin==4
				this.atoms = atoms(fnameOrSize, boxLengths, temperature);
			end

			this.natoms = this.atoms.natoms; 
			this.lbox = this.atoms.lbox;
			this.volume = this.lbox(1)*this.lbox(2)*this.lbox(3);
			this.nthreads = 4;

			this.integrator = integrator();
			this.pairforce = prforce(); 
			
			this.thermostat = thermostat(this.atoms);
			if nargin == 3
				this.thermostat.temperature = temperature;
			end

		end

		function setbonds(this, fname)
			
			_bonds = load(fname); 
			_nbonds = rows(_bonds);
			
			this.bonds = bonds(_nbonds);
			this.bonds.pidx = _bonds(:,2:3) + 1;
			this.bonds.btypes = _bonds(:,4);			
 
		end

		function setangles(this, fname)
			
			_angles = load(fname); 
			_nangles = rows(_angles);
			
			this.angles = angles(_nangles);
			this.angles.pidx = _angles(:,2:4) + 1;
			this.angles.atypes = _angles(:,5);			
 
		end

		function setdihedrals(this, fname)
			
			_dihedrals = load(fname); 
			_ndihedrals = rows(_dihedrals);
			
			this.dihedrals = dihedrals(_ndihedrals);
			this.dihedrals.pidx = _dihedrals(:,2:5) + 1;
			this.dihedrals.dtypes = _dihedrals(:,6);			
 
		end

		function scalebox(this, dens0, prefac=0.999)

			densnow = this.natoms/this.volume;

			if densnow < dens0
				this.lbox = prefac*this.lbox;
			else 
				this.lbox = this.lbox./prefac;
			end	
			
			this.volume = prod(this.lbox);
			this.atoms.lbox = this.lbox;
				
		end

		function set_nthreads(this, nthreads)
			ms_setomp(nthreads);
			this.nthreads = nthreads;
		end

	end

end

