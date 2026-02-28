This file describes the different resource files in the directory

xyz-files: Single molecule configuration files. Parameter descriptions below
            
top-files: Corresponding topology files defining bonds, angles, and torsion angles


1: water.xyz/water.top -----> H-O-H (2 part. types, 1 bond type, 1 angle) 
	SPC/Fw water model from Wu et al. J.Chem.Phys.124:024503 (2006) 	  
	Part. parameters: 
					sigma = 3.16 Å, epsilon/kB = 78.2 K, m = 16 g/mol, q/e=26.30 (H->10.782, O->-21.566)
	Mol. parameters: 
					bond, type 0: lbond = 1.012 Å, ks = 1060 kcal/(Å^2*mol)
				 	angle, type 0: angle = 1.97 rad., ka = 75.60 kcal/(rad^2*mol)

					 
2: toluene.xyz/toluene.top -----> (C)6-C (1 part. type, 2 bond types, 1 angle, 2 torsion angles)
	Toluene model from Hansen Mol. Sims. 47:1391 (2021)
	Part. parameters: 
			  		sigma = 3.675 Å, epsilon/kB = 60 K, m = 13.143 g/mol (zero charge)
	Mol. parameters:	
					bond, type 1: lbond = 1.4 Å, ks = 431  kcal/(Å^2*mol) [in phenyl grp]
					bond, type 0: lbond = 1.5 Å, ks = 431  kcal/(Å^2*mol) [phenyl grp-methyl grp]
					angle, type 0: angle = 2.09 rad. ka = 139  kcal/(rad^2*mol) 
					torsion, type 0: RB params (0, 15, 0, 0, 0, 0)  kcal/(rad^2*mol) [ensures flat ring structure]
					torsion, type 1: RB params (0, -15, 0, 0, 0, 0)  kcal/(rad^2*mol) [ensures in plane methyl grp]
	
					
3: butane.xyz/butane.top -----> C-C-C-C (1 part. type, 1 bond type, 1 angle type, 1 torsion angle)
	Flexible version of the Ryckaert-Bellemann butane model
	Part. parameters: 
					sigma = 3.9 Å, epsilon/kB = 72.1 K, m = 14.5 g/mol 
	Mol. parameters:
					bond, type 0: lbond = 1.58 Å, ks = 317  kcal/(Å^2*mol)
					angle, type 0: angle = 1.9 rad., ka = 124 kcal/(rad^2*mol)

