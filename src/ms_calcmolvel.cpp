#include <octave/oct.h>
#include "ms_misc.h"

#define HELPTXT "ms_calcmolvel is a part of the molsim package \n\
It is not inteded that it is used directly by the user\n\
but only through the package class methods\n\
For usage see the package source https://github.com/jesperschmidthansen/molsim"


DEFUN_DLD(ms_calcmolvel, args, ,HELPTXT){
	octave_value_list retval;

	Matrix atomVel(args(0).matrix_value());
	NDArray atomMass = args(1).array_value();
	const unsigned natoms = (unsigned)args(2).scalar_value();
	Matrix atomIdx(args(3).matrix_value());
	const unsigned nuau = (unsigned)args(4).scalar_value();

	int nmols = natoms/nuau;
	
	Matrix molVel(nmols, 3);
	Matrix Pkin(3, 3, 0.0);

	for ( int n=0; n<nmols; n++ ){
		
		for ( int k=0; k<3; k++ ) molVel(n,k) = 0.0;
		double mass = 0.0;

		for ( unsigned i = 0; i<nuau; i++ ){
			int iadx = atomIdx(n, i) - 1;
			mass += atomMass(iadx);
			for ( int k=0; k<3; k++ ){
				molVel(n, k) += atomVel(iadx, k)*atomMass(iadx);			
			}
		}	
		
		for ( int k=0; k<3; k++ ) {
			molVel(n,k) = molVel(n,k)/mass;
			for ( int kk=0; kk<3; kk++ ) 
				Pkin(k, kk) += mass*molVel(n,k)*molVel(n,kk);
		}

	}	
	
	retval.append(molVel);
	retval.append(Pkin);

	return retval;
}

