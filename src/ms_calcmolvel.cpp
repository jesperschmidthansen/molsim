#include <octave/oct.h>
#include "ms_misc.h"

#define HELPTXT ("")

DEFUN_DLD(ms_calcmolvel, args, ,HELPTXT){
	octave_value_list retval;

	Matrix atomVel(args(0).matrix_value());
	NDArray atomMass = args(1).array_value();
	const unsigned natoms = (unsigned)args(2).scalar_value();
	Matrix atomIdx(args(3).matrix_value());
	const unsigned nuau = (unsigned)args(4).scalar_value();

	int nmols = natoms/nuau;
	
	Matrix molVel(nmols, 3);

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
		
		for ( int k=0; k<3; k++ ) molVel(n,k) = molVel(n,k)/mass;
	}	
	
	retval.append(molVel);

	return retval;
}

