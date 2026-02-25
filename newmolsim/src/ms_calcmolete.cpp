#include <octave/oct.h>
#include "ms_misc.h"
#include <stdio.h>

#define HELPTXT ("")

DEFUN_DLD(ms_calcmolete, args, ,HELPTXT){
	octave_value_list retval;

	Matrix atomPos(args(0).matrix_value());
	const unsigned natoms = (unsigned)args(1).scalar_value();
	Matrix atomIdx(args(2).matrix_value());
	const unsigned nuau = (unsigned)args(3).scalar_value();
	NDArray lbox = args(4).array_value();	

	int nmols = natoms/nuau;
	
	Matrix molEte(nmols, 3);

	for ( int n=0; n<nmols; n++ ){
		int aidx = atomIdx(n, 0) - 1;
		int bidx = atomIdx(n, nuau-1) - 1;
				
		for ( int k=0; k<3; k++ ){
			molEte(n, k) = atomPos(aidx, k) - atomPos(bidx,k);
			_Wrap(molEte(n,k), lbox(k));
		}
	}	
	
	retval.append(molEte);

	return retval;
}

