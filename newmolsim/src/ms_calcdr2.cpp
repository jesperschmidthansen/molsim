#include <octave/oct.h>
#include "ms_misc.h"

#define HELPTXT ("Usage: ms_calcdr2(r, r0, lbox) \n\nCalculates the maximum distance squared from r0 to r \n")

DEFUN_DLD(ms_calcdr2, args, ,HELPTXT){
	octave_value_list retval;

	if ( args.length() != 3 ) {
		error(HELPTXT);
		return retval;
	}

	Matrix r(args(0).matrix_value());
	Matrix r0(args(1).matrix_value());
	RowVector lbox(args(2).vector_value());

	const unsigned natoms = r.rows();

	double maxdr2 = 0.0f;
	for ( unsigned n=0; n<natoms; n++ ){
		double dr2 = 0.0f; 
		for ( unsigned k=0; k<3; k++ ){
			double dr = r(n,k)-r0(n,k);	
			_Wrap( dr, lbox(k) );
			dr2 += dr*dr;
		}
		if ( dr2 > maxdr2 ) maxdr2 = dr2;
	}	

	retval.append(maxdr2);

	return retval;	
}
