#include <octave/oct.h>


DEFUN_DLD(evcorr_2, args, , "usage:c = evcorr(a,b)"){
	octave_value_list retval;

	if ( args.length() != 2 ) {
		error("evcorr: Incorrect input\n");
		return retval;
	}

	ColumnVector a( args(0).vector_value() );
	ColumnVector b( args(1).vector_value() );

	const unsigned lvec = a.numel();

	ColumnVector c(lvec);

	const double *aptr = a.fortran_vec();
	const double *bptr = b.fortran_vec();
	double *cptr = c.fortran_vec();

	for ( unsigned n=0; n<lvec; n++ ){
		cptr[n] = 0.0;
		for ( unsigned nn=0; nn<lvec-n; nn++ ){
	  		cptr[n] += aptr[nn]*bptr[n+nn];
		}
	}

	for ( unsigned n=0; n<lvec; n++ ){
		double fac = 1.0/(lvec-n);
		cptr[n] *= fac;
	}

	retval.append(c);

	return retval;

}
