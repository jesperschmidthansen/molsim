#include <octave/oct.h>


DEFUN_DLD(evcorr, args, , "usage:c = evcorr(a,b)"){
  octave_value_list retval;

  if ( args.length() != 2 ) {
    error("evcorr: Incorrect input\n");
    return retval;
  }
  
  ComplexColumnVector a( args(0).complex_column_vector_value() );
  ComplexColumnVector b( args(1).complex_column_vector_value() );

  size_t lvec = a.numel();
  
  ComplexColumnVector c(a.numel());

  for ( unsigned n=0; n<lvec; n++ ){
    c(n) = 0.0;
    for ( unsigned nn=0; nn<lvec-n; nn++ ){
      c(n) += a(nn)*b(n+nn);
    }
  }

  for ( unsigned n=0; n<lvec; n++ ){
    double fac = 1.0/(lvec-n);
    c(n) *= fac;
  }
  
  retval.append(c);

  return retval;
  
}
