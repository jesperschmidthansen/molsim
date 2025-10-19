#include <octave/oct.h>
#include <omp.h>

#define HELPTXT ("Usage: ms_setomp(numthreads) \n\n Sets the number of threads to numthreads\n")

DEFUN_DLD(ms_setomp, args, ,HELPTXT){
	octave_value_list retval;

	if ( args.length() != 1 ){
		error(HELPTXT);
		return retval;
	}

	const int nthreads = args(0).int_value();
	omp_set_num_threads(nthreads);

	retval.append(nthreads);
	return retval;
}
