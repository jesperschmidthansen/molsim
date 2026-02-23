#include <octave/oct.h>
#include "ms_misc.h"

#define HELPTXT ("")

DEFUN_DLD(ms_calcmolpos, args, ,HELPTXT){
	octave_value_list retval;

	Matrix atomPos(args(0).matrix_value());
	NDArray atomMass = args(1).array_value();
	const unsigned natoms = (unsigned)args(2).scalar_value();
	Matrix atomIdx(args(3).matrix_value());
	const unsigned nuau = (unsigned)args(4).scalar_value();
	NDArray lbox = args(5).array_value();	

	int nmols = natoms/nuau;
	
	Matrix molPos(nmols, 3);
	ColumnVector molMass(nmols);

	double dr[3];

	for ( int n=0; n<nmols; n++ ){
		double x[100], y[100], z[100]; // trying stack...
		
		int aidx = atomIdx(n, 0);
		x[0] = atomPos(aidx, 0);
		y[0] = atomPos(aidx, 1);
		z[0] = atomPos(aidx, 2);

		for ( unsigned i = 1; i<nuau; i++ ){
			for ( int k=0; k<3; k++ ){
				dr[k] = atomPos(aidx+i, k) - atomPos(aidx+i-1,k);
				_Wrap(dr[k], lbox(k));
			}
			x[i] = x[i-1] + dr[0];
			y[i] = y[i-1] + dr[1];
			z[i] = z[i-1] + dr[2];
		}	
		
		double cm[3] = {0.0}; double mass = 0.0;
		for ( unsigned i = 0; i<nuau; i++ ){
			double amass = atomMass(aidx+i);
			cm[0] += amass*x[i];
			cm[1] += amass*y[i];
			cm[2] += amass*z[i];
			mass += amass;
		}	
			
		
		for ( int k=0; k<3; k++ ){
			molPos(n, k) = cm[k]/mass;
			_Periodic0(cm[k], lbox(k));
		}
		molMass(n) = mass;
	
	}	
	
	retval.append(molPos);
	retval.append(molMass);

	return retval;
}

