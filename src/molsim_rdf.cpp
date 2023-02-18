
#include <iostream>
#include <octave/oct.h>

#define HELP ("--- [rdf, r] = molsim_rdf(rdf pair types, particle positions, particle types, box dimensions, lvec) ")


double getDisplacement(double xi, double xj, double lbox){

	double dx = xi - xj; 
	
	if ( dx > 0.5*lbox ) dx -= lbox;  
	else if ( dx < -0.5*lbox ) 	dx += lbox;
	
	return dx;
}


DEFUN_DLD(molsim_rdf, args, , HELP){
	octave_value_list retval;
	
	if ( args.length() != 5 ){
		print_usage();
		return retval;
	}
		
	charMatrix pairtypes(args(0).char_matrix_value());
	Matrix positions(args(1).matrix_value());
	charMatrix types(args(2).char_matrix_value());
	ColumnVector lbox(args(3).vector_value());
	unsigned int lvec(args(4).scalar_value());
	
	double minlbox = 0.5*lbox(0);
	for ( int k=1; k<3; k++ ) {
		if ( 0.5*lbox(k) < minlbox ) minlbox = 0.5*lbox(k);
	}
	double dr = minlbox/lvec;
	
	ColumnVector rdf(lvec);
	for ( int n=0; n<lvec; n++ ) rdf(n)=0.0;
	
	unsigned npart = positions.rows();
	
	int ntype1 = 0; int ntype2 = 0;
	unsigned i;
	for ( i=0; i<npart-1; i++ ){
		
		if ( types(i) == pairtypes(0) ) ntype1++;
		if ( types(i) == pairtypes(1) ) ntype2++;
		
		for ( unsigned j=i+1; j<npart; j++ ){
			
			if ( (types(i)==pairtypes(0) && types(j)==pairtypes(1)) ||
					(types(i)==pairtypes(1) && types(j)==pairtypes(0)) ){
				
				double dist = pow(getDisplacement(positions(i,0), positions(j,0), lbox(0)), 2.0);
				dist = dist + pow(getDisplacement(positions(i,1), positions(j,1), lbox(1)), 2.0);
				dist = dist + pow(getDisplacement(positions(i,2), positions(j,2), lbox(2)), 2.0);
			
				dist = sqrt(dist);
				if ( dist < minlbox ){
					unsigned idx = dist/dr;
					rdf(idx) += 2.0;
				}
				
			}
		}
	}
	
	i++;
	if ( types(i) == pairtypes(0) ) ntype1++;
	if ( types(i) == pairtypes(1) ) ntype2++;
	
	unsigned ntypes = ntype1 + ntype2;
	
	double dens = ntypes/(lbox(0)*lbox(1)*lbox(2));
	const double fac = 4.0/3.0*3.14159; 
	ColumnVector r(lvec);
	for ( unsigned n=0; n<lvec; n++ ){
		double r1 = n*dr; double r2 = (n+1)*dr;
		double volshell = fac*(pow(r2,3.0)-pow(r1,3.0));
		double ideal = volshell*dens;
		
		rdf(n) = 4*rdf(n)/(ntypes*ideal);
		r(n) = (n+0.5)*dr;
	}
	
	retval.append(rdf);
	retval.append(r);
	
	return retval;
}

