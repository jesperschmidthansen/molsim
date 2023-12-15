#include <octave/oct.h>
#include "cmolsim.h"

#define HELP "usage: cmolsim(<action>, <specifier>, <specifier args>) - see cmolsim documentation for more help"

enum {
  RESET=547, CALCFORCE=930, INTEGRATE=963,
  THERMOSTAT=1099, SAMPLE=642, ADD=297,
  GET=320, PRINT=557, SAVE=431,
  TASK=435, COMPRESS=876, CLEAR=519,
  SET=332, HELLO=532, LOAD=416,
  HASHVALUE=961, BAROSTAT=864, CONVERT=769 
};

unsigned hashfun(const std::string key);

void action_load(const octave_value_list& args);
void action_reset(void);
void action_calcforce(const octave_value_list& args);
void action_integrate(const octave_value_list& args);
void action_save(const octave_value_list& args);
void action_clear(void);
void action_get(octave_value_list& retval, const octave_value_list& args);
void action_set(const octave_value_list& args);
void action_thermostat(const octave_value_list& args);


DEFUN_DLD(cmolsim, args, , HELP){
	octave_value_list retval;
 
	const std::string action = args(0).string_value();

	switch ( hashfun(action) ){
		case LOAD:
			action_load(args);	break;
		case RESET:
			action_reset();	break;
		case CALCFORCE:
			action_calcforce(args);	break;
		case INTEGRATE:
			action_integrate(args); break;
		case GET:
			action_get(retval, args); break;	
		case SAVE:
			action_save(args); break;
		case CLEAR:
			action_clear();	break;
		case SET:
			action_set(args); break; 
		case THERMOSTAT:
			action_thermostat(args); break;
		default:
			octave_stdout << "Not a valid action\n";
	}

	
	return retval;	
}


unsigned hashfun(const std::string key){

  const size_t lvec = key.length();	
  unsigned sum_char = 0;
  for ( size_t n=0; n<lvec; n++ ) sum_char += (unsigned)key[n];
  
  return sum_char;

}

void action_load(const octave_value_list& args){
	
	const std::string specifier = args(1).string_value();
	const std::string filename = args(2).string_value();

	if ( strcmp(specifier.c_str(), "xyz")==0 ){
		load_xyz(filename.c_str());
	}
	else if ( strcmp(specifier.c_str(), "top")==0 ){
		load_top(filename.c_str());
	}
	else 
		error("Something went wrong with input for action 'load'");

}

void action_reset(void){
	reset_iteration();
}

void action_calcforce(const octave_value_list& args){

	const std::string specifier = args(1).string_value();

	if ( strcmp(specifier.c_str(), "lj")==0 && args.length()==7 ){
      const std::string types  =  args(2).string_value();
	  
      float cf = args(3).scalar_value();
      float sigma = args(4).scalar_value();
      float epsilon = args(5).scalar_value();
	  float aw = args(6).scalar_value();      
      
	  float ljparam[4]={sigma, epsilon, cf, aw};

	  force_lj(types.c_str(), ljparam);
	}
	else if ( strcmp(specifier.c_str(), "coulomb")==0 && (args.length()==4 || args.length()==5) ) {
		const std::string method = args(2).string_value();	// Not yet supported
		float cf = args(3).scalar_value();
		force_coulomb(cf);
	}
	else if ( strcmp(specifier.c_str(), "bond")==0 && args.length()==5 ){
		int type = args(2).int_value();
		float lbond = args(3).scalar_value();
		float ks = args(4).scalar_value();

		force_bond(type, lbond, ks);			
	}
	else if ( strcmp(specifier.c_str(), "angle")==0 && args.length()==5 ){
		int type = args(2).int_value();
		float angle = args(3).scalar_value();
		float kangle = args(4).scalar_value();

		force_angle(type, angle, kangle);
	}
	else if ( (strcmp(specifier.c_str(), "torsion")==0 || strcmp(specifier.c_str(), "dihedral")==0) && args.length()==4 ){
		int type = args(2).int_value();
		RowVector octparams = args(3).vector_value();
		
		if ( octparams.numel() != 6 )
			error("Torsion/dihedral potential needs 6 parameters");

		float params[6]; for ( int n=0; n<6; n++ ) params[n] = octparams(n);
			
		force_torsion(type, params);
	}
	else	
		error("Something went wrong with input to specifier 'calcforce'");

}

void action_integrate(const octave_value_list& args){

	const std::string specifier = args(1).string_value();

	if ( strcmp(specifier.c_str(), "leapfrog")==0 ){
		integrate_leapfrog();		
	}
	else 
		error("Something went wrong with input for action 'integrate'");

}

void action_get(octave_value_list& retval, const octave_value_list& args){

	const std::string specifier = args(1).string_value();

	if ( strcmp(specifier.c_str(), "energies")==0 )	{
		double energies[2];
		get_energies(energies);
		RowVector en(2); en(0)=energies[0]; en(1)=energies[1];
		retval.append(en);
	}	
	else if ( strcmp(specifier.c_str(), "pressure")==0 ){
		double press[4];
		get_pressure(press);
		RowVector pr(4);for ( int k=0; k<4; k++ ) pr(k)=press[k];
		
		retval.append(pr);	
	}
	else if ( strcmp(specifier.c_str(), "positions")==0 ){
		unsigned npart = get_npart();
		Matrix positions(npart, 3);
		double *tmp = (double*)malloc(3*npart*sizeof(double));	
		if ( tmp==NULL )
			error("Memory allocation error for action 'get', specifier 'positions'");

		get_positions(tmp);
		
		for ( unsigned n=0; n<npart; n++ )
			for ( int k=0; k<3; k++ ) positions(n, k) = tmp[n+k]; 

		free(tmp);
		retval.append(positions);
	}
	else if ( strcmp(specifier.c_str(), "velocities")==0 ){
		unsigned npart = get_npart();
		Matrix velocities(npart, 3);
		double *tmp = (double*)malloc(3*npart*sizeof(double));	
		if ( tmp==NULL )
			error("Memory allocation error for action 'get', specifier 'velocities'");

		get_velocities(tmp);
		
		for ( unsigned n=0; n<npart; n++ )
			for ( int k=0; k<3; k++ ) velocities(n, k) = tmp[n+k]; 

		free(tmp);
		retval.append(velocities);
	}
	else if ( strcmp(specifier.c_str(), "masses")==0 ){
		unsigned npart = get_npart();
		double *tmp = (double *) malloc(npart*sizeof(double));
		if ( tmp==NULL )
			error("Memory allocation error for action 'get', specifier 'masses'");

		get_masses(tmp);
		ColumnVector mass(npart);
		for ( unsigned n=0; n<npart; n++ ) mass(n) = tmp[n];

		free(tmp);
		retval.append(mass);
	}
	else if ( strcmp(specifier.c_str(), "charges")==0 ){
		unsigned npart = get_npart();
		double *tmp = (double *) malloc(npart*sizeof(double));
		if ( tmp==NULL )
			error("Memory allocation error for action 'get', specifier 'charges'");

		get_charges(tmp);
		ColumnVector charges(npart);
		for ( unsigned n=0; n<npart; n++ ) charges(n) = tmp[n];

		free(tmp);
		retval.append(charges);
	}
	else if ( strcmp(specifier.c_str(), "types")==0 ){
		unsigned npart = get_npart();
		char *tmp = (char*) malloc(npart*sizeof(char));
		if ( tmp==NULL )
			error("Memory allocation error for action 'get', specifier 'types'");

		get_types(tmp);
		ColumnVector charges(npart);
		for ( unsigned n=0; n<npart; n++ ) charges(n) = tmp[n];

		free(tmp);
		retval.append(charges);
	}
	else 
		error("Something went wrong with input for action 'get'");

}

void action_set(const octave_value_list& args){

	const std::string specifier = args(1).string_value();
	
	if ( strcmp(specifier.c_str(), "exclusion")==0 && args.length()==3 ){
		const std::string exclusion = args(2).string_value();
		if ( strcmp(exclusion.c_str(), "molecule") == 0  )
			set_exlusion_molecule(exclusion.c_str());
		else	
			error("Something went wrong with input for specifier 'exclusion'");
	}
	else if ( strcmp(specifier.c_str(), "timestep")==0 && args.length()==3 ){
		float dt = args(2).scalar_value();
		set_timestep(dt);
	}
	else if ( strcmp(specifier.c_str(), "momresetfrq")==0 ){
		int resetfreq = args(2).int_value();
		reset_momentum(resetfreq);
	}
	else if ( strcmp(specifier.c_str(), "cutoff")==0 ){
		float cf = args(2).scalar_value();
		set_maximum_cf(cf);	
	}
	else 
		error("Something went wrong with input for action 'set'");
}

void action_thermostat(const octave_value_list& args){

	const std::string specifier = args(1).string_value();

	if ( strcmp(specifier.c_str(), "nosehoover")==0 && args.length()==5 ){
		const std::string type  =  args(2).string_value();
	  	float temp0 = args(3).scalar_value();
		float thermostatmass = args(4).scalar_value();

		thermostat_nh(type.c_str() , temp0, thermostatmass);
	}
	else 
		error("Something went wrong with input for action 'thermostat'");
}

void action_save(const octave_value_list& args){

	const std::string specifier = args(1).string_value();

	save_xyz(specifier.c_str());		
}

void action_clear(void){ 

	free_memory();

}



