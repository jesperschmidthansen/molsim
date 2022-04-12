
#include "molsim.h"


// Main mex
void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
 
  if ( nrhs == 0 ){
    mexPrintf("molsim - a wrapper for seplib. Check documentation. \n");
    return;
  }
  
  // ...and ACTION!!
  char *action = mxArrayToString(prhs[0]);
	 
  switch ( hashfun(action) ){

  case SET: action_set(nrhs, prhs); break;
    
  case LOAD: action_load(nrhs, prhs); break; 
   
  case RESET: action_reset(nrhs); break;

  case CALCFORCE: action_calcforce(nrhs, prhs); break;

  case INTEGRATE: action_integrate(nrhs, prhs); break;

  case THERMOSTATE: action_thermostate(nrhs, prhs); break;

  case BAROSTATE: action_barostate(nrhs, prhs); break;
    
  case SAVE: action_save(nrhs, prhs); break;

  case PRINT: action_print(); break;
    
  case GET: action_get(plhs, nrhs, prhs); break;

  case SAMPLE: action_sample(nrhs, prhs); break; 

  case TASK: action_task(nrhs, prhs); break;
    
  case COMPRESS: action_compress(nrhs, prhs); break;

  case CLEAR: action_clear(nrhs, prhs); break;

  case ADD:  action_add(nrhs, prhs); break;

  case HELLO:  mexPrintf("Hello. \n"); break;

  case HASHVALUE: action_hash(nrhs, prhs); break;
    
  default:
    mexPrintf("Action %s given -> ", action);
    mexErrMsgTxt("Not a valid action\n");
  
    break;
  }

#ifdef OCTAVE
  free(action);
#endif
    
}




/* 
 * Function definitions
 */


double spring_x0(double r2, char opt){
  double ret = 0.0;
  
  switch (opt){
  case 'f':
    ret = - SPRING_X0;
    break;
  case 'u':
    ret = 0.5*r2*SPRING_X0;
    break;
  }
  
  return ret;
}


void inputerror(const char* funstr){

  mexPrintf("From %s: ", funstr);
  mexErrMsgTxt("An input went wrong - exiting. Check carefully \n");    
      
}


unsigned hashfun(const char *key){

  const size_t len_key = strlen(key);

  unsigned sum_char = 0;
  for ( size_t n=0; n<len_key; n++ ) sum_char += (unsigned)key[n];
  
  return sum_char;

}

void action_reset(int nrhs){

  if ( nrhs != 1 ) inputerror(__func__);
  
  sep_reset_retval(&ret);
  sep_reset_force(atoms, &sys);
  
  if ( initmol )
    sep_reset_force_mol(&sys);
  
}

void action_set(int nrhs, const mxArray* prhs[]){

  char *specifier;
  
  if ( nrhs < 2 ) inputerror(__func__);
  
  specifier = mxArrayToString(prhs[1]);

  if ( strcmp(specifier, "timestep")==0 ){
    if ( nrhs != 3 ) inputerror(__func__);   
    dt = sys.dt = mxGetScalar(prhs[2]);
  }
  else if ( strcmp(specifier, "cutoff")==0 ){
    if ( nrhs != 3 ) inputerror(__func__);
    maxcutoff = mxGetScalar(prhs[2]);
  }
  else if ( strcmp(specifier, "temperature")==0 ){
    if ( nrhs != 3 ) inputerror(__func__);
    temperature = mxGetScalar(prhs[2]);
    sep_set_vel_seed(atoms, temperature, 42, sys);
    tempflag = true;
  }
  else if ( strcmp(specifier, "omp")==0 ){
    if ( nrhs != 3 ) inputerror(__func__);
    int nthreads = (int)mxGetScalar(prhs[2]);
    sep_set_omp(nthreads, &sys);
  }
  else if ( strcmp(specifier, "exclusion")==0 ){
    if ( nrhs != 3 ) inputerror(__func__);
    char *ex = mxArrayToString(prhs[2]);
    if ( strcmp(ex, "bonded" )==0 )
      exclusionflag = SEP_EXCL_BONDED;
    else if ( strcmp(ex, "molecule" )==0 )
      exclusionflag = SEP_EXCL_SAME_MOL;
    else if ( strcmp(ex, "all" )==0 )
      exclusionflag = SEP_ALL;
    else
      mexErrMsgTxt("Not valid option for exclude specifier\n");
#ifdef OCTAVE
    free(ex);
#endif
  }
  else if (strcmp(specifier, "lattice")==0 ){
    if ( nrhs != 4 ) inputerror(__func__);
    
    char str1[256]; 
    double *nxyz = mxGetPr(prhs[2]); 
    sprintf(str1, "-n=%d,%d,%d", (int)nxyz[0], (int)nxyz[1], (int)nxyz[2]);
    
    char str2[256];
    double *lbox = mxGetPr(prhs[3]);
    sprintf(str2, "-l=%f,%f,%f", lbox[0], lbox[1], lbox[2]);
    
    char str3[3]="-b";
    
    char* argv[4]; 
    argv[0] = NULL; argv[1] = str1; argv[2]=str2; argv[3]=str3; 
    
    sep_lattice(4, argv); 
  }
  else if (strcmp(specifier, "virtualsites")==0 ){
    if ( nrhs != 2 ) inputerror(__func__);
    sep_set_x0(atoms, natoms);
  }
  else if (strcmp(specifier, "types")==0 ) {
    if ( nrhs != 3 ) inputerror(__func__);
    char *types = mxArrayToString(prhs[2]);      
    for ( int n=0; n<natoms; n++ ) atoms[n].type =  types[n];
#ifdef OCTAVE
    free(types);
#endif
  }
  else if (strcmp(specifier, "force")==0 ) {
    if ( nrhs != 3 ) inputerror(__func__);
    double *force = mxGetPr(prhs[2]);      
    for ( int k=0; k<3; k++ )
      for ( int n=0; n<natoms; n++ )
	atoms[n].f[k] += force[k*natoms + n];
  }
  else if (strcmp(specifier, "compressionfactor")==0 ) {
    if ( nrhs != 3 ) inputerror(__func__);
    compressionfactor = mxGetScalar(prhs[2]);      
  }
  else if (strcmp(specifier, "temperaturerelax")==0 ) {
    if ( nrhs != 3 ) inputerror(__func__);
    taufactor = mxGetScalar(prhs[2]);      
  }
  else if (strcmp(specifier, "skin")==0 ) {
    if ( nrhs != 3 ) inputerror(__func__);
    sys.skin = mxGetScalar(prhs[2]);
  }
  else if ( strcmp(specifier, "charges")==0 ) {
    if ( nrhs != 3 ) inputerror(__func__);
    double *charge = mxGetPr(prhs[2]);      
    for ( int n=0; n<natoms; n++ ) atoms[n].z =  charge[n];
  }
  else if ( strcmp(specifier, "molconfig")==0 ) {

    int ntypes = (nrhs-2)/3;
    if ( ntypes != 1 ) inputerror(__func__);
    
    if ( ntypes == 1 ){
    
      char *xyzfile1 = mxArrayToString(prhs[2]);      
      char *topfile1 = mxArrayToString(prhs[3]);
      char nmol1[32]; sprintf(nmol1, "%d",(int)mxGetScalar(prhs[4]));
      char dens[32]; sprintf(dens, "%.6f", mxGetScalar(prhs[5]));
      char seed[32]; sprintf(seed, "%d", (int)mxGetScalar(prhs[6]));
      
      char* argv[6]; 
      argv[0] = NULL; argv[1] = xyzfile1; argv[2]=topfile1;
      argv[3] = nmol1; argv[4] = dens; argv[5] = seed; 
      
      sep_sfg(6, argv);
#ifdef OCTAVE
      free(xyzfile1); free(topfile1); 
#endif
    }
    else {
      mexErrMsgTxt("Molecular mixtures not supported just yet\n");
    }
  }
  else
    mexErrMsgTxt("Action 'set' -> not valid specifier\n");
#ifdef OCTAVE
  free(specifier);
#endif
  
}

void action_load(int nrhs, const mxArray **prhs){
  
  char *specifier = mxArrayToString(prhs[1]);
    
  if ( strcmp(specifier, "xyz")==0 ){
    
    if ( nrhs != 3 ) inputerror(__func__);
    
    if ( initflag ){
      mexPrintf("From 'load': One instant of a system already exists");
      mexPrintf(" - clear first. ");
      mexErrMsgTxt("Exiting \n");	
    }
    
    char *file = mxArrayToString(prhs[2]);
    
    atoms = sep_init_xyz(lbox, &natoms, file, 'v');
    
    sys = sep_sys_setup(lbox[0], lbox[1], lbox[2],
			maxcutoff, dt, natoms, SEP_LLIST_NEIGHBLIST);
    
    if ( !tempflag )
      sep_set_vel(atoms, temperature, sys);
    
    initflag = true;
    
#ifdef OCTAVE
    free(file);
#endif
  }
  else if ( strcmp(specifier, "top")==0 ){
    
    if ( nrhs != 3 ) inputerror(__func__);
    
    if ( !initflag ){
      mexErrMsgTxt("From 'load': You must call specifier 'xyz' before 'top'");
      return;
    }
    
    char *file = mxArrayToString(prhs[2]);
    sep_read_topology_file(atoms, file, &sys, 'v');
    mols = sep_init_mol(atoms, &sys);
    
    initmol = true;
  }
  else {
    mexErrMsgTxt("Action 'load' -> not valid specifier");
  }
#ifdef OCTAVE
  free(specifier);
#endif
  
}


void action_calcforce(int nrhs, const mxArray **prhs){

    if ( nrhs < 2 ) inputerror(__func__);

    char *specifier = mxArrayToString(prhs[1]);
    
    // van der Waal
    if ( strcmp("lj", specifier)==0 ){
      if ( nrhs != 7 ) inputerror(__func__);
      char *types =  mxArrayToString(prhs[2]);
      double cf = mxGetScalar(prhs[3]);
      double sigma = mxGetScalar(prhs[4]);
      double epsilon = mxGetScalar(prhs[5]);
      double postfac = mxGetScalar(prhs[6]);
      
      double param[4]={cf, epsilon, sigma, postfac};
      
      sep_force_lj(atoms, types, param, &sys, &ret, exclusionflag);
#ifdef OCTAVE
      free(types);
#endif
    }
    // Harmonic bond force 
    else if ( strcmp(specifier, "bond")==0 ){
      if ( nrhs != 5 ) inputerror(__func__);

      int type =  (int)mxGetScalar(prhs[2]);
      double length =  mxGetScalar(prhs[3]);
      double constant =  mxGetScalar(prhs[4]);

      sep_stretch_harmonic(atoms, type, length, constant, &sys, &ret);
    }
    // Angle force 
    else if ( strcmp(specifier, "angle")==0 ) {
      if ( nrhs != 5 ) inputerror(__func__);
      
      int type =  (int)mxGetScalar(prhs[2]);
      double angle =  mxGetScalar(prhs[3]);
      double constant =  mxGetScalar(prhs[4]);

      sep_angle_harmonic(atoms, type, angle, constant, &sys, &ret);
    }
    // Torsion force
    else if ( strcmp(specifier, "torsion")==0 ) {
      if ( nrhs != 4 ) inputerror(__func__);
 
      int type =  (int)mxGetScalar(prhs[2]);
      double *param =  mxGetPr(prhs[3]);
      
      sep_torsion_Ryckaert(atoms, type, param, &sys, &ret);
    }
    //Virtual lattice - for walls
    else if ( strcmp(specifier, "lattice")==0 ) {
      if ( nrhs != 4 ) inputerror(__func__);

      char *type =  mxArrayToString(prhs[2]);
      SPRING_X0 = mxGetScalar(prhs[3]);
      
      sep_force_x0(atoms, type[0], spring_x0, &sys);
#ifdef OCTAVE
      free(type);
#endif
    }
    // Coulomb forces
    else if ( strcmp(specifier, "coulomb")==0 ) {
      if ( nrhs < 4 && nrhs > 5 ) inputerror(__func__);
 
      char *algorithm =  mxArrayToString(prhs[2]);
      if ( strcmp(algorithm, "sf") == 0 ){
	double cf = mxGetScalar(prhs[3]);
	sep_coulomb_sf(atoms, cf, &sys, &ret, exclusionflag);
      }
      else if ( strcmp(algorithm, "wolf") == 0 ) {
	double cf = mxGetScalar(prhs[3]);
	double alpha = mxGetScalar(prhs[4]);

	sep_coulomb_wolf(atoms, alpha, cf, &sys, &ret, exclusionflag);
      }
#ifdef OCTAVE
      free(algorithm);
#endif
    }
    // DPD (only standard linear force)
    else if ( strcmp(specifier, "dpd")==0 ) {

      if ( nrhs != 6 ) inputerror(__func__);
      
      char *types =  mxArrayToString(prhs[2]);
      double cf = mxGetScalar(prhs[3]);
      double a = mxGetScalar(prhs[4]);
      double sigma = mxGetScalar(prhs[5]);

      sep_dpdforce_neighb(atoms, types, cf, a, temperature, sigma,
			  &sys, &ret, exclusionflag);
#ifdef OCTAVE
      free(types);
#endif
    }
    else {
      mexErrMsgTxt("Action 'calcforce' - not valid valid specifier\n");
    }

#ifdef OCTAVE
    free(specifier);
#endif

}


void action_integrate(int nrhs, const mxArray **prhs){

    char *specifier = mxArrayToString(prhs[1]);
    
    if ( strcmp(specifier, "leapfrog") == 0 ){
      if ( nrhs != 2 ) inputerror(__func__);
      sep_leapfrog(atoms, &sys, &ret);
    }
    else if ( strcmp(specifier, "dpd") == 0 ){
      if ( nrhs != 3 ) inputerror(__func__);
      double lambda = mxGetScalar(prhs[2]);
      sep_verlet_dpd(atoms, lambda, iterationNumber, &sys, &ret);
    }
    else
      mexErrMsgTxt("Action 'integrate' - not valid valid specifier\n");

    iterationNumber  ++;

#ifdef OCTAVE
    free(specifier);
#endif
}

void action_thermostate(int nrhs, const mxArray **prhs){
  
    if ( nrhs != 5 ) inputerror(__func__);
    
    char *specifier = mxArrayToString(prhs[1]);
    char *types =  mxArrayToString(prhs[2]);
    double Temp0 = mxGetScalar(prhs[3]);
    double tauQ =  mxGetScalar(prhs[4]);
    
    if ( strcmp(specifier, "relax") == 0 )
      sep_relax_temp(atoms, types[0], Temp0, tauQ, &sys);
    else if ( strcmp(specifier, "nosehoover")==0 )
      sep_nosehoover_type(atoms, types[0], Temp0, alpha, tauQ, &sys);
    else
      mexErrMsgTxt("Action 'thermostate' - not valid valid specifier\n");

#ifdef OCTAVE
    free(types); free(specifier);
#endif

}


void action_barostate(int nrhs, const mxArray **prhs){
  
    if ( nrhs != 4 ) inputerror(__func__);
        
    double Pd = mxGetScalar(prhs[2]);
    double beta =  mxGetScalar(prhs[3]);

    sep_berendsen(atoms, Pd, beta, &ret, &sys);
}


void action_save(int nrhs, const mxArray **prhs){

    if ( nrhs != 3 ) inputerror(__func__);
    
    char *types =  mxArrayToString(prhs[1]);
    char *file =  mxArrayToString(prhs[2]);
    
    sep_save_xyz(atoms, types, file, "w", sys);
#ifdef OCTAVE
    free(types);  free(file);
#endif
}

void action_print(void){
  
    double sump = sep_eval_mom(atoms, sys.npart);

    mexPrintf("\r");
    mexPrintf("Iteration no: %d  Epot: %.2f  Ekin: %.2f  Etot: %.4e ",
	      iterationNumber, ret.epot/sys.npart, ret.ekin/sys.npart,
	      (ret.epot+ret.ekin)/sys.npart);
    mexPrintf("Tkin: %.2f  Tot. momentum:  %1.3e  ", ret.ekin*2/(3*sys.npart-3), sump);

}

void action_get(mxArray **plhs, int nrhs, const mxArray **prhs){

  if ( nrhs != 2 && nrhs != 3 ) inputerror(__func__);
  
  char *specifier =  mxArrayToString(prhs[1]);
    
  // Positions
  if ( strcmp("positions", specifier)==0 ){
    plhs[0] = mxCreateDoubleMatrix(natoms, 3, mxREAL);
    
    double *pos = mxGetPr(plhs[0]);
    
    for ( int k=0; k<3; k++ )
      for ( int n=0; n<natoms; n++ )
	pos[k*natoms + n] = atoms[n].x[k];
  }
  // Velocities
  else if ( strcmp("velocities", specifier)==0 ){
    plhs[0] = mxCreateDoubleMatrix(natoms, 3, mxREAL);
    
    double *vel = mxGetPr(plhs[0]);
    
    for ( int k=0; k<3; k++ )
      for ( int n=0; n<natoms; n++ )
	vel[k*natoms + n] = atoms[n].v[k];
  }
  // Forces
  else if ( strcmp("forces", specifier)==0 ){
    plhs[0] = mxCreateDoubleMatrix(natoms, 3, mxREAL);
    
    double *force = mxGetPr(plhs[0]);
    
    for ( int k=0; k<3; k++ )
      for ( int n=0; n<natoms; n++ )
	  force[k*natoms + n] = atoms[n].f[k];
  }
  // Types
  else if ( strcmp("types", specifier)==0 ){
    char *types = malloc(sizeof(char)*(natoms + 1));
    if ( types == NULL )
      mexErrMsgTxt("Memory allocation error\n");
    
    for ( int n=0; n<natoms; n++ ) types[n] = atoms[n].type;
    
    plhs[0] = mxCreateString(types);
    
    free(types);
  }
  // Energies
  else if ( strcmp("energies", specifier)==0 ){
    plhs[0] = mxCreateDoubleMatrix(1, 2, mxREAL);
    
    double *energies = mxGetPr(plhs[0]);
    energies[0] = ret.ekin;
    energies[1] = ret.epot;
  }
  // Pressure
  else if ( strcmp("pressure", specifier)==0 ){
    
    sep_pressure_tensor(&ret, &sys);
    plhs[0] = mxCreateDoubleScalar(ret.p);
    
  }
  // Pressure
  else if ( strcmp("molpressure", specifier)==0 ){

    sep_eval_mol_pressure_tensor(atoms, mols, &ret, &sys);
    plhs[0] = mxCreateDoubleScalar(ret.p_mol);
    
  }
  // Number of particles 
  else if ( strcmp("numbpart", specifier)==0 ){
      
    plhs[0] = mxCreateDoubleScalar(sys.npart);
    
  }
  // Simulation box dimensions
  else if ( strcmp("box", specifier)==0 ){
    plhs[0] = mxCreateDoubleMatrix(1, 3, mxREAL);
    
    double *lbox = mxGetPr(plhs[0]);
    for ( int k=0; k<3; k++ ) lbox[k] = sys.length[k];      
  }
  // Molecular cm
  else if ( strcmp("molpositions", specifier)==0 ){
    const int nmols = sys.molptr->num_mols;
    
    plhs[0] = mxCreateDoubleMatrix(nmols, 3, mxREAL);
    double *pos = mxGetPr(plhs[0]);
    
    sep_mol_cm(atoms, mols, &sys);
    
    for ( int k=0; k<3; k++ )
      for ( int n=0; n<nmols; n++ )
	  pos[k*nmols + n] = mols[n].x[k];
  }
  // Molecular dipoles
  else if ( strcmp("moldipoles", specifier)==0 ){
    const int nmols = sys.molptr->num_mols;
	
    plhs[0] = mxCreateDoubleMatrix(nmols, 3, mxREAL);
    double *dipole = mxGetPr(plhs[0]);

    sep_mol_dipoles(atoms, mols, &sys);
    
    for ( int k=0; k<3; k++ )
      for ( int n=0; n<nmols; n++ )
	dipole[k*nmols + n] = mols[n].pel[k];
  }
  // Molecular end-to-end
  else if ( strcmp("endtoend", specifier)==0 ){
    const int nmols = sys.molptr->num_mols;
    
    plhs[0] = mxCreateDoubleMatrix(nmols, 3, mxREAL);
    double *ete = mxGetPr(plhs[0]);

    sep_mol_ete(atoms, mols, 'A', 0, mols[0].nuau-1, sys);
    
    for ( int k=0; k<3; k++ )
      for ( int n=0; n<nmols; n++ )
	ete[k*nmols + n] = mols[n].ete[k];
  }
  else if ( strcmp("numbmol", specifier)==0 ){
    const int nmols = sys.molptr->num_mols;
    plhs[0] = mxCreateDoubleScalar(nmols);
  }
  else if ( strcmp("distance", specifier)==0 ){
    double i = (int)mxGetScalar(prhs[2]);
    double j = (int)mxGetScalar(prhs[3]);
    
    plhs[1] = mxCreateDoubleMatrix(1, 3, mxREAL);
    double *r = mxGetPr(plhs[1]);
    
    double dist = sep_dist_ij(r, atoms, i, j, &sys);
    plhs[0] = mxCreateDoubleScalar(dist);
  }
  else if ( strcmp("shearpressure", specifier)==0 ){
    
    int nk =(int) mxGetScalar(prhs[2]);
    if ( nk == 0 )
      mexErrMsgTxt("Wavenumber cannot be zero\n");
    
    double k = 2*M_PI*nk/sys.length[1];
    
    sep_eval_xtrue(atoms, &sys);
    
    complex double sump = 0.0 + 0.0*I;
    complex double sump_m = 0.0 + 0.0*I;
    
    for ( int n=0; n<sys.npart; n++ ){
      double mass = atoms[n].m;
      complex double kfac = cexp(I*k*atoms[n].xtrue[1]);
      complex double kfac_m = cexp(-I*k*atoms[n].xtrue[1]);
      
      complex double a = I*atoms[n].f[0]/k - mass*atoms[n].v[0]*atoms[n].v[1];
      complex double a_m = -I*atoms[n].f[0]/k - mass*atoms[n].v[0]*atoms[n].v[1];
      
      sump += a*kfac;
      sump_m += a_m*kfac_m;
    }
    
    plhs[0] = mxCreateDoubleMatrix(1, 2, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 2, mxREAL);
    
    double *r = mxGetPr(plhs[0]);
    r[0] = creal(sump); r[1] = cimag(sump);
    
    r = mxGetPr(plhs[1]);
    r[0] = creal(sump_m); r[1] = cimag(sump_m);
    // plhs[1].real = creal(sump_m); plhs[1].imag = cimag(sump_m);
  }
  else {
    mexErrMsgTxt("Action 'get' -> Not a valid specifier\n");
  }
#ifdef OCTAVE
  free(specifier);
#endif
  
 
}

void action_sample(int nrhs, const mxArray **prhs){

      
    if ( initsampler==false ){
      initsampler = true;
      sampler = sep_init_sampler();

      if ( initmol ) 
	sep_add_mol_sampler(&sampler, mols);  
    }
    
    char *specifier =  mxArrayToString(prhs[1]);

    if ( strcmp(specifier, "do")==0 && initsampler ){
      if ( nrhs != 2 ) inputerror(__func__);
      sep_sample(atoms, &sampler, &ret, sys, iterationNumber);
    }
    else if ( strcmp(specifier, "vacf")==0 ){
      if ( nrhs != 4 ) inputerror(__func__);
      int lvec = (int)mxGetScalar(prhs[2]);
      double time = mxGetScalar(prhs[3]);
      sep_add_sampler(&sampler, "vacf", sys, lvec, time);
    }
    else if ( strcmp(specifier, "sacf")==0 ){
      if ( nrhs != 4 ) inputerror(__func__);
      int lvec = (int)mxGetScalar(prhs[2]);
      double time = mxGetScalar(prhs[3]);
      sep_add_sampler(&sampler, "sacf", sys, lvec, time);
    }
    else if ( strcmp(specifier, "hydrocorrelations")==0 ){
      if ( nrhs != 5 ) inputerror(__func__);
      int lvec = (int)mxGetScalar(prhs[2]);
      double time = mxGetScalar(prhs[3]);
      int nk = (int)mxGetScalar(prhs[4]);
      sep_add_sampler(&sampler, "gh", sys, lvec, time, nk, 0, 1);
    }
    else if ( strcmp(specifier, "profiles")==0 ){
      if ( nrhs != 5 ) inputerror(__func__);
      char *type =  mxArrayToString(prhs[2]);
      int lvec = (int)mxGetScalar(prhs[3]);
      int interval = (int)mxGetScalar(prhs[4]);
      sep_add_sampler(&sampler, "profs", sys, lvec, type[0], interval);
#ifdef OCTAVE
      free(type);
#endif
    }
    else if ( strcmp(specifier, "msacf")==0 ){
      if ( nrhs != 4 ) inputerror(__func__);
      int lvec = (int)mxGetScalar(prhs[2]);
      double time = mxGetScalar(prhs[3]);
      sep_add_sampler(&sampler, "msacf", sys, lvec, time);
    }
    else if ( strcmp(specifier, "mvacf")==0 ){
      if ( nrhs != 4 ) inputerror(__func__);
      int lvec = (int)mxGetScalar(prhs[2]);
      double time = mxGetScalar(prhs[3]);
      sep_add_sampler(&sampler, "mvacf", sys, lvec, time);
    }
    else if ( strcmp(specifier, "radial")==0 ){
      if ( nrhs != 5 ) inputerror(__func__);
      int lvec = (int)mxGetScalar(prhs[2]);
      int sampleinterval = (int)mxGetScalar(prhs[3]);
      char *types =  mxArrayToString(prhs[4]);
      sep_add_sampler(&sampler, "radial", sys, lvec, sampleinterval, types);
#ifdef OCTAVE
      free(types);
#endif
    }
    else if ( strcmp(specifier, "msd")==0 ){
      if ( nrhs != 6 ) inputerror(__func__);
      int lvec = (int)mxGetScalar(prhs[2]);
      double time = mxGetScalar(prhs[3]);
      int nwave = (int)mxGetScalar(prhs[4]);
      char *type = mxArrayToString(prhs[5]);   
      sep_add_sampler(&sampler, "msd", sys, lvec, time, nwave, type[0]);
#ifdef OCTAVE
      free(type);
#endif
    }
    else if ( strcmp(specifier, "mhydrocorrelations")==0 ){
      if ( !initmol )
	mexErrMsgTxt("Molecules are not set...\n");

      if ( nrhs != 5 ) inputerror(__func__);
    
      int lvec = (int)mxGetScalar(prhs[2]);
      double time = mxGetScalar(prhs[3]);
      int nk = (int)mxGetScalar(prhs[4]);
      sep_add_sampler(&sampler, "mgh", sys, lvec, time, nk, 0, 1);
    }
    else {
      mexErrMsgTxt("Activator 'sample' -> not valid specifier\n");
    }
#ifdef OCTAVE
    free(specifier);
#endif
}

void action_task(int nrhs, const mxArray **prhs){


  if ( inittasks == false ){
    inittasks = true;
    tasks = malloc(sizeof(taskmanager)*MAXNUMBTASK);
      if ( tasks==NULL )
	mexErrMsgTxt("Error allocating memory for task manager"); 
  }

  char *specifier =  mxArrayToString(prhs[1]);
  
  if ( strcmp("do", specifier)==0 ) {
    
    if ( nrhs != 3 ) inputerror(__func__);
      
      int numtasks = (int)mxGetScalar(prhs[2]);
      if ( numtasks==2 )
	dotask2(atoms, tasks, ntasks, &sys, exclusionflag);
      else if ( numtasks==3 )
	dotask3(atoms, tasks, ntasks, &sys, exclusionflag);
      else if ( numtasks==4 )
	dotask4(atoms, tasks, ntasks, &sys, exclusionflag);
    }
    else if ( strcmp("print", specifier)==0 ){
      if ( nrhs != 3 ) inputerror(__func__);
      int tasknr = (int)mxGetScalar(prhs[2]);
      printtask(tasks, tasknr-1);
    }
    else if ( strcmp("lj", specifier)==0 ){

      if ( nrhs != 7 ) inputerror(__func__);
    
      char *types =  mxArrayToString(prhs[2]);
      double cf = mxGetScalar(prhs[3]);
      double sigma = mxGetScalar(prhs[4]);
      double eps = mxGetScalar(prhs[5]);
      int block = (int)mxGetScalar(prhs[6]) - 1;

      settask(tasks, block, "lj", types, cf, sigma, eps);
	
      ntasks ++;
      if ( ntasks>MAXNUMBTASK )
	mexErrMsgTxt("Maximum no. of tasks exceeded \n");
#ifdef OCTAVE
      free(types);
#endif
    }
    else if ( strcmp("bond", specifier)==0 ){

      if ( nrhs != 6 ) inputerror(__func__);
    
      int type =  (int) mxGetScalar(prhs[2]);
      double bondlength = mxGetScalar(prhs[3]);
      double springconst = mxGetScalar(prhs[4]);
      int block = (int)mxGetScalar(prhs[5]) - 1;
     
      settask(tasks, block, "bond", type, bondlength, springconst);

      ntasks ++;
      if ( ntasks>MAXNUMBTASK )
	mexErrMsgTxt("Maximum no. of tasks exceeded \n");
    }
    else if ( strcmp("angle", specifier)==0 ){

      if ( nrhs != 6 ) inputerror(__func__);
    
      int type =  (int) mxGetScalar(prhs[2]);
      double angle = mxGetScalar(prhs[3]);
      double springconst = mxGetScalar(prhs[4]);
      int block = (int)mxGetScalar(prhs[5]) - 1;

      settask(tasks, block, "angle", type, angle, springconst);

      ntasks ++;
      if ( ntasks>MAXNUMBTASK )
	mexErrMsgTxt("Maximum no. of tasks exceeded \n");
    }
    else if ( strcmp("torsion", specifier)==0 ){
      if ( nrhs != 5 ) inputerror(__func__);
    
      int type =  (int) mxGetScalar(prhs[2]);
      const double *param =  mxGetPr(prhs[3]);
      int block = (int)mxGetScalar(prhs[4]) - 1;
      
      settask(tasks, block, "torsion", type, param);

      ntasks ++;
      if ( ntasks>MAXNUMBTASK )
	mexErrMsgTxt("Maximum no. of tasks exceeded \n");

    }
#ifdef OCTAVE
    free(specifier);
#endif
    
}

void action_compress(int nrhs, const mxArray **prhs){

  if ( nrhs != 2 ) inputerror(__func__);
  
  double targetdens = mxGetScalar(prhs[1]);
  sep_compress_box(atoms, targetdens, compressionfactor, &sys);  

}

void action_clear(int nrhs, const mxArray **prhs){

    if ( initflag ){
      
      if ( initsampler ) {
	initsampler = false;
	sep_close_sampler(&sampler);
      }

      if ( initmol ) {
	initmol = false;
        sep_free_mol(mols, &sys);
      }

      if ( inittasks ){
	free(tasks);
	inittasks = false;
	ntasks = 0;
      }
	   
      sep_free_sys(&sys);
      sep_close(atoms, natoms);
      
      initflag = false;
    }
    else {
      mexPrintf("No memory is allocated\n");
    }

    iterationNumber = 0;

}

void action_add(int nrhs, const mxArray **prhs){

  
    if ( nrhs < 2 ) inputerror(__func__);
    
    char *specifier =  mxArrayToString(prhs[1]);

    if ( strcmp(specifier, "force")==0 ){

      if ( nrhs != 4 ) inputerror(__func__);
    
      double *force = mxGetPr(prhs[2]);
      int dir = (int)mxGetScalar(prhs[3]);

      for ( int n=0; n<natoms; n++ )  atoms[n].f[dir] += force[n];
    }
    else if ( strcmp(specifier, "tolattice")==0 ) {

      if ( nrhs != 4 ) inputerror(__func__);
      
      double *dx = mxGetPr(prhs[2]);
      int dir = (int)mxGetScalar(prhs[3]);

      for ( int n=0; n<natoms; n++ ){
	atoms[n].x0[dir] += dx[n];
	if ( atoms[n].x0[dir] > sys.length[dir] )
	  atoms[n].x0[dir] = atoms[n].x0[dir] - sys.length[dir];
      }

    }
    else {
      mexErrMsgTxt("Not a valid specifier for action add\n");
    }
#ifdef OCTAVE
    free(specifier);
#endif

}

void action_hash(int nrhs, const mxArray **prhs){

    char *specifier = mxArrayToString(prhs[1]);    
    mexPrintf("%u \n", hashfun(specifier));

#ifdef OCTAVE
    free(specifier);
#endif

}
