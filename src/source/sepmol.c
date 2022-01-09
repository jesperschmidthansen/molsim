
#include "sepmol.h"

FILE *sep_set_file_pointer(FILE *fptr, const char *section){
  char line[256];
  
  do {

    if ( fgets(line, 256, fptr) == NULL ) break;

    if ( strcmp(line, section) == 0 ){
      if ( fgets(line, 256, fptr)==NULL )
	sep_error("%s at line %d: Read error", __func__, __LINE__);
      return fptr;
    }

  }  while ( !feof(fptr) );

  return NULL;
}

void sep_read_bonds_top(sepatom *aptr, sepmolinfo *ptr, 
			const char *file, int npart, char opt){
  const char section[] = {'[', ' ', 'b', 'o', 'n', 'd', 's', ' ', ']', 
			  '\n', '\0'};
  char line[256];
  fpos_t pos_file;
  unsigned moli, a, b, type;
 
  int *bindex = sep_vector_int(npart);

  FILE *fptr = fopen(file, "r");
  if ( fptr == NULL ) 
    sep_error("%s at line %d: Couldn't open file", __func__, __LINE__);
   
  ptr->flag_bonds = 1;
  ptr->num_bonds = 0;
  ptr->num_btypes = 0;
  ptr->num_mols = 0;

  // We *must* init the pointers since they will be freed
  // no matter if the read is sucessful or not
  ptr->blist = malloc(0);
  ptr->blengths = malloc(0);
  if ( ptr->blist == NULL ||  ptr->blengths == NULL ) 
    sep_error("%s at line %d: Couldn't allocate memory", 
	      __func__, __LINE__);

  
  // Find the 'bonds' section 
  FILE *tmp_fptr = sep_set_file_pointer(fptr, section);
  
  if ( tmp_fptr == NULL ) {
    fclose(fptr);
    return;
  }

  fptr = tmp_fptr;
  do {

    fgetpos(fptr, &pos_file); 
    if ( fgets(line, 256, fptr) == NULL  )
      sep_error("%s at line %d: Read error", __func__, __LINE__);
    
    if ( line[0] == '[' ){ 
      break;
    }
    else {
      
      fsetpos(fptr, &pos_file); 
      
      int sc = fscanf(fptr, "%u%u%u%u\n", &moli, &a, &b, &type);
      if ( sc != 4 )
	sep_error("%s at line %d: Format in top file not correct", 
		  __func__, __LINE__);
     
      (ptr->num_bonds) ++;
     
      ptr->blist = realloc(ptr->blist, sizeof(unsigned)*3*ptr->num_bonds);
      if ( ptr->blist == NULL )
	sep_error("%s at line %d: Couldn't allocate memory", 
		  __func__, __LINE__);
      
      int index0 = (ptr->num_bonds-1)*3;
      
      ptr->blist[index0] = a;
      ptr->blist[index0+1] = b;
      ptr->blist[index0+2] = type;

      // Set molecule index that atoms 'a' and 'b' are in
      aptr[a].molindex = moli;
      aptr[b].molindex = moli;

      // Set bonding infor for atoms 'a' and 'b' 
      aptr[a].bond[bindex[a]] = b; 
      aptr[b].bond[bindex[b]] = a; 
      bindex[a]++; bindex[b]++;

      if ( bindex[a] > SEP_BOND || bindex[b] > SEP_BOND )
	sep_error("%s at %d: Index exceeds the allowed number of bonds",
		    __func__, __LINE__);
	   
      if ( type > ptr->num_btypes ) ptr->num_btypes = type; 
      if ( moli > ptr->num_mols ) ptr->num_mols = moli;
    }

  } while ( !feof(fptr) ); 
  
  fclose(fptr);

  (ptr->num_btypes)++;
  (ptr->num_mols)++;

  if ( opt == 'v' ){
    printf("Succesfully read 'bond' section in file %s -> ", file);
    printf("Found %d molecules, %d bond(s) and %d bond type(s).\n", 
	   ptr->num_mols, ptr->num_bonds, ptr->num_btypes);
  }

  free(bindex);


  ptr->blengths = realloc(ptr->blengths, sizeof(double)*ptr->num_bonds);
  if ( ptr->blengths == NULL )
    sep_error("%s at %s: Memory allocation error");

}

void sep_free_bonds(sepmolinfo *ptr){

  if ( ptr->flag_bonds == 1 ){
    free(ptr->blist);
    free(ptr->blengths);
  }
  
}


void sep_read_angles_top(sepmolinfo *ptr, const char *file, char opt){
  const char section[] = {'[', ' ', 'a', 'n', 'g', 'l', 'e', 's', ' ', ']', 
			  '\n', '\0'};
  char line[256];
  fpos_t pos_file;
  unsigned moli, a, b, c, type;

  FILE *fptr = fopen(file, "r");
  if ( fptr == NULL ) 
    sep_error("%s at line %d: Couldn't open file", __func__, __LINE__);
   
  ptr->flag_angles = 1;
  ptr->num_angles = 0;
  ptr->num_atypes = 0;

  // We *must* init the pointers since they will be freed
  // no matter if the read is sucessful or not
  ptr->alist = malloc(0);
  ptr->angles = malloc(0);
  if ( ptr->alist == NULL || ptr->angles == NULL ) 
    sep_error("%s at line %d: Couldn't allocate memory", __func__, __LINE__);

  
  // Find the 'bonds' section 
  FILE *tmp_fptr = sep_set_file_pointer(fptr, section);
  
  if ( tmp_fptr == NULL ) {
    fclose(fptr);
    return;
  }

  fptr = tmp_fptr;
  do {

    fgetpos(fptr, &pos_file); 
    if ( fgets(line, 256, fptr) == NULL  )
      sep_error("%s at line %d: Read error", __func__, __LINE__);
    
    if ( line[0] == '[' ){ 
      break;
    }
    else {
      
      fsetpos(fptr, &pos_file); 
      
      int sc = fscanf(fptr, "%u%u%u%u%u\n", &moli, &a, &b, &c, &type);
      if ( sc != 5 )
	sep_error("%s at line %d: Format in top file not correct", __func__, __LINE__);
     
      (ptr->num_angles) ++;
     
      ptr->alist = realloc(ptr->alist, sizeof(unsigned)*4*ptr->num_angles);
      if ( ptr->alist == NULL )
	sep_error("%s at line %d: Couldn't allocate memory", __func__, __LINE__);
      
      int index0 = (ptr->num_angles-1)*4;
      
      ptr->alist[index0] = a;
      ptr->alist[index0+1] = b;
      ptr->alist[index0+2] = c;
      ptr->alist[index0+3] = type;

      if ( type > ptr->num_atypes ) ptr->num_atypes = type; 
      
    }

  } while ( !feof(fptr) ); 
  
  fclose(fptr);

  (ptr->num_atypes)++;
  
  if ( opt == 'v' ){
    printf("Succesfully read 'angles' section in file %s -> ", file);
    printf("Found %d angles(s) and %d bond angles(s).\n", 
	   ptr->num_angles, ptr->num_atypes);
  }

  ptr->angles = realloc(ptr->angles, sizeof(double)*ptr->num_angles);
  if ( ptr->angles == NULL )
    sep_error("%s at %s: Memory allocation error");

}

void sep_free_angles(sepmolinfo *ptr){

  if ( ptr->flag_angles == 1 ){
    free(ptr->alist);
    free(ptr->angles);
  }
}


void sep_read_dihedrals_top(sepmolinfo *ptr, const char *file, char opt){
  const char section[] = {'[', ' ', 'd', 'i', 'h', 'e', 'd', 'r', 'a', 
			  'l', 's', ' ', ']', '\n', '\0'};
  char line[256];
  fpos_t pos_file;
  unsigned moli, a, b, c, d, type;

  FILE *fptr = fopen(file, "r");
  if ( fptr == NULL ) 
    sep_error("%s at line %d: Couldn't open file", __func__, __LINE__);
 
  ptr->flag_dihedrals = 1;
  ptr->num_dihedrals = 0;
  ptr->num_dtypes = 0;
 
  // We *must* init the pointers since they will be freed
  // no matter if the read is sucessful or not
  ptr->dlist = malloc(0);
  ptr->dihedrals = malloc(0);
  if ( ptr->dlist == NULL || ptr->dihedrals == NULL ) 
    sep_error("%s at line %d: Couldn't allocate memory", 
	      __func__, __LINE__);

  // Find the 'bonds' section 
  FILE *tmp_fptr = sep_set_file_pointer(fptr, section);
  
  if ( tmp_fptr == NULL ) {
    fclose(fptr);
    return;
  }

  fptr = tmp_fptr;
  do {

    fgetpos(fptr, &pos_file); 
    if ( fgets(line, 256, fptr) == NULL  )
      sep_error("%s at line %d: Read error", __func__, __LINE__);
    
    if ( line[0] == '[' ){ 
      break;
    }
    else {
      
      fsetpos(fptr, &pos_file); 
      
      int sc = fscanf(fptr, "%u%u%u%u%u%u\n", &moli, &a, &b, &c, &d, &type);
      if ( sc != 6 )
	sep_error("%s at line %d: Format in top file not correct", __func__, __LINE__);
     
      (ptr->num_dihedrals) ++;
     
      ptr->dlist = realloc(ptr->dlist, sizeof(unsigned)*5*ptr->num_dihedrals);
      if ( ptr->dlist == NULL )
	sep_error("%s at line %d: Couldn't allocate memory", __func__, __LINE__);
      
      int index0 = (ptr->num_dihedrals-1)*5;
      
      ptr->dlist[index0] = a;
      ptr->dlist[index0+1] = b;
      ptr->dlist[index0+2] = c;
      ptr->dlist[index0+3] = d;
      ptr->dlist[index0+4] = type;

      if ( type > ptr->num_dtypes ) ptr->num_dtypes = type; 
      
    }

  } while ( !feof(fptr) ); 
  
  fclose(fptr);

  (ptr->num_dtypes)++;

  if ( opt == 'v' ){
    printf("Succesfully read 'dihedrals' section in file %s -> ", file);
    printf("Found %d dihedrals(s) and %d dihedral types(s).\n", 
	   ptr->num_dihedrals, ptr->num_dtypes);
  }

  ptr->dihedrals = realloc(ptr->dihedrals, sizeof(double)*ptr->num_dihedrals);
  if ( ptr->dihedrals == NULL )
    sep_error("%s at %s: Memory allocation error");

}

void sep_free_dihedrals(sepmolinfo *ptr){

  if ( ptr->flag_dihedrals == 1 ){
    free(ptr->dlist);
    free(ptr->dihedrals);
  }

}



void sep_read_topology_file(sepatom *aptr, const char *file, sepsys *sysptr, char opt){


  sep_read_bonds_top(aptr, sysptr->molptr, file, sysptr->npart, opt);
  sep_read_angles_top(sysptr->molptr, file, opt);
  sep_read_dihedrals_top(sysptr->molptr, file, opt);

}


void sep_stretch_harmonic(sepatom *aptr, int type, 
			  const double lbond, const double ks, 
			  sepsys *sys, sepret *ret){
  double r[3], f[3];

  unsigned num_bonds = sys->molptr->num_bonds;

  for ( unsigned n=0; n<num_bonds; n++ ){
   
    int this_type = sys->molptr->blist[3*n+2];
   
    if ( this_type == type ){
      unsigned a = sys->molptr->blist[3*n];
      unsigned b = sys->molptr->blist[3*n+1];

      double r2 = 0.0;
      for ( int k=0; k<3; k++ ){
	r[k] = aptr[a].x[k] - aptr[b].x[k];
	sep_Wrap( r[k], sys->length[k] );
	r2 += r[k]*r[k];
      }
      
      double dist = sqrt(r2);
      double ft = -ks*(dist - lbond)/dist;

      for ( int k=0; k<3; k++ ){
        f[k] = ft*r[k];
	aptr[a].f[k] += f[k];
	aptr[b].f[k] -= f[k];
      }
      
      for (int k=0; k<3; k++){
        for ( int kk=0; kk<3; kk++ ){
	  ret->pot_P[k][kk]      += f[k]*r[kk];
	  ret->pot_P_bond[k][kk] += f[k]*r[kk];
	}
      }
      ret->epot += 0.5*ks*sep_Sq(dist - lbond);
      sys->molptr->blengths[n] = dist;
    }
  }

}



void sep_angle_cossq(sepatom *ptr, int type, 
		     const double angle0, const double k, 
		     sepsys *sys, sepret *ret){
  double dr1[3], dr2[3];
  const double cCon = cos(SEP_PI - angle0);
 
  unsigned num_angles = sys->molptr->num_angles;

  for ( unsigned n=0; n<num_angles; n++ ){
   
    int this_type = sys->molptr->alist[4*n+3];
   
    if ( this_type == type ){
      unsigned a = sys->molptr->alist[4*n];
      unsigned b = sys->molptr->alist[4*n+1];
      unsigned c = sys->molptr->alist[4*n+2];

      for ( int k=0; k<3; k++ ){
	dr1[k] = ptr[b].x[k] - ptr[a].x[k];
	sep_Wrap(dr1[k], sys->length[k]);
				
	dr2[k] = ptr[c].x[k] - ptr[b].x[k];
	sep_Wrap(dr2[k], sys->length[k]);
      }
      
      double c11 = sep_dot(dr1, dr1, 3);
      double c12 = sep_dot(dr1, dr2, 3);
      double c22 = sep_dot(dr2, dr2, 3);
      
      double cD = sqrt(c11*c22); 
      double cc = c12/cD; 

      double f = -k*(cc - cCon);
      
      for ( int k=0; k<3; k++ ){
	double f1 = f*((c12/c11)*dr1[k] - dr2[k])/cD;
	double f2 = f*(dr1[k] - (c12/c22)*dr2[k])/cD;
						
	ptr[a].f[k] += f1;
	ptr[b].f[k] += (-f1-f2);
	ptr[c].f[k] += f2;
      }

      ret->epot += 0.5*k*sep_Sq(cc - cCon);
      sys->molptr->angles[n] = SEP_PI - acos(cc);
     
    }
  }

}

void sep_angle_harmonic(sepatom *ptr, int type, 
			const double angle0, const double k, 
			sepsys *sys, sepret *ret){
  double dr1[3], dr2[3];
  unsigned num_angles = sys->molptr->num_angles;

  for ( unsigned n=0; n<num_angles; n++ ){
   
    int this_type = sys->molptr->alist[4*n+3];
   
    if ( this_type == type ){
      unsigned a = sys->molptr->alist[4*n];
      unsigned b = sys->molptr->alist[4*n+1];
      unsigned c = sys->molptr->alist[4*n+2];

      for ( int k=0; k<3; k++ ){
	dr1[k] = ptr[b].x[k] - ptr[a].x[k];
	sep_Wrap(dr1[k], sys->length[k]);
				
	dr2[k] = ptr[c].x[k] - ptr[b].x[k];
	sep_Wrap(dr2[k], sys->length[k]);
      }
      
      double c11 = sep_dot(dr1, dr1, 3);
      double c12 = sep_dot(dr1, dr2, 3);
      double c22 = sep_dot(dr2, dr2, 3);
      
      double cD = sqrt(c11*c22);
      double angle = SEP_PI - acos(c12/cD); 
      
      double f = -k*(angle - angle0);
      
      for ( int k=0; k<3; k++ ){
	double f1 = f*((c12/c11)*dr1[k] - dr2[k])/cD;
	double f2 = f*(dr1[k] - (c12/c22)*dr2[k])/cD;
						
	ptr[a].f[k] += f1;
	ptr[b].f[k] += (-f1-f2);
	ptr[c].f[k] += f2;
      }

      ret->epot += 0.5*k*sep_Sq(angle - angle0);
      sys->molptr->angles[n] = angle;
     
    }
  }

}


// From Rapapport
void sep_torsion_Ryckaert(sepatom *ptr, int type, 
			  const double g[6],  sepsys *sys, 
			  sepret *ret){
  double dr1[3], dr2[3], dr3[3];
  unsigned num_dihedrals = sys->molptr->num_dihedrals;
  
  for ( unsigned n=0; n<num_dihedrals; n++ ){
   
    int this_type = sys->molptr->dlist[5*n+4];
   
    if ( this_type == type ){
      
      unsigned a = sys->molptr->dlist[5*n];
      unsigned b = sys->molptr->dlist[5*n+1];
      unsigned c = sys->molptr->dlist[5*n+2];
      unsigned d = sys->molptr->dlist[5*n+3];


      for ( int k=0; k<3; k++ ){
	dr1[k] = ptr[b].x[k] - ptr[a].x[k];
	sep_Wrap(dr1[k], sys->length[k]);
				
	dr2[k] = ptr[c].x[k] - ptr[b].x[k];
	sep_Wrap(dr2[k], sys->length[k]);
				
	dr3[k] = ptr[d].x[k] - ptr[c].x[k];
	sep_Wrap(dr3[k], sys->length[k]);
      }

      double c11 = sep_dot(dr1, dr1, 3);
      double c12 = sep_dot(dr1, dr2, 3);
      double c13 = sep_dot(dr1, dr3, 3);
      double c22 = sep_dot(dr2, dr2, 3);
      double c23 = sep_dot(dr2, dr3, 3);
      double c33 = sep_dot(dr3, dr3, 3);
      
      double cA = c13*c22 - c12*c23;
      double cB1 = c11*c22 - c12*c12;
      double cB2 = c22*c33 - c23*c23;
      double cD = sqrt(cB1*cB2); 
      double cc = cA/cD;
                  
      double f = -(g[1]+(2.*g[2]+(3.*g[3]+(4.*g[4]+5.*g[5]*cc)*cc)*cc)*cc);
      double t1 = cA; 
      double t2 = c11*c23 - c12*c13;
      double t3 = -cB1; 
      double t4 = cB2;
      double t5 = c13*c23 - c12*c33; 
      double t6 = -cA;
      double cR1 = c12/c22; 
      double cR2 = c23/c22;
      
      for ( int k=0; k<3; k++ ){
	double f1 = f*c22*(t1*dr1[k] + t2*dr2[k] + t3*dr3[k])/(cD*cB1);
	double f2 = f*c22*(t4*dr1[k] + t5*dr2[k] + t6*dr3[k])/(cD*cB2);
							
	ptr[a].f[k] += f1;
	ptr[b].f[k] += (-(1.0 + cR1)*f1 + cR2*f2);
	ptr[c].f[k] += (cR1*f1 - (1.0 + cR2)*f2);
	ptr[d].f[k] += f2;
      }

      ret->epot += g[0]+(g[1]+(g[2]+(g[3]+(g[4]+g[5]*cc)*cc)*cc)*cc)*cc;
      sys->molptr->dihedrals[n] = SEP_PI - acos(cc);
    }
  }

}


sepmol *sep_init_mol(sepatom *atom, sepsys *sys){
  int *index;
  unsigned num_mols = sys->molptr->num_mols;
	
  index = sep_vector_int(num_mols);

  sepmol *ptr = malloc(sizeof(sepmol)*num_mols);
  if ( ptr==NULL )
    sep_error("%s at %d: Couldn't allocate memory", __func__, __LINE__);

  for ( int n=0; n<sys->npart; n++ ){
    if ( atom[n].molindex > -1 )
      index[atom[n].molindex]++;
  }

  for ( unsigned n=0; n<num_mols; n++ ){
    ptr[n].nuau = index[n];
    ptr[n].shake_flag = 0;
    ptr[n].type = 'A'; //default type
    ptr[n].index = malloc(sizeof(int)*index[n]);
    if ( ptr[n].index == NULL ) 
      sep_error("%s at line %d: Couldn't allocate memory", 
		__func__, __LINE__);
  }


  sys->molptr->max_nuau = 0;
  for ( unsigned n=0; n<num_mols; n++ )
    if ( ptr[n].nuau > sys->molptr->max_nuau )
      sys->molptr->max_nuau = ptr[n].nuau;
		
  for ( unsigned n=0; n<num_mols; n++ ) index[n] = 0;

  for ( int m=0; m<sys->npart; m++ ){
    int a = atom[m].molindex;
    if (  a > -1 ){
      ptr[a].index[index[a]] = m;
      index[a]++;
    }
  }

  for ( unsigned n=0; n<num_mols; n++ ){
    ptr[n].m = 0.0;
    for ( unsigned m=0; m<ptr[n].nuau; m++ ){
      int a = ptr[n].index[m];
      ptr[n].m += atom[a].m;
    }
    for ( int k=0; k<3; k++ ) ptr[n].pel[k] = 0.0;
  }

  free(index);
  
  unsigned int a = 0;
  for ( unsigned i=1; i<(num_mols-1); i++ ) a += i;

  if ( num_mols <= SEP_MAX_NUM_MOL ){ 
    sys->molptr->flag_Fij = 1;
    sys->molptr->Fij   = sep_tensor_float(num_mols, num_mols, 3);
#ifdef COUPLE_TENSOR
    sys->molptr->Fiajb = sep_tensor_float(sys->npart, num_mols, 3);
#endif
  }
  else {
    sys->molptr->flag_Fij = 0;
    sep_warning("Molecular force/torque evaluation disabled");
  }
	
	 	
  return ptr;
}

void sep_free_mol(sepmol *ptr, sepsys *sys){

  unsigned int num_mols = sys->molptr->num_mols;
  for ( unsigned n=0; n<num_mols; n++ ){
    free(ptr[n].index);
    if ( ptr[n].shake_flag == 1 ) {  
	free(ptr[n].blength);
    }
  }

  free(ptr);
  
  if ( sys->molptr->flag_Fij == 1 ){
    sep_free_tensor_float(sys->molptr->Fij, num_mols, num_mols);
#ifdef COUPLE_TENSOR
    sep_free_tensor_float(sys->molptr->Fiajb, sys->npart, num_mols);
#endif
  }	

}


void sep_mol_cm(seppart *ptr, sepmol *mol, sepsys *sys){
  double *x, *y, *z, r[3], cm[3];
  unsigned n, m, k, i;

  unsigned nmol = sys->molptr->num_mols;

  for ( n=0; n<nmol; n++ ){

    x = sep_vector(mol[n].nuau);
    y = sep_vector(mol[n].nuau);
    z = sep_vector(mol[n].nuau);

    i = mol[n].index[0];
    
    x[0] = ptr[i].x[0];
    y[0] = ptr[i].x[1];
    z[0] = ptr[i].x[2];

    for ( m=1; m<mol[n].nuau; m++ ){
      for ( k=0; k<3; k++ ){
        r[k] = ptr[i+m].x[k] - ptr[i+m-1].x[k];
        sep_Wrap(r[k], sys->length[k]);
      }
      x[m] = x[m-1] + r[0];
      y[m] = y[m-1] + r[1];
      z[m] = z[m-1] + r[2];
    }   
      
    for ( k=0; k<3; k++ ) cm[k] = 0.0;
    
    for ( m=0; m<mol[n].nuau; m++ ){
      cm[0] += ptr[i+m].m*x[m];
      cm[1] += ptr[i+m].m*y[m];
      cm[2] += ptr[i+m].m*z[m];
    }

    for ( k=0; k<3; k++ ){
      cm[k] = cm[k]/mol[n].m;
      sep_Periodic( cm[k], sys->length[k]);
    }
      
    for ( k=0; k<3; k++ ) mol[n].x[k] = cm[k];
    
    free(x); free(y); free(z);
  }

}

void sep_mol_velcm(seppart *atom, sepmol *mol, sepsys *sys){
  unsigned i, a, n, k;
    
  for ( i=0; i<sys->molptr->num_mols; i++ ){

    sep_vector_set(&mol[i].v[0], 3, 0.0);

    for ( n=0; n<mol[i].nuau; n++ ){
      a = mol[i].index[n];
      for ( k=0; k<3; k++ )
	mol[i].v[k] += atom[a].v[k]*atom[a].m;
    } 

    for ( k=0; k<3; k++ )  mol[i].v[k] /= mol[i].m;
  }
  
}

// All this should be shiftet to coordiantes of xtrue?
void sep_mol_spin(sepatom *atom, sepmol *mol, sepsys *sys, bool safe){
  unsigned i, n, k, nmol=sys->molptr->num_mols;
  int ia;
  double dcm[3], p[3], s[3], sums[3], **inertia, mia, w[3]={0.0};
  
  sep_mol_cm(atom, mol, sys);  

  for ( i=0; i<nmol; i++ ){
    
    inertia = sep_matrix(3,3); 
    sep_vector_set(sums, 3, 0.0);
  
    for ( n=0; n<mol[i].nuau; n++ ){

      ia = mol[i].index[n];
      if ( ia == -1 ) break; // Just in case
      
      mia = atom[ia].m;
      
      for ( k=0; k<3; k++ ){
	dcm[k] = atom[ia].x[k] - mol[i].x[k];
	sep_Wrap(dcm[k], sys->length[k]);
	p[k] = atom[ia].v[k]*mia;
      } 

      sep_cross3(s, dcm, p);
      for ( k=0; k<3; k++ ) sums[k] += s[k];

      inertia[0][0] += mia*(sep_Sq(dcm[1]) + sep_Sq(dcm[2]));
      inertia[1][1] += mia*(sep_Sq(dcm[0]) + sep_Sq(dcm[2]));
      inertia[2][2] += mia*(sep_Sq(dcm[0]) + sep_Sq(dcm[1]));
      inertia[0][1] -= mia*dcm[0]*dcm[1];
      inertia[0][2] -= mia*dcm[0]*dcm[2];
      inertia[1][2] -= mia*dcm[1]*dcm[2];
    }

    // Moment of inertia is a symmetric tensor
    inertia[1][0] = inertia[0][1];
    inertia[2][0] = inertia[0][2];
    inertia[2][1] = inertia[1][2];
  
    // Angular velocity of mol i 
    if ( !safe ){
      // Check if inertia tensor is singular (can be for linear mol.)
      // This is a bit dodgy...
      double det = sep_det3(inertia);
      if ( fabs(det) <  DBL_EPSILON ){ 
	double eig[3];
	sep_eig_real_symmetric(eig, inertia);
	double Ip = (eig[0]+eig[1]+eig[2])/3.0;
	for ( k=0; k<3; k++ ) w[k] = sums[k]/Ip;
	mol[i].method_w = 0;
      }
      else {
	sep_solvelineq1(w, sums, inertia, 3);
	mol[i].method_w = 1;
      }
    }
      
    // Copy data to the structure
    for ( int k=0; k<3; k++ ){
      if ( !safe ) mol[i].w[k] = w[k];
      mol[i].s[k] = sums[k];
      for ( int kk=0; kk<3; kk++ )
	mol[i].inertia[k][kk] = inertia[k][kk];
    }

    sep_free_matrix(inertia, 3);
  }

}

void sep_save_mol_config(sepatom *atom, sepmol *mol, 
			 double time, const char *file, sepsys sys){

  unsigned nmol = sys.molptr->num_mols;

  
  sep_mol_velcm(atom, mol, &sys);
  sep_mol_spin(atom, mol, &sys, 1); //cm is evaluated here

  FILE *fout = fopen(file, "w");
  if ( fout == NULL ) 
    sep_error("%s at %s: Couldn't open file");

  fprintf(fout, "%u\n", nmol);
  fprintf(fout, "%f %f %f %f\n", sys.length[0], sys.length[1], 
	  sys.length[2], time);
	  
  for ( unsigned n=0; n<nmol; n++ )
    fprintf(fout, "%f %f %f %f %f %f %f %f %f %f %f %f\n", 
	    mol[n].x[0], mol[n].x[1], mol[n].x[2], 
	    mol[n].v[0], mol[n].v[1], mol[n].v[2],
	    mol[n].w[0], mol[n].w[1], mol[n].w[2],
	    mol[n].s[0], mol[n].s[1], mol[n].s[2]);

  fclose(fout);
 
}


void sep_FENE(sepatom *aptr, int type, const double R0, const double K,
	      sepsys *sys, sepret *ret){
  double r[3], f[3];
  const double R02 = R0*R0, fac = 0.5*K*R02;

  unsigned num_bonds = sys->molptr->num_bonds;

  for ( unsigned n=0; n<num_bonds; n++ ){
   
    int this_type = sys->molptr->blist[3*n+2];
   
    if ( this_type == type ){
      unsigned a = sys->molptr->blist[3*n];
      unsigned b = sys->molptr->blist[3*n+1];

      double r2 = 0.0;
      for ( int k=0; k<3; k++ ){
	r[k] = aptr[a].x[k] - aptr[b].x[k];
	sep_Wrap( r[k], sys->length[k] );
	r2 += r[k]*r[k];
      }
      
      double rr = r2/R02;
      double ft = - K/(1.0 - rr);
	
      for ( int k=0; k<3; k++ ){
        f[k] = ft*r[k];
	aptr[a].f[k] += f[k];
	aptr[b].f[k] -= f[k];
      }
      
      for (int k=0; k<3; k++)
	for ( int kk=0; kk<3; kk++ )
	  ret->pot_P[k][kk] += f[k]*r[kk];

      ret->epot += -fac*log(1.0 - rr);
      sys->molptr->blengths[n] = sqrt(r2);
    }
  }

}


double sep_average_bondlengths(int type, sepsys *sys){

  unsigned num_bonds = sys->molptr->num_bonds, ntype =0;
  double lbond = 0.0;

  for ( unsigned n=0; n<num_bonds; n++ ){

    int this_type = sys->molptr->blist[3*n+2];
   
    if ( this_type == type ){
      lbond += sys->molptr->blengths[n];
      ntype ++;
    }
  }

  return lbond/ntype;
}


void sep_eval_mol_pressure_tensor(sepatom *atoms, sepmol *mols, 
			     sepret *ret, sepsys *sys){
  
  if ( sys->molptr->flag_Fij == 0 ) return ;

  double rij[3];	
  int num_mol = sys->molptr->num_mols;
	
  sep_mol_cm(atoms, mols, sys);
  sep_mol_velcm(atoms, mols, sys);
	
  for ( int k=0; k<3; k++ )	
    for ( int kk=0; kk<3; kk++ )
      ret->kin_P_mol[k][kk] = ret->pot_P_mol[k][kk] = 0.0;
      
  for ( int i=0; i<num_mol; i++ ){
    for ( int k=0; k<3; k++ ){
      for ( int kk=0; kk<3; kk++ ){	
	ret->kin_P_mol[k][kk] += mols[i].m*mols[i].v[k]*mols[i].v[kk];	
      }
    }
  }
  
  for ( int i=0; i<num_mol-1; i++ ){
    for ( int j=i+1; j<num_mol; j++ ){			

      for ( int k=0; k<3; k++ ){
	rij[k] = mols[i].x[k] - mols[j].x[k];
	sep_Wrap( rij[k], sys->length[k] );
      }

      for ( int k=0; k<3; k++ )	
	for ( int kk=0; kk<3; kk++ )
	  ret->pot_P_mol[k][kk] += sys->molptr->Fij[i][j][k]*rij[kk];	
    }
  }
	
  double ivol	= 1.0/sys->volume;

  for ( int k=0; k<3; k++ )
    for ( int kk=0; kk<3; kk++ )
      ret->P_mol[k][kk] = (ret->kin_P_mol[k][kk] + 
			   ret->pot_P_mol[k][kk])*ivol; 

  ret->p_mol=0.0;
  for ( int k=0; k<3; k++ )
    ret->p_mol += ret->P_mol[k][k];

  ret->p_mol /= 3.0;
	
}


void sep_eval_mol_couple_tensor(sepatom *atoms, sepmol *mols, 
				sepret *ret, sepsys *sys){

  if ( sys->molptr->flag_Fij == 0 ) return ;

  double rij[3], Ria[3], tau[3];	
  int num_mol = sys->molptr->num_mols;
	
  sep_mol_spin(atoms, mols, sys, 0);
  sep_mol_velcm(atoms, mols, sys);
	
	
  // Kinetic part 
  for ( int i=0; i<num_mol; i++ ){
    for ( int k=0; k<3; k++ ){
      for ( int kk=0; kk<3; kk++ ){	
	ret->kin_T_mol[k][kk] += mols[i].m*mols[i].v[k]*mols[i].s[kk];	
      }
    }
  }
	
  // Potential part
  for ( int i=0; i<num_mol-1; i++ ){

    for ( int j=i+1; j<num_mol; j++ ){			
			
      for ( int k=0; k<3; k++ ) tau[k] = 0.0;
			
      for ( unsigned n=0; n<mols[i].nuau; n++ ){

	int ia = mols[i].index[n];
	// Center of mass distance vector
	for ( int k=0; k<3; k ++ ){
	  Ria[k] = atoms[ia].x[k] - mols[i].x[k];
	  sep_Wrap( Ria[k], sys->length[k] );
	}
	// Torque on i due to j	
	tau[0] += Ria[1]*sys->molptr->Fiajb[ia][j][2] - Ria[2]*sys->molptr->Fiajb[ia][j][1];
	tau[1] += Ria[2]*sys->molptr->Fiajb[ia][j][0] - Ria[0]*sys->molptr->Fiajb[ia][j][2];
	tau[2] += Ria[0]*sys->molptr->Fiajb[ia][j][1] - Ria[1]*sys->molptr->Fiajb[ia][j][0];
      }
			
      for ( int k=0; k<3; k++ ){
	rij[k] = mols[i].x[k] - mols[j].x[k];
	sep_Wrap( rij[k], sys->length[k] );
      }

      for ( int k=0; k<3; k++ )	
	for ( int kk=0; kk<3; kk++ )
	  ret->pot_T_mol[k][kk] += rij[k]*tau[kk];		
    }
  }
	
  double ivol	= 1.0/sys->volume;


  for ( int k=0; k<3; k++ )
    for ( int kk=0; kk<3; kk++ )
      ret->T_mol[k][kk] = (ret->kin_T_mol[k][kk] + 
			   ret->pot_T_mol[k][kk])*ivol; 

  ret->t_mol=0.0;
  for ( int k=0; k<3; k++ )
    ret->t_mol += ret->T_mol[k][k];
  ret->t_mol /= 3.0;
	
}


// neutral molecules only
void sep_mol_dipoles(seppart *atom, sepmol *mol, sepsys *sys){

  const int nmol = sys->molptr->num_mols;

  for ( int i=0; i<nmol; i++ ){
		
    double sumz = 0.0;
    double rpos[3]={0.0, 0.0, 0.0}, rneg[3]={0.0, 0.0, 0.0};
    int pos_counter=0, neg_counter=0;

    for ( unsigned n=0; n<mol[i].nuau; n++ ){
    
      int ia = mol[i].index[n];

      double offset[3]= {0.0};
      for ( int k=0; k<3; k++ ){
	if ( fabs(atom[ia].x[k] - mol[i].x[k]) > 0.5*sys->length[k] ){
	  if ( atom[ia].x[k] - mol[i].x[k] > 0.0 )
	    offset[k] = -sys->length[k];
	  else  
	    offset[k] = sys->length[k];
	}
      }
	 
      if ( atom[ia].z > 0.0 ){
	for ( int k=0; k<3; k++ ) rpos[k] += atom[ia].x[k] + offset[k];
	sumz += atom[ia].z;
	pos_counter++;
      }
      else if ( atom[ia].z < 0.0 ){
	for ( int k=0; k<3; k++ ) rneg[k] += atom[ia].x[k] + offset[k];
	neg_counter++;
      }

    }

    if ( neg_counter > 0 && pos_counter > 0 ){
      for ( int k=0; k<3; k++ ){
	double d = rpos[k]/pos_counter - rneg[k]/neg_counter;
	mol[i].pel[k] = sumz*d;
      }
    }

  }// nmol

}



double sep_calc_mol_temp(sepatom *atoms, sepmol *mols, sepsys sys){

  double sumpp = 0.0;

  for ( unsigned i=0; i<sys.molptr->num_mols; i++ ){
    double p[3]={0.0};
    for ( unsigned n=0; n<mols[i].nuau; n++ ){
      unsigned a = mols[i].index[n];
      for ( int k=0; k<3; k++ )
	p[k] += atoms[a].v[k]*atoms[a].m;
    }

    sumpp += (p[0]*p[0] + p[1]*p[1] + p[2]*p[2])/mols[i].m;
  }

  return sumpp/(3.0*sys.molptr->num_mols - 4.0);
}


void sep_mol_ete(sepatom *atoms, sepmol *mols, char type, 
		      unsigned a, unsigned b, sepsys sys){
  
  for ( unsigned i=0; i<sys.molptr->num_mols; i++ ){
    if ( mols[i].type == type ){
      int ia = mols[i].index[a];
      int ib = mols[i].index[b]; 
      
      for ( int k=0; k<3; k++ ){
	double dr = atoms[ia].x[k] - atoms[ib].x[k];
	sep_Wrap(dr, sys.length[k]);
	mols[i].ete[k] = dr;
      }
      
      mols[i].re2 = sep_Sq(mols[i].ete[0]) + sep_Sq(mols[i].ete[1]) + 
	sep_Sq(mols[i].ete[2]);
    }
  }
}


void sep_mol_eval_xtrue(seppart *ptr, sepmol *mol, sepsys sys){

  sep_eval_xtrue(ptr, &sys);

  unsigned nmol = sys.molptr->num_mols;

  for ( unsigned n=0; n<nmol; n++ ){

    for ( int k=0; k<3; k++ )  
      mol[n].xtrue[k] = 0.0;

    for ( unsigned m=0; m<mol[n].nuau; m++ ){
      int i = mol[n].index[m];
      for ( int k=0; k<3; k++ )  
	mol[n].xtrue[k] += ptr[i].m*ptr[i].xtrue[k];
    }

    for ( int k=0; k<3; k++ )  
	mol[n].xtrue[k] /= mol[n].m;
  }

}


int sep_count_mol_type(sepmol *mols, const char type, sepsys sys){

  int counter = 0;
  for ( size_t i=0; i<sys.molptr->num_mols; i++ )
    if ( mols[i].type == type) counter ++;

  return counter;
}



void sep_write_molconf(seppart *atoms, sepmol *mols, 
		       const char *file, sepsys sys){

  sep_mol_velcm(atoms, mols, &sys);
  sep_mol_spin(atoms, mols, &sys, 1);
  
  FILE *fout = fopen(file, "w");
  if ( fout == NULL )
    sep_error("%s at line %d: Couldn't open file", __func__, __LINE__);
  
  fprintf(fout, "%f %f %f %f\n", sys.tnow, sys.length[0], 
	  sys.length[1], sys.length[2]);

  for ( size_t i=0; i<sys.molptr->num_mols; i++ ){
    fprintf(fout, 
	    "%c %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", 
	    mols[i].type, mols[i].nuau, mols[i].m, 
	    mols[i].x[0], mols[i].x[1], mols[i].x[2],
	    mols[i].v[0], mols[i].v[1], mols[i].v[2],
	    mols[i].s[0], mols[i].s[1], mols[i].s[2],
	    mols[i].w[0], mols[i].w[1], mols[i].w[2],
	    mols[i].inertia[0][0], mols[i].inertia[1][1], mols[i].inertia[2][2],
	    mols[i].inertia[0][1], mols[i].inertia[0][2], mols[i].inertia[1][2]);
  }

  fclose(fout);
}


void sep_mol_write_config(seppart *atoms, sepmol *mols, sepsys sys){
  char file_name[256];
  static int file_counter = 0;

  sprintf(file_name, "%05d_mol.conf", file_counter);
  file_counter++;

  sep_write_molconf(atoms, mols, file_name, sys);      

}
