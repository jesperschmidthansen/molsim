/* 
* sepinit.c - This file is a part of the sep-library 
*
* Copyright (C) 2008 Jesper Schmidt Hansen 
* 
* License: GPL - see COPYING for copying conditions.
* There is ABSOLUTELY NO WARRANTY, not even for MERCHANTIBILITY or
* FITNESS FOR A PARTICULAR PURPOSE.
*
* Contact: schmidt@zigzak.net
*/

#include "sepinit.h"
#include "omp.h"

seppart *sep_init(size_t npart, size_t nneighb){
  seppart *ptr;
  size_t n, k;

  ptr = malloc(npart*sizeof(seppart));
  if (ptr == NULL)
    sep_error("%s at line %d: Couldn't allocate memory\n",
	      __func__, __LINE__);

  for (n=0; n<npart; n++){
    ptr[n].neighb = malloc(nneighb*sizeof(int));
    if (ptr[n].neighb == NULL)
      sep_error("%s at line %d: Couldn't allocate memory\n",
		__func__, __LINE__);
     
  }

  /* Setting default values for mass, type etc.*/
  for (n=0; n<npart; n++){
    ptr[n].type = 'A';
    ptr[n].m = 1.0;
    ptr[n].molindex = -1;
    ptr[n].z = 0.0;
    ptr[n].ldiff = 1.0;

    for (k=0; k<3; k++)
      ptr[n].x[k] = 0.0;
    for ( k=0; k<4; k++ )
      ptr[n].bond[k] = -1;
    for (k=0; k<3; k++)
      ptr[n].crossings[k] = 0;
    for (k=0; k<3; k++)
      ptr[n].cross_neighb[k] = 0;
  }
  
  return ptr;
}

void sep_close(seppart *ptr, size_t npart){
  size_t n;

  for (n=0; n<npart; n++)
    free(ptr[n].neighb);
  

  free(ptr);
}



seppart *sep_init_xyz(double *lbox, int *npart, const char *file,  
		      char verbose){
 
  if ( verbose == 'v' )
    fprintf(stdout, "Opening %s\n", file);
  
  
  FILE *fin = fopen(file, "r");
  if ( fin==NULL )
    sep_error("%s at line %d: Couldn't open file\n",
	      __func__, __LINE__);
  

  if ( fscanf(fin, "%d", npart) != 1 )
    sep_error("%s at line %d: Error reading xyz file\n",
	      __func__, __LINE__);

  seppart *ptr = sep_init(*npart, SEP_NUM_NEIGHB);
  
  if ( fscanf(fin, "%lf%lf%lf\n", &lbox[0], &lbox[1], &lbox[2]) != 3 )
    sep_error("%s at line %d: Error reading xyz file\n",
	      __func__, __LINE__);

  for ( int n=0; n<*npart; n++ ){
    
    if ( fscanf(fin, "%c%lf%lf%lf%lf%lf%lf%lf%lf\n",
		&(ptr[n].type), 
		&(ptr[n].x[0]), &(ptr[n].x[1]), &(ptr[n].x[2]), 
		&(ptr[n].v[0]), &(ptr[n].v[1]), &(ptr[n].v[2]),
		&(ptr[n].m), &(ptr[n].z)) != 9 ){
	sep_error("%s at line %d: Error reading xyz file\n",
		  __func__, __LINE__);
    }
  }
  
  if ( verbose == 'v' ){
    fprintf(stdout, "Closing file\n");
    fprintf(stdout, "Number of particles: %d\n", *npart);
    SEP_FLUSH;
  }

  fclose(fin);

  return ptr;
}



 
void sep_set_vel(seppart*ptr, double temp, sepsys sys){
  int n, m;
  double smom[3], scale, sekin;

  srand(time(NULL));
  
  for (m=0; m<3; m++) smom[m]=0.0;
  sekin = 0.0;
  for (n=0; n<sys.npart; n++){
    for (m=0; m<3; m++){
      ptr[n].v[m] = sep_rand() - 0.5;
      smom[m]    += ptr[n].v[m]*ptr[n].m;
      sekin      += ptr[n].v[m]*ptr[n].v[m]*ptr[n].m;
    }
  }
  
  scale = sqrt(3*sys.npart*temp/sekin);
  for (n=0; n<sys.npart; n++){
    for (m=0; m<3; m++){
      ptr[n].v[m]  = (ptr[n].v[m]-smom[m]/(ptr[n].m*sys.npart))*scale;
    }
  }

  // DPD-init appears to work
  for ( int n=0; n<sys.npart; n++ ){
    for ( int k=0; k<3; k++ ) {
      ptr[n].px[k] = ptr[n].x[k] + (sep_rand()-0.5)*0.01;
      ptr[n].pv[k] = ptr[n].v[k] + (sep_rand()-0.5)*0.01;
    }
  }


}

void sep_set_vel_seed(seppart*ptr, double temp, unsigned int seed,
                      sepsys sys){

  int n, m, ndim=3, npart=sys.npart;
  double smom[3], scale, sekin;

  srand(seed);
  
  for (m=0; m<ndim; m++) smom[m]=0.0;
  sekin = 0.0;
  for (n=0; n<npart; n++){
    for (m=0; m<ndim; m++){
      ptr[n].v[m] = sep_rand() - 0.5;
      smom[m]    += ptr[n].v[m]*ptr[n].m;
      sekin      += ptr[n].v[m]*ptr[n].v[m]*ptr[n].m;
    }
  }
  
  scale = sqrt(npart*ndim*temp/sekin);
  for (n=0; n<npart; n++){
    for (m=0; m<ndim; m++){
      ptr[n].v[m]  = (ptr[n].v[m]-smom[m]/(ptr[n].m*npart))*scale;
    }
  }

  // Check on momentum and temperature
  for (m=0; m<ndim; m++) smom[m]=0.0;
  sekin = 0.0;

  for (n=0; n<npart; n++){
    for (m=0; m<ndim; m++){
      smom[m]    += ptr[n].v[m]*ptr[n].m;
      sekin      += ptr[n].v[m]*ptr[n].v[m]*ptr[n].m;
    }
  }

  // DPD-init appears to work
  for ( int n=0; n<sys.npart; n++ ){
    for ( int k=0; k<3; k++ ) {

      // additional noise within initialization
      ptr[n].px[k] = ptr[n].x[k] + (sep_rand()-0.5)*0.01;
      ptr[n].pv[k] = ptr[n].v[k] + (sep_rand()-0.5)*0.01;

      ptr[n].px[k] = ptr[n].x[k];
      ptr[n].pv[k] = ptr[n].v[k];

      ptr[n].pa[k] = 0.0;
    }
  }
 
}

void sep_set_vel_type(seppart*ptr, char type, double temp, unsigned int seed,
		      sepsys sys){
  int n, k, ntype, npart=sys.npart, ndim=3;
  double smom[3], scale, sekin;
	
  ntype = sep_count_type(ptr, type, npart);
  if ( ntype == 0 ) return;
	
  srand(seed);
	
  for (k=0; k<ndim; k++) smom[k]=0.0;
  sekin = 0.0;
  for (n=0; n<npart; n++){
    if ( ptr[n].type == type ){
      for (k=0; k<ndim; k++){
	ptr[n].v[k] = sep_rand() - 0.5;
	smom[k]  += ptr[n].v[k]*ptr[n].m;
	sekin += ptr[n].v[k]*ptr[n].v[k]*ptr[n].m;
      }
    }
  }
	
  scale = sqrt(ntype*ndim*temp/sekin);
  for (n=0; n<npart; n++){
    if ( ptr[n].type == type ){
      for (k=0; k<ndim; k++)
	ptr[n].v[k]  = (ptr[n].v[k]-smom[k]/(ptr[n].m*ntype))*scale;
	
    }
  }


  // DPD-init -  appears to work 
  for ( int n=0; n<sys.npart; n++ ){
    if ( ptr[n].type == type ){
      for ( int k=0; k<3; k++ ) {
	ptr[n].px[k] = ptr[n].x[k] + (sep_rand()-0.5)*0.01;
	ptr[n].pv[k] = ptr[n].v[k] + (sep_rand()-0.5)*0.01;
      }
    }
  }
	
}


sepsys sep_sys_setup(double lengthx, double lengthy, double lengthz, 
		     double cf, double dt, size_t npart, size_t update){
  sep3D sys;
  const double skinparam = 0.25;

  sys.length[0] = lengthx; 
  sys.length[1] = lengthy; 
  sys.length[2] = lengthz;  

  sys.volume = lengthx*lengthy*lengthz;
  //sys.nsubbox[0] = sep_nsubbox(cf, 0.0, lengthx);
  sys.nsubbox[0] = sep_nsubbox(cf, skinparam, lengthx);
  if ( update > 1 && sys.nsubbox[0] < 3 ) 
    sep_warning("%s at %d: Number of subboxes in x direction are less than three",
		__func__, __LINE__);

  //sys.nsubbox[1] = sep_nsubbox(cf, 0.0, lengthy);
  sys.nsubbox[1] = sep_nsubbox(cf, skinparam, lengthy);
  if ( update > 1 && sys.nsubbox[1] < 3 ) 
    sep_warning("%s at %d: Number of subboxes in y direction are less than three",
		__func__, __LINE__);

  //sys.nsubbox[2] = sep_nsubbox(cf, 0.0, lengthz);
  sys.nsubbox[2] = sep_nsubbox(cf, skinparam, lengthz);
  if ( update > 1 && sys.nsubbox[2] < 3 ) 
    sep_warning("%s at %d: Number of subboxes in z direction are less than three",
		__func__, __LINE__);

  sys.lsubbox[0] = lengthx/sys.nsubbox[0];
  sys.lsubbox[1] = lengthy/sys.nsubbox[1];
  sys.lsubbox[2] = lengthz/sys.nsubbox[2];

  sys.npart = npart; 
  sys.cf = cf; 
  sys.dt = dt;
  sys.ndof = 3*npart-3;
  sys.tnow = 0.0;

  sys.nupdate_neighb = 0;
  sys.neighb_update = update;
  sys.neighb_flag = 1;
  sys.skin = skinparam;

  sys.molptr = malloc(sizeof(sepmolinfo));
  if ( sys.molptr == NULL )
    sep_error("%s at line %d: Couldn't allocate memory", __func__, __LINE__);

  sys.molptr->flag_bonds = 0;
  sys.molptr->flag_angles = 0;
  sys.molptr->flag_dihedrals = 0;

  sys.molptr->flag_Fij = false;
  sys.omp_flag = false;
  
  omp_set_num_threads(1);

  return sys;
}


void sep_free_sys(sepsys *ptr){

  sep_free_bonds(ptr->molptr);
  sep_free_angles(ptr->molptr);
  sep_free_dihedrals(ptr->molptr);
  //free(ptr->interactionptr);

}

void sep_set_lattice(seppart *ptr, sepsys sys){
  double gab[3], c[3], lbox[3], dens;
  int numb[3], nX, nY, nZ, n, k;

  for (n=0; n<3; n++) lbox[n] = sys.length[n];
 
  dens = sys.npart/(lbox[0]*lbox[1]*lbox[2]);

  numb[0] = (int)ceil(pow(dens, 1./3.)*lbox[0]);
  numb[1] = (int)ceil(pow(dens, 1./3.)*lbox[1]); 
  numb[2] = (int)ceil(pow(dens, 1./3.)*lbox[2]); 

  gab[0] = lbox[0]/numb[0];
  gab[1] = lbox[1]/numb[1];
  gab[2] = lbox[2]/numb[2];
  
  n=0;
  for (nZ = 0; nZ < numb[2]; nZ++){
    c[2] = nZ*gab[2] + 1.0;
    // c[2] = nZ*gab[2] + 0.01;
    for (nY = 0; nY < numb[1]; nY++){
      c[1] = nY*gab[1] + 1.0;
      //c[1] = nY*gab[1] + 0.01;
      for (nX = 0; nX < numb[0]; nX++){
      c[0] = nX*gab[0] + 1.0;
      //c[0] = nX*gab[0] + 0.01;
	for (k=0; k<3; k++)
	  ptr[n].x[k]=c[k]; 
	n++;
	if (n == sys.npart)
	  return;
      }
    }
  }

}





