#include "task.h"

// Local
// Local functions

int tasktoint(char *taskstr){
  int retval = -1;

  if ( strcmp(taskstr, "lj") == 0 )
    retval = 0;
  else if ( strcmp(taskstr, "bond") == 0 )
    retval = 1;
  else if ( strcmp(taskstr, "angle") == 0 )
    retval = 2;
  else if ( strcmp(taskstr, "torsion") == 0 )
    retval = 3;
  else if ( strcmp(taskstr, "coulomb") == 0 )
    retval = 4;
  else if ( strcmp(taskstr, "dpd") == 0 )
    retval = 5;
  
  return retval;

}

void printtask(taskmanager *ptask, int tasknr){

  
  printf("Block: %d \n", ptask[tasknr].block);

  printf("Maximum cutoff: %f\n", ptask[0].max_cutoff);
  
  if ( ptask[tasknr].type == 0 ){
    printf("Task type: %d (lj)\n", ptask[tasknr].type);

    printf("Types: %s \n", ptask[tasknr].ptypes);
    printf("Cut-off: %f \n", ptask[tasknr].cutoff);
    printf("Sigma: %f \n", ptask[tasknr].sigma);
    printf("Epsilon: %f \n", ptask[tasknr].epsilon);
  }
  else if ( ptask[tasknr].type == 1 ) {
    printf("Task type: %d (bond)\n", ptask[tasknr].type);
    
    printf("Bond type: %d \n", ptask[tasknr].btype);
    printf("Bond length: %f \n", ptask[tasknr].lbond);
    printf("Spring bond: %f \n", ptask[tasknr].kbond);
  }
  else if ( ptask[tasknr].type == 2 ) {
    printf("Task type: %d (angle)\n", ptask[tasknr].type);
    
    printf("Angle type: %d \n", ptask[tasknr].btype);
    printf("Angle: %f \n", ptask[tasknr].angle);
    printf("Angle constant: %f \n", ptask[tasknr].kangle);
  }
  else if ( ptask[tasknr].type == 3 ) {
    printf("Task type: %d (torsion)\n", ptask[tasknr].type);
    
    printf("Torsion type: %d \n", ptask[tasknr].btype);
    printf("Angle constant: ");
    for ( int n=0; n<6; n++ ) printf("%f ", ptask[tasknr].dparam[n]);
    printf("\n");
  }

  printf("------------------ooo-------------------\n");

}

void settask(taskmanager *ptask, int block, char *taskopt, ... ){
  static int taskc = 0; // Task counter
  va_list arglist;
  char* types;
  double *param;

  if ( taskc == 0 ){
    ptask[0].max_cutoff = 2.5f;
    ptask[0].pairflag = false;
  }
  
  int taskt = tasktoint(taskopt); //Tasktype

  ptask[taskc].type = taskt;
  ptask[taskc].block = block;

  va_start(arglist, taskopt);
  
  switch (taskt) {
    
  case 0: // LJ task
    types = va_arg(arglist, char *);
    
    ptask[taskc].ptypes[0] = types[0];
    ptask[taskc].ptypes[1] = types[1];
    ptask[taskc].ptypes[2] = '\0';
    
    ptask[taskc].cutoff = va_arg(arglist, double);
    ptask[taskc].sigma = va_arg(arglist, double);
    ptask[taskc].epsilon = va_arg(arglist, double);

    if ( ptask[taskc].cutoff > ptask[0].max_cutoff )
      ptask[0].max_cutoff = ptask[taskc].cutoff;

    ptask[0].pairflag = true;
    
    break;

  case 1: // Bond task

    ptask[taskc].btype = va_arg(arglist, int);
    ptask[taskc].lbond =  va_arg(arglist, double);
    ptask[taskc].kbond =  va_arg(arglist, double);
    
    break;

  case 2: // Angle

    ptask[taskc].atype = va_arg(arglist, int);
    ptask[taskc].angle =  va_arg(arglist, double);
    ptask[taskc].kangle =  va_arg(arglist, double);
    
    break;

  case 3: // Dihedral
    
    ptask[taskc].dtype = va_arg(arglist, int);
    param = va_arg(arglist, double *);

    for ( int n=0; n<6; n++ )
      ptask[taskc].dparam[n] = param[n];
    
    break;

  case 4: // Coulomb
    
    ptask[taskc].dtype = va_arg(arglist, int);

    if ( ptask[taskc].cutoff > ptask[0].max_cutoff )
      ptask[0].max_cutoff = ptask[taskc].cutoff;

    ptask[0].pairflag = true;
    
    break;
  }
  // Need dpd
  
  va_end(arglist);
  
  taskc ++;
}


void dotask(double **f, seppart *atoms, taskmanager *ptask,
	    int tid, sepsys *sys){

  
  if ( ptask[tid].type == 0 ){
    const double param[3]={ptask[tid].cutoff, ptask[tid].sigma, ptask[tid].epsilon};
    sep_omp_pairs_lj(f, atoms, ptask[tid].ptypes, param, sys);
  }
  else if ( ptask[tid].type == 1 ){
    sep_omp_bond(f, atoms, ptask[tid].btype, ptask[tid].lbond,
		 ptask[tid].kbond, sys);
  }
  else if ( ptask[tid].type == 2 ){
    sep_omp_angle(f, atoms, ptask[tid].atype, ptask[tid].angle,
		  ptask[tid].kangle, sys);
  }
  else if ( ptask[tid].type == 3 ){
    sep_omp_torsion(f, atoms, ptask[tid].dtype, ptask[tid].dparam, sys);
  }
  else if ( ptask[tid].type == 4 ){
    sep_omp_coulomb(f, atoms, ptask[tid].sf_cutoff, sys);
  }

}

void dotask2(seppart *pptr, taskmanager *ptask, int ntasks, sepsys *sys,
	     const unsigned exclopt){

  const int natoms = sys->npart;

  double **f_0 = sep_matrix(natoms, 3);
  double **f_1 = sep_matrix(natoms, 3);

  if ( ptask[0].pairflag == true ){
    const double maxcf = ptask[0].max_cutoff;
    sep_omp_make_neighb(pptr, maxcf, sys, exclopt);
  }

  
#pragma omp parallel sections
  {
#pragma omp section
    {
          
      for ( int n=0; n<ntasks; n++ )
	if ( ptask[n].block == 0 ) dotask(f_0, pptr, ptask, n, sys);

    }
#pragma omp section
     {
       
       for ( int n=0; n<ntasks; n++ )
	 if ( ptask[n].block == 1 ) dotask(f_1, pptr, ptask, n, sys);

    }
  } 

  for ( int n=0; n<sys->npart; n++ )
    for ( int k=0; k<3; k++ )  pptr[n].f[k] += f_0[n][k] + f_1[n][k];
  
  sep_free_matrix(f_0, natoms);
  sep_free_matrix(f_1, natoms);

}



void dotask3(seppart *pptr, taskmanager *ptask, int ntasks, sepsys *sys,
	     const unsigned exclopt){

  
  const int natoms = sys->npart;

  double **f_0 = sep_matrix(natoms, 3);
  double **f_1 = sep_matrix(natoms, 3);
  double **f_2 = sep_matrix(natoms, 3);

  if ( ptask[0].pairflag == true ){
    const double maxcf = ptask[0].max_cutoff;
    sep_omp_make_neighb(pptr, maxcf, sys, exclopt);
  }
    
#pragma omp parallel sections
  {
#pragma omp section
    {

      for ( int n=0; n<ntasks; n++ )
	if ( ptask[n].block == 0 ) dotask(f_0, pptr, ptask, n, sys);
	
    }
#pragma omp section
    {
            
      for ( int n=0; n<ntasks; n++ )
	if ( ptask[n].block == 1 ) dotask(f_1, pptr, ptask, n, sys);
      
    }
#pragma omp section
    {

      for ( int n=0; n<ntasks; n++ )	
	if ( ptask[n].block == 2 ) dotask(f_2, pptr, ptask, n, sys);
      
    }
  } 
  
  for ( int n=0; n<sys->npart; n++ )
    for ( int k=0; k<3; k++ )
      pptr[n].f[k] += f_0[n][k] + f_1[n][k] + f_2[n][k];

  sep_free_matrix(f_0, natoms);
  sep_free_matrix(f_1, natoms);
  sep_free_matrix(f_2, natoms);

}



void dotask4(seppart *pptr, taskmanager *ptask, int ntasks, sepsys *sys,
	     const unsigned exclopt){

  const int natoms = sys->npart;

  double **f_0 = sep_matrix(natoms, 3);
  double **f_1 = sep_matrix(natoms, 3);
  double **f_2 = sep_matrix(natoms, 3);
  double **f_3 = sep_matrix(natoms, 3);

  
  if ( ptask[0].pairflag == true ){
    const double maxcf = ptask[0].max_cutoff;
    sep_omp_make_neighb(pptr, maxcf, sys, exclopt);
  }
    
#pragma omp parallel sections
  {
#pragma omp section
    {

      for ( int n=0; n<ntasks; n++ )
	if ( ptask[n].block == 0 ) dotask(f_0, pptr, ptask, n, sys);
	      
    }
#pragma omp section
    {
      
      for ( int n=0; n<ntasks; n++ )
	if ( ptask[n].block == 1 ) dotask(f_1, pptr, ptask, n, sys);
	      
    }
   
#pragma omp section
    {
      
      for ( int n=0; n<ntasks; n++ )
	if ( ptask[n].block == 2 ) dotask(f_2, pptr, ptask, n, sys);
	      
    }
    
#pragma omp section
    {
      
      for ( int n=0; n<ntasks; n++ )
	if ( ptask[n].block == 3 ) dotask(f_3, pptr, ptask, n, sys);
	      
    }
  } 

  for ( int n=0; n<sys->npart; n++ )
    for ( int k=0; k<3; k++ )
      pptr[n].f[k] += f_0[n][k] + f_1[n][k] + f_2[n][k] + f_3[n][k];

  sep_free_matrix(f_0, natoms);
  sep_free_matrix(f_1, natoms);
  sep_free_matrix(f_2, natoms);
  sep_free_matrix(f_3, natoms);
  
}
