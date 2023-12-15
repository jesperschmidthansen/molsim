/* 
* sepprfrc.c - This file is a part of the sep-library 
*
* Copyright (C) 2008 Jesper Schmidt Hansen 
* 
* License: GPL - see COPYING for copying conditions.
* There is ABSOLUTELY NO WARRANTY, not even for MERCHANTIBILITY or
* FITNESS FOR A PARTICULAR PURPOSE.
*
* Contact: schmidt@zigzak.net
*/


#include "sepprfrc.h"
#include "sepinit.h"



void sep_force_pair_brute(seppart *ptr, const char *types, double cf,
			  double (*fun)(double, char), sepsys *sys, 
			  sepret *retval, const int opt){
  register int n,m,k;
  double r[3], r2, ft, f[3]; 

  register const double cf2 = cf*cf;

  for ( n=0; n<sys->npart-1; n++ ){
		
    int i = ptr[n].molindex;
					
    for ( m=n+1; m<sys->npart; m++ ){

      if ( opt == SEP_NEIGHB_EXCL_BONDED && sep_bonded_direct(ptr, n, m) == 1 ){
	continue;
      }
      else if ( opt == SEP_NEIGHB_EXCL_SAME_MOL && 
		ptr[n].molindex == ptr[m].molindex &&
		ptr[n].molindex != -1 )
	continue;


      if ( (ptr[n].type == types[0] && ptr[m].type == types[1]) || 
	   (ptr[n].type == types[1] && ptr[m].type == types[0]) ){
				
	r2 = 0.0;
	for ( k=0; k<3; k++ ){
	  r[k]  = ptr[n].x[k] - ptr[m].x[k];
	  sep_Wrap( r[k], sys->length[k] );
	  r2   += r[k]*r[k];
	}
	if (r2 < cf2){ 

	  // Force between particles	
	  ft = (*fun)(r2, 'f'); 
					
	  for (k=0; k<3; k++){
	    f[k] = ft*r[k];

	    ptr[n].f[k] += f[k];
	    ptr[m].f[k] -= f[k];
	  }
				
	  // Energy 
	  retval->epot += (*fun)(r2, 'u');

	  // Config. part of the stress/pressure tensor 
	  for (k=0; k<3; k++)
	    for ( int kk=0; kk<3; kk++ ) retval->pot_P[k][kk] += f[k]*r[kk];
	  
	  // Force between molecules
	  if ( sys->molptr->flag_Fij ){
	    int j = ptr[m].molindex;
	    if ( i != -1 && j != -1 ){
						
	      for ( k=0; k<3;k++ ){
		// pressure tensor	
		sys->molptr->Fij[i][j][k] += f[k];
		sys->molptr->Fij[j][i][k] -= f[k];
#ifdef COUPLE_TENSOR
		sys->molptr->Fiajb[n][j][k] += f[k];
		sys->molptr->Fiajb[m][i][k] -= f[k];
#endif
	      }
	    }	
	  }
					
	} // End of if ( r2 < rcut )
      }
    }
  }
}
	

void sep_force_pair_neighb(seppart *ptr, const char *types, double cf,
			   double (*fun)(double, char), sepsys *sys, 
			   sepret *retval, bool parallel) {
  int i1, i2, n, k, kk;
  double r2, ft, f[3], r[3];
  const double cf2 = cf*cf;
  size_t lvec = sys->npart*3;
  double *force_array = sep_vector(lvec);
 
  double epot = 0.0;
  int moli_i1, moli_i2;

  if ( parallel ){
   
    double *pconf = sep_vector(9);

#pragma omp parallel for schedule(dynamic)			\
  private(i1, n, i2, k, kk, r, r2, ft, f, moli_i1, moli_i2)	\
  reduction(+:epot, force_array[:lvec], pconf[:9]) 
    for (i1=0; i1<sys->npart; i1++){
      
      if ( ptr[i1].type != types[0] && ptr[i1].type != types[1] )
	continue;
      
      moli_i1 = ptr[i1].molindex;
      
      n = 0;
      while (1) {
	
	i2 = ptr[i1].neighb[n];
	if ( i2 == -1 ) break; 
	
	if ( (ptr[i1].type == types[0] && ptr[i2].type == types[1]) || 
	     (ptr[i1].type == types[1] && ptr[i2].type == types[0]) ){
	  r2 = 0.0;
	  for ( k=0; k<3; k++ ){
	    r[k] = ptr[i1].x[k] - ptr[i2].x[k];
	    sep_Wrap( r[k], sys->length[k] );
	    r2 += r[k]*r[k];
	  }
	  
	  if ( r2 < cf2 ){
	    
	    ft = (*fun)(r2, 'f');
	    for ( k=0; k<3; k++ ){
	      f[k] = ft*r[k];
	      
	      force_array[i1*3 + k] += f[k];
	      force_array[i2*3 + k] += -f[k];
	    }
	    epot += (*fun)(r2, 'u');
	    
	    for (k=0; k<3; k++)
	      for ( kk=0; kk<3; kk++ )
		pconf[k*3+kk] += f[k]*r[kk];
	  }
	}
	n++;
      }
    }

    for ( int k=0; k<3; k++ )
      for ( int kk=0; kk<3; kk++ ) retval->pot_P[k][kk] += pconf[3*k + kk];

    free(pconf);
  }
  
  else {
   
    for (i1=0; i1<sys->npart; i1++){
    
      if ( ptr[i1].type != types[0] && ptr[i1].type != types[1] )
	continue;
      
      moli_i1 = ptr[i1].molindex;
      n = 0;
      while (1) {
	
	i2 = ptr[i1].neighb[n];
	if ( i2 == -1 ) break; 
	
	if ( (ptr[i1].type == types[0] && ptr[i2].type == types[1]) || 
	     (ptr[i1].type == types[1] && ptr[i2].type == types[0]) ){
	  r2 = 0.0;
	  for ( k=0; k<3; k++ ){
	    r[k] = ptr[i1].x[k] - ptr[i2].x[k];
	    sep_Wrap( r[k], sys->length[k] );
	    r2 += r[k]*r[k];
	  }
	  
	  if ( r2 < cf2 ){
	    
	    ft = (*fun)(r2, 'f');
	    for ( k=0; k<3; k++ ){
	      f[k] = ft*r[k];
	      
	      force_array[i1*3 + k] += f[k];
	      force_array[i2*3 + k] += -f[k];
	    }
	    epot += (*fun)(r2, 'u'); 
	    
	    for (k=0; k<3; k++)
	      for ( kk=0; kk<3; kk++ )
		retval->pot_P[k][kk] += f[k]*r[kk];
	    
	    if ( sys->molptr->flag_Fij == 1 ){
	      moli_i2 = ptr[i2].molindex;
	      if ( (moli_i1 != -1 && moli_i2 != -1) && (moli_i1 != moli_i2) ){
		for ( k=0; k<3;k++ ){
		  sys->molptr->Fij[moli_i1][moli_i2][k] += f[k];
		  sys->molptr->Fij[moli_i2][moli_i1][k] -= f[k];	
		}
	      }	
	    }
	    
	  }
	}
	n++;
      }
    }
  }
  
  for ( int n=0; n<sys->npart; n++ )
    for ( int k=0; k<3; k++ )
      ptr[n].f[k] += force_array[n*3 + k];
  
  free(force_array);

  retval->epot = epot;

}

int sep_force_pairs(seppart *ptr, const char *types, double cf,
		     double (*fun)(double, char), sepsys *sys, 
		     sepret *retval, const unsigned opt){

	if ( cf > sys->cf ){
		#ifdef OCTAVE
			sep_warning("cutoff for an interaction cannot be larger than maximum cutoff");
			return SEP_FAILURE;
		#else
			sep_error("cutoff for an interaction cannot be larger than maximum cutoff");
		#endif
	}
    
  
  if ( sys->neighb_update == SEP_BRUTE ) {

    if (  sys->omp_flag ){
      sep_warning("omp flag set, SEP_BRUTE does not support threads.");
      sep_warning("Resetting omp flag");
      sys->omp_flag = false;
    }
    
    sep_force_pair_brute(ptr, types, cf, fun, sys, retval, opt);    
  }
  else {
    if ( sys->neighb_flag == 1 ){
      
      if ( sys->neighb_update == SEP_NEIGHBLIST ){
	sep_make_neighblist(ptr, sys, opt);
      }
      else if ( sys->neighb_update == SEP_LLIST_NEIGHBLIST ){
	
	if ( opt == SEP_ALL )
	  sep_neighb(ptr, sys);
	else if ( opt == SEP_NEIGHB_EXCL_BONDED )
	  sep_neighb_nonbonded(ptr, sys);
	else if ( opt == SEP_NEIGHB_EXCL_SAME_MOL )
	  sep_neighb_excl_same_mol(ptr, sys);
      }

      sys->neighb_flag = 0;
    }

    sep_force_pair_neighb(ptr, types, cf, fun, sys, retval, sys->omp_flag);
    
  }

  return SEP_SUCCESS;
}



void sep_force_dpd(seppart *ptr, const char *types, 
		   const double cf, const double aij, 
		   const double temp_desired, 
		   const double sigma, sepsys *sys, sepret *retval,
		   const unsigned opt){

  if ( sys->neighb_update == SEP_BRUTE ) {

    if (  sys->omp_flag ){
      sep_warning("omp flag set, SEP_BRUTE does not support threads.");
      sep_warning("Resetting omp flag");
      sys->omp_flag = false;
    }

    sep_dpdforce_brute(ptr, types, cf, aij, temp_desired, sigma, sys, retval);
    
  }
  else {
    
    sep_dpdforce_neighb(ptr, types, cf, aij, temp_desired, sigma, sys, retval, opt);
    
  }
  
}
  
  

// Making a neighb list - but not from linked-list
void sep_make_neighblist(seppart *ptr, sepsys *sys, const unsigned opt){
  int n, m, k, index;
  double r[3], r2, cf2; 
  
  cf2 = sep_Sq(sys->cf + sys->skin); 

  for ( n=0; n<sys->npart; n++ )
    for ( m=0; m<SEP_NEIGHB; m++ ) ptr[n].neighb[m] = -1;

  for ( n=0; n<sys->npart-1; n++ ){
    index = 0;
    for ( m=n+1; m<sys->npart; m++ ){
      
      if ( opt == SEP_NEIGHB_EXCL_BONDED && 
           sep_bonded_direct(ptr, n, m) == 1 ){
        continue;
      }
      else if ( opt == SEP_NEIGHB_EXCL_SAME_MOL && 
                ptr[n].molindex == ptr[m].molindex &&
                ptr[n].molindex != -1 ){
        continue;
      }
      
      r2 = 0.0;
      for ( k=0; k<3; k++ ){
        r[k]  = ptr[n].x[k]-ptr[m].x[k];
        sep_Wrap( r[k], sys->length[k] );
        r2   += r[k]*r[k];
      }

      if ( r2 < cf2 ){ 
        ptr[n].neighb[index] = m;
        index++;
      }
      
    }
  }

}

// LLIST + neigbour list

void sep_neighb(seppart *ptr, sepsys *sys){
  int *list;

  list = sep_allocate_celllist(sys);
  sep_make_celllist(ptr, list, sys);
  sep_make_neighblist_from_llist(ptr, SEP_NEIGHB, list, sys);

  free(list);

}

void sep_neighb_excl_same_mol(seppart *ptr, sepsys *sys){
  int *list;

  list = sep_allocate_celllist(sys);
  sep_make_celllist(ptr, list, sys);
  sep_make_neighblist_from_llist_excl_same_mol(ptr, SEP_NEIGHB, list, sys);
  free(list);

}


void sep_neighb_nonbonded(seppart *ptr, sepsys *sys){
  int *list;

  list = sep_allocate_celllist(sys);
  sep_make_celllist(ptr, list, sys);
  sep_make_neighblist_from_llist_nonbonded(ptr, SEP_NEIGHB, list, sys);

  free(list);

}


int *sep_allocate_celllist(sepsys *sys){
  int *ptr;

  ptr = sep_vector_int(sys->npart+
		       sys->nsubbox[0]*sys->nsubbox[1]*sys->nsubbox[2]);

  if ( ptr == NULL )
    sep_error("%s at %d: Memory allocation error");

  return ptr;
}


void sep_make_celllist(seppart *ptr, int *list, sepsys *sys){
  int n, i, nsubbox3, nsubbox2, length;

  nsubbox2 = sys->nsubbox[0]*sys->nsubbox[1];
  nsubbox3 = nsubbox2*sys->nsubbox[2];
  length = sys->npart + nsubbox3;

  for ( n=0; n<length; n++ ) list[n] = -1;

  for (n=0; n<sys->npart; n++){
    i = (int)(ptr[n].x[0]/sys->lsubbox[0]) + 
      (int)(ptr[n].x[1]/sys->lsubbox[1])*sys->nsubbox[0] +
      (int)(ptr[n].x[2]/sys->lsubbox[2])*nsubbox2;

    if ( i > length - 1 )
      sep_error("%s at %d: Index larger than array length", __func__, __LINE__);

    list[n+nsubbox3] = list[i];
    list[i] = n;  
  }

}



void sep_make_neighblist_from_llist(seppart *ptr,  int nneighb, 
				    int *list, sepsys *sys){
  double dr, r2, cf2;
  int j1, j2, m1, m1X, m1Y, m1Z, m2, m2X, m2Y, m2Z,
    n, k, i, offset, nsubbox3, nsubbox2;
  static const int iofX[] = {0,1,1,0,-1,0,1,1,0,-1,-1,-1, 0, 1}; 
  static const int iofY[] = {0,0,1,1, 1,0,0,1,1, 1, 0,-1,-1,-1};
  static const int iofZ[] = {0,0,0,0, 0,1,1,1,1, 1, 1, 1, 1, 1}; 
  int *index, *icc, nccell;

  index = sep_vector_int(sys->npart);
  icc = sep_vector_int(sys->npart/3);  
  nsubbox2 = sys->nsubbox[1]*sys->nsubbox[0]; 
  nsubbox3 = sys->nsubbox[2]*nsubbox2;
  cf2 = sep_Sq(sys->cf + sys->skin);

#pragma omp parallel for private(k)
  for (n=0; n<sys->npart; n++)
    for (k=0; k<nneighb; k++)  ptr[n].neighb[k] = -1;

  
  for (m1Z = 0; m1Z < sys->nsubbox[2]; m1Z++){
    for (m1Y = 0; m1Y < sys->nsubbox[1]; m1Y++) {
      for (m1X = 0; m1X < sys->nsubbox[0]; m1X++) {

	// cell index
	m1 = m1Z*nsubbox2 + m1Y*sys->nsubbox[0] + m1X;

	// Re-organize
	if ( list[m1] == -1 ) continue;
	else j1 = list[m1];

	nccell = 0;
	while ( j1 != -1 ) {
	  icc[nccell] = j1;
	  nccell++;
	  j1 = list[j1+nsubbox3];
	}

#pragma omp parallel for schedule(dynamic)	\
  private(offset, j1, m2X, m2Y, m2Z, m2, j2, r2, dr, k)
	for ( i=0; i<nccell; i++ ){
	  
	  j1 = icc[i];
	  
	  // Over neighb. cells
	  for ( offset = 0; offset < 14; offset++ ) { 
	  
	    m2X = m1X + iofX[offset];    
	    if (m2X == sys->nsubbox[0] ) m2X = 0;    
	    else if (m2X == -1 ) m2X = sys->nsubbox[0]-1;    
	    
	    m2Y = m1Y + iofY[offset];    
	    if (m2Y == sys->nsubbox[1] )  m2Y = 0;    
	    else if (m2Y == -1 ) m2Y = sys->nsubbox[1]-1;    
	    
	    m2Z = m1Z + iofZ[offset];    
	    if ( m2Z == sys->nsubbox[2] ) m2Z = 0;    
	    
	    // Neighb. cell index
	    m2 = m2Z*nsubbox2 + m2Y*sys->nsubbox[0] + m2X;
	    
	    // Head-index for nieghb. cell
	    j2 = list[m2];
 
	    // Loop over particles in neighb. cell
	    while ( j2 != -1 ) {
	      
	      if ( m1 != m2 || j2 > j1 ) {
		r2 = 0.0;
		for ( k = 0; k<3; k++ ){
		  dr = ptr[j1].x[k] - ptr[j2].x[k];
		  sep_Wrap( dr, sys->length[k] );
		  r2 += dr*dr;
		}
		
		if ( r2 < cf2 ) {
		  ptr[j1].neighb[index[j1]] = j2;
		  index[j1]++;
		}	  
		
		if ( index[j1] == nneighb )
		  sep_error("%s at %d: Too many neighbours\n",
			    __func__, __LINE__);
		
	      }
	      // Get next particle in list for cell m2 
	      j2 = list[j2+nsubbox3]; 
	    }
	  } 
	}
      } } }
  
  free(index);
  free(icc);
}



void sep_make_neighblist_from_llist_nonbonded(seppart *ptr, int nneighb, 
					      int *list, sepsys *sys){
  double dr, r2, cf2;
  int j1, j2, m1, m1X, m1Y, m1Z, m2, m2X, m2Y, m2Z,
    i, n, k, offset, nsubbox3, nsubbox2, nccell;
  static int iofX[] = {0,1,1,0,-1,0,1,1,0,-1,-1,-1,0,1}; 
  static int iofY[] = {0,0,1,1,1,0,0,1,1,1,0,-1,-1,-1};
  static int iofZ[] = {0,0,0,0,0,1,1,1,1,1,1,1,1,1}; 
  int *index, *icc;

  icc = sep_vector_int(sys->npart/3);  
  index = sep_vector_int(sys->npart);  
  nsubbox2 = sys->nsubbox[1]*sys->nsubbox[0]; 
  nsubbox3 = sys->nsubbox[2]*nsubbox2;
  cf2 = sep_Sq(sys->cf + sys->skin);
    
  for (n=0; n<sys->npart; n++)
    for (k=0; k<nneighb; k++)   ptr[n].neighb[k] = -1;
    
  for (m1Z = 0; m1Z < sys->nsubbox[2]; m1Z++){
    for (m1Y = 0; m1Y < sys->nsubbox[1]; m1Y++) {
      for (m1X = 0; m1X < sys->nsubbox[0]; m1X++) {

	m1 = m1Z*nsubbox2 + m1Y*sys->nsubbox[0] + m1X;

	// Re-organize
	if ( list[m1] == -1 ) continue;
	else j1 = list[m1];

	nccell = 0;
	while ( j1 != -1 ) {
	  icc[nccell] = j1;
	  nccell++;
	  j1 = list[j1+nsubbox3];
	}

#pragma omp parallel for schedule(dynamic)	\
  private(offset, j1, m2X, m2Y, m2Z, m2, j2, r2, dr, k)
	for ( i=0; i<nccell; i++ ){
	  
	  j1 = icc[i];
	
	  for (offset = 0; offset < 14; offset++) {
	    
	    m2X = m1X + iofX[offset];    
	    if (m2X == sys->nsubbox[0] ) m2X = 0;    
	    else if (m2X == -1 ) m2X = sys->nsubbox[0]-1;    
	    
	    
	    m2Y = m1Y + iofY[offset];    
	    if (m2Y == sys->nsubbox[1] )  m2Y = 0;    
	    else if (m2Y == -1 ) m2Y = sys->nsubbox[1]-1;    
	    
	    m2Z = m1Z + iofZ[offset];    
	    if ( m2Z == sys->nsubbox[2] ) m2Z = 0;    

	    m2 = m2Z*nsubbox2 + m2Y*sys->nsubbox[0] + m2X;
	    
	    j2 = list[m2];
	    while ( j2 != -1 ) {
	      if ( m1 != m2 || j2 > j1 ) {
		r2 = 0.0;
		for ( k = 0; k<3; k++ ){
		  dr = ptr[j1].x[k] - ptr[j2].x[k];
		  sep_Wrap( dr, sys->length[k] );
		  r2 += dr*dr;
		}
		if ( r2 < cf2 && sep_bonded_direct(ptr, j1, j2) == 0 ) {
		  ptr[j1].neighb[index[j1]] = j2;
		  index[j1]++;
		  if ( index[j1] == nneighb ){
		    sep_error("%s at %d: Too many neighbours\n",
			      __func__, __LINE__);
		  }
		   
		}	
	      }
	      j2 = list[j2+nsubbox3];
	    }
	    
	  }
	}
      } } }
    

  free(index);
  free(icc);
  
}


void sep_make_neighblist_from_llist_excl_same_mol(seppart *ptr,  int nneighb, 
						  int *list, sepsys *sys){
  double dr, r2, cf2;
  int j1, j2, m1, m1X, m1Y, m1Z, m2, m2X, m2Y, m2Z,
    n, k, i, offset, nsubbox3, nsubbox2;
  static const int iofX[] = {0,1,1,0,-1,0,1,1,0,-1,-1,-1, 0, 1}; 
  static const int iofY[] = {0,0,1,1, 1,0,0,1,1, 1, 0,-1,-1,-1};
  static const int iofZ[] = {0,0,0,0, 0,1,1,1,1, 1, 1, 1, 1, 1}; 
  int *index, *icc, nccell;

  index = sep_vector_int(sys->npart);
  icc = sep_vector_int(sys->npart/3);  
  nsubbox2 = sys->nsubbox[1]*sys->nsubbox[0]; 
  nsubbox3 = sys->nsubbox[2]*nsubbox2;
  cf2 = sep_Sq(sys->cf + sys->skin);

#pragma omp parallel for private(k)
  for (n=0; n<sys->npart; n++)
    for (k=0; k<nneighb; k++)  ptr[n].neighb[k] = -1;

  
  for (m1Z = 0; m1Z < sys->nsubbox[2]; m1Z++){
    for (m1Y = 0; m1Y < sys->nsubbox[1]; m1Y++) {
      for (m1X = 0; m1X < sys->nsubbox[0]; m1X++) {

	// cell index
	m1 = m1Z*nsubbox2 + m1Y*sys->nsubbox[0] + m1X;

	// Re-organize
	if ( list[m1] == -1 ) continue;
	else j1 = list[m1];

	nccell = 0;
	while ( j1 != -1 ) {
	  icc[nccell] = j1;
	  nccell++;
	  j1 = list[j1+nsubbox3];
	}

#pragma omp parallel for schedule(dynamic)	\
  private(offset, j1, m2X, m2Y, m2Z, m2, j2, r2, dr, k)
	for ( i=0; i<nccell; i++ ){
	  
	  j1 = icc[i];
	  
	  // Over neighb. cells
	  for ( offset = 0; offset < 14; offset++ ) { 
	  
	    m2X = m1X + iofX[offset];    
	    if (m2X == sys->nsubbox[0] ) m2X = 0;    
	    else if (m2X == -1 ) m2X = sys->nsubbox[0]-1;    
	    
	    m2Y = m1Y + iofY[offset];    
	    if (m2Y == sys->nsubbox[1] )  m2Y = 0;    
	    else if (m2Y == -1 ) m2Y = sys->nsubbox[1]-1;    
	    
	    m2Z = m1Z + iofZ[offset];    
	    if ( m2Z == sys->nsubbox[2] ) m2Z = 0;    
	    
	    // Neighb. cell index
	    m2 = m2Z*nsubbox2 + m2Y*sys->nsubbox[0] + m2X;
	    
	    // Head-index for nieghb. cell
	    j2 = list[m2];
 
	    // Loop over particles in neighb. cell
	    while ( j2 != -1 ) {
	      
	      if ( (m1 != m2 || j2 > j1) && ( ptr[j1].molindex == -1 || 
			ptr[j1].molindex != ptr[j2].molindex ) ){
		r2 = 0.0;
		for ( k = 0; k<3; k++ ){
		  dr = ptr[j1].x[k] - ptr[j2].x[k];
		  sep_Wrap( dr, sys->length[k] );
		  r2 += dr*dr;
		}
		
		if ( r2 < cf2 ) {
		  ptr[j1].neighb[index[j1]] = j2;
		  index[j1]++;
		}	  
		
		if ( index[j1] == nneighb )
		  sep_error("%s at %d: Too many neighbours\n",
			    __func__, __LINE__);
		
	      }
	      j2 = list[j2+nsubbox3];
	    }
	  }
	}
      } } }
  
  free(index);
  free(icc);
}


unsigned int sep_bonded_direct(seppart *ptr, int j1, int j2) {

  for ( int k=0; k<SEP_BOND; k++ ) {   
    if ( ptr[j1].bond[k] == j2 || ptr[j2].bond[k] == j1 ) 
      return 1;
  }
  
  return 0;
}



	
////// LJ specific pair calculations

void sep_force_lj(seppart *ptr, const char *types, 
		  const double *p, sepsys *sys, 
		  sepret *retval, const unsigned opt){

  const double cf = p[0];
  if ( cf > sys->cf )
    sep_error("cutoff for an interaction cannot be larger than maximum cutoff");
  
  if ( sys->neighb_update == SEP_BRUTE ) {
    if (  sys->omp_flag == 1 ){
      sep_warning("omp flag set, SEP_BRUTE does not support threads.");
      sep_warning("Resetting omp flag");
      sys->omp_flag = 0;
    }
    sep_lj_pair_brute(ptr, types, p, sys, retval, opt);    
  }
  else {
    if ( sys->neighb_flag == 1 ){
      
      if ( sys->neighb_update == SEP_NEIGHBLIST ){
	sep_make_neighblist(ptr, sys, opt);
      }
      else if ( sys->neighb_update == SEP_LLIST_NEIGHBLIST ){
	
	if ( opt == SEP_ALL )
	  sep_neighb(ptr, sys);
	else if ( opt == SEP_NEIGHB_EXCL_BONDED )
	  sep_neighb_nonbonded(ptr, sys);
	else if ( opt == SEP_NEIGHB_EXCL_SAME_MOL )
	  sep_neighb_excl_same_mol(ptr, sys);
      }

      sys->neighb_flag = 0;
    }
    
    sep_lj_pair_neighb(ptr, types, p, sys, retval, sys->omp_flag);
      
  }

}

void sep_lj_pair_neighb(seppart *ptr, const char *types,
			const double *param, sepsys *sys, 
			sepret *retval, bool parallel) {
  int i1, i2, n, k, kk;
  double r2, ft, f[3], r[3];
  const double cf = param[0], eps=param[1], sigma=param[2], aw=param[3];
  const double shift = 4.0*eps*(pow(sigma/cf, 12.) - aw*pow(sigma/cf,6.));
  const double cf2 = cf*cf;
  size_t lvec = sys->npart*3;
 
  double epot = 0.0;
  int moli_i1, moli_i2;

  const double eps48 = 48.0*eps;
  const double eps4 = 4.0*eps;
  const double awh = 0.5*aw;
  const double sigmasqr = sigma*sigma;
  
  if ( parallel ){

    double *force_array = sep_vector(lvec);
    double *pconf = sep_vector(9);

#pragma omp parallel for schedule(dynamic)			\
  private(i1, n, i2, k, kk, r, r2, ft, f, moli_i1, moli_i2)	\
  reduction(+:epot, force_array[:lvec], pconf[:9]) 
    for (i1=0; i1<sys->npart; i1++){
      
      if ( ptr[i1].type != types[0] && ptr[i1].type != types[1] )
	continue;
      
      moli_i1 = ptr[i1].molindex;
      
      n = 0;
      while (1) {
	
	i2 = ptr[i1].neighb[n];
	if ( i2 == -1 ) break; 
	
	if ( (ptr[i1].type == types[0] && ptr[i2].type == types[1]) || 
	     (ptr[i1].type == types[1] && ptr[i2].type == types[0]) ){
	  r2 = 0.0;
	  for ( k=0; k<3; k++ ){
	    r[k] = ptr[i1].x[k] - ptr[i2].x[k];
	    sep_Wrap( r[k], sys->length[k] );
	    r2 += r[k]*r[k];
	  }
	  
	  if ( r2 < cf2 ){

	    double rri = sigmasqr/r2; double rri3 = rri*rri*rri;
	    ft = eps48*rri3*(rri3 - awh)*rri;
	 
	    for ( k=0; k<3; k++ ){
	      f[k] = ft*r[k];
	      
	      force_array[i1*3 + k] += f[k];
	      force_array[i2*3 + k] += -f[k];
	    }

	    epot += eps4*rri3*(rri3 - aw) - shift; 
	    	    
	    for (k=0; k<3; k++)
	      for ( kk=0; kk<3; kk++ )
		pconf[k*3+kk] += f[k]*r[kk];
	  }
	     
	}
	n++;
      }
    }

    for ( int k=0; k<3; k++ )
      for ( int kk=0; kk<3; kk++ ) retval->pot_P[k][kk] += pconf[3*k + kk];

    free(pconf);

    for ( int n=0; n<sys->npart; n++ )
      for ( int k=0; k<3; k++ )
	ptr[n].f[k] += force_array[n*3 + k];
    
    free(force_array);
  }
  
  else {
   
    for (i1=0; i1<sys->npart; i1++){
    
      if ( ptr[i1].type != types[0] && ptr[i1].type != types[1] )
	continue;
      
      moli_i1 = ptr[i1].molindex;
      n = 0;
      while (1) {
	
	i2 = ptr[i1].neighb[n];
	if ( i2 == -1 ) break; 
	
	if ( (ptr[i1].type == types[0] && ptr[i2].type == types[1]) || 
	     (ptr[i1].type == types[1] && ptr[i2].type == types[0]) ){

	  r2 = 0.0;
	  for ( k=0; k<3; k++ ){
	    r[k] = ptr[i1].x[k] - ptr[i2].x[k];
	    sep_Wrap( r[k], sys->length[k] );
	    r2 += r[k]*r[k];
	  }
	  
	  if ( r2 < cf2 ){

	    double rri = sigmasqr/r2; double rri3 = rri*rri*rri;
	    ft = eps48*rri3*(rri3 - awh)*rri;

	    for ( k=0; k<3; k++ ){
	      f[k] = ft*r[k];
	      
	      ptr[i1].f[k] += f[k];
	      ptr[i2].f[k] += -f[k];
	    }
	    epot += eps4*rri3*(rri3 - aw) - shift; 
	    
	    for (k=0; k<3; k++)
	      for ( kk=0; kk<3; kk++ )
		retval->pot_P[k][kk] += f[k]*r[kk];
	    
	    if ( sys->molptr->flag_Fij == 1 ){
	      moli_i2 = ptr[i2].molindex;
	      if ( (moli_i1 != -1 && moli_i2 != -1) && (moli_i1 != moli_i2) ){
		for ( k=0; k<3;k++ ){
		  sys->molptr->Fij[moli_i1][moli_i2][k] += f[k];
		  sys->molptr->Fij[moli_i2][moli_i1][k] -= f[k];	
		}
	      }	
	    }
	    
	  }
	}
	n++;
      }
    }
  }

  retval->epot += epot;
}
 
void sep_lj_pair_brute(seppart *ptr, const char *types, const double *p, 
		       sepsys *sys, sepret *retval, const int opt){
  int n,m,k;
  double r[3], r2, ft, f[3]; 
  const double cf = p[0], eps=p[1], sigma=p[2], aw=p[3];
 
  const double shift = 4.0*eps*(pow(sigma/cf, 12.) - aw*pow(sigma/cf,6.));
  const double cf2 = cf*cf;

  for (n=0; n<sys->npart-1; n++){
		
    int i = ptr[n].molindex;
					
    for (m=n+1; m<sys->npart; m++){

      if ( opt == SEP_NEIGHB_EXCL_BONDED && sep_bonded_direct(ptr, n, m) == 1 ){
	continue;
      }
      else if ( opt == SEP_NEIGHB_EXCL_SAME_MOL && 
		ptr[n].molindex == ptr[m].molindex &&
		ptr[n].molindex != -1 )
	continue;


      if ( (ptr[n].type == types[0] && ptr[m].type == types[1]) || 
	   (ptr[n].type == types[1] && ptr[m].type == types[0]) ){
				
	r2 = 0.0;
	for (k=0; k<3; k++){
	  r[k]  = ptr[n].x[k]-ptr[m].x[k];
	  sep_Wrap( r[k], sys->length[k] );
	  r2   += r[k]*r[k];
	}
	if (r2 < cf2){ 

	  // Force between particles
	  double rri = sigma*sigma/r2; double rri3 = rri*rri*rri;
	  ft = 48.0*eps*rri3*(rri3 - 0.5*aw)*rri;
	 
	  for (k=0; k<3; k++){
	    f[k] = ft*r[k];
						
	    ptr[n].f[k] += f[k];
	    ptr[m].f[k] -= f[k];
	  }
				
	  // Energy
	  retval->epot += 4.0*eps*rri3*(rri3 - aw) - shift; 
 	 					
	  // Config. part of the stress/pressure tensor 
	  for (k=0; k<3; k++)
	    for ( int kk=0; kk<3; kk++ ) retval->pot_P[k][kk] += f[k]*r[kk];
	
	  // Force between molecules
	  if ( sys->molptr->flag_Fij == 1 ){
	    int j = ptr[m].molindex;
	    if ( i != -1 && j != -1 ){
						
	      for ( k=0; k<3;k++ ){
		// pressure tensor	
		sys->molptr->Fij[i][j][k] += f[k];
		sys->molptr->Fij[j][i][k] -= f[k];
#ifdef COUPLE_TENSOR
		sys->molptr->Fiajb[n][j][k] += f[k];
		sys->molptr->Fiajb[m][i][k] -= f[k];
#endif
	      }
	    }	
	  }
					
	} // End of if ( r2 < rcut )
      }
    }
  }

}


/*
 * DPD Method from Groot and Warren 
 */  	

void sep_dpdforce_neighb(seppart *ptr, const char *types, 
			 const double cf, const double aij, 
			 const double temp_desired, 
			 const double sigma, sepsys *sys, sepret *retval,
			 const unsigned opt){
  int i1, i2, n, k;
  double r2,  r[3], rhat[3], vij[3], fC[3], fD[3], fR[3], dij, one_dij, randnum;
  const double cf2 = cf*cf, isqrtdt = 1.0/sqrt(sys->dt);
  const double facchk = 2.0*sqrt(3.0);
  const double gamma = sigma*sigma/(2.0*temp_desired);
  const size_t lvec = sys->npart*3;
  double *force_array = sep_vector(lvec);
  double dotrv, epot=.0;
  
  if ( sys->neighb_flag == 1 ){
    
    if ( opt == SEP_ALL )	
      sep_neighb(ptr, sys);
    else if ( opt == SEP_NEIGHB_EXCL_BONDED )
      sep_neighb_nonbonded(ptr, sys);
    else if ( opt == SEP_NEIGHB_EXCL_SAME_MOL )
      sep_neighb_excl_same_mol(ptr, sys);
  
    sys->neighb_flag = 0;
  }

#pragma omp parallel for schedule(dynamic)			\
  private(i1, n, i2, k, r, r2, dij, one_dij, dotrv, rhat, vij, randnum, fC, fD, fR) \
  reduction(+:epot, force_array[:lvec]) 
  for (i1=0; i1< sys->npart; i1++){

    if ( ptr[i1].type != types[0] && ptr[i1].type != types[1] ){
      continue;
    }

    n = 0;
    while (1){
      i2 = ptr[i1].neighb[n];
      if ( i2 == -1 ) break; 

      if ( (ptr[i1].type == types[0] && ptr[i2].type == types[1]) || 
	   (ptr[i1].type == types[1] && ptr[i2].type == types[0]) ){

	r2 = 0.0;
	for ( k=0; k<3; k++ ){
	  r[k] = ptr[i1].x[k] - ptr[i2].x[k];
	  sep_Wrap( r[k], sys->length[k] );
	  r2 += r[k]*r[k];
	}
	
	if ( r2 < cf2 ){

	  dij = sqrt(r2);
	  one_dij = (1.0-dij);

	  dotrv = 0.0;
	  for ( k=0; k<3; k++ ){
	    rhat[k] = r[k]/dij;
	    vij[k] = ptr[i1].pv[k] - ptr[i2].pv[k];
	    dotrv += rhat[k]*vij[k];
	  }

	  randnum = (sep_rand()-0.5)*facchk;

	  for ( k=0; k<3; k++ ){
	    // Conservative force
	    fC[k] = aij*one_dij*rhat[k];

	    // Dissipative force
	    fD[k] = -gamma*one_dij*one_dij*dotrv*rhat[k];

	    // Random force
	    fR[k] = sigma*one_dij*rhat[k]*isqrtdt*(randnum); 

	    // Summing up
	    force_array[i1*3 + k] += fC[k] + fD[k] + fR[k];
	    force_array[i2*3 + k] -= fC[k] + fD[k] + fR[k];
	  }

	  // potential energy (conservative force)
	  epot += 0.5*aij*one_dij*one_dij;

	  /*
	  pressure
	  for (k=0; k<3; k++){
	    for ( int kk=0; kk<3; kk++ ){

	      retval->pot_P[k][kk] += fC[k]*r[kk];

	      retval->pot_P_conservative[k][kk] += fC[k]*r[kk];
	      retval->pot_P_random[k][kk]       += fR[k]*r[kk];
	      retval->pot_P_dissipative[k][kk]  += fD[k]*r[kk];
	    }
	  }
	  
	  if ( sys->molptr->flag_Fij == 1 ){
	    int moli_i2 = ptr[i2].molindex;
	    if ( moli_i1 != -1 && moli_i2 != -1 ){
	      for ( int k=0; k<3;k++ ){
		sys->molptr->Fij[moli_i1][moli_i2][k] += fC[k];
		sys->molptr->Fij[moli_i2][moli_i1][k] -= fC[k];	

#ifdef COUPLE_TENSOR
		sys->molptr->Fiajb[i1][moli_i2][k] += fC[k];
		sys->molptr->Fiajb[i2][moli_i1][k] -= fC[k];
#endif
	      }
	    }	
	  }
	  */

	}  //  if r^2 < cutoff
	
      }    //  if atom type is considered
      n++;
    }
  }        //  neighbor list loops


  for ( int n=0; n<sys->npart; n++ )
    for ( int k=0; k<3; k++ )
      ptr[n].f[k] += force_array[n*3 + k];
  
  free(force_array);

  retval->epot = epot;
}

void sep_dpdforce_brute(seppart *ptr, const char *types, 
			 const double cf, const double aij, 
			 const double temp_desired, 
                         const double sigma, sepsys *sys, sepret *retval){
  int i1, i2, k;
  double r2,  r[3], rhat[3], vij[3], fC[3], fD[3], fR[3], dij, one_dij, randnum;
  const double cf2 = cf*cf, isqrtdt = 1.0/sqrt(sys->dt);
  const double facchk = 2.0*sqrt(3.0);
      
  const double gamma = sigma*sigma/(2.0*temp_desired);

  for (i1=0; i1< sys->npart-1; i1++){

    if ( ptr[i1].type != types[0] && ptr[i1].type != types[1] ){
      continue;
    }

    int moli_i1 = ptr[i1].molindex;
    
    for (i2=i1+1; i2< sys->npart; i2++){
      if ( i2 == -1 ) break; 

      if ( (ptr[i1].type == types[0] && ptr[i2].type == types[1]) || 
	   (ptr[i1].type == types[1] && ptr[i2].type == types[0]) ){

	r2 = 0.0;
	for ( k=0; k<3; k++ ){
	  r[k] = ptr[i1].x[k] - ptr[i2].x[k];
	  sep_Wrap( r[k], sys->length[k] );
	  r2 += r[k]*r[k];
	}
	
	if ( r2 < cf2 ){

	  dij = sqrt(r2);
	  one_dij = (1.0-dij);
      
	  double dotrv = 0.0;
	  for ( int k=0; k<3; k++ ){
	    // r[k] /= dij;
	    rhat[k] = r[k]/dij;
	    vij[k] = ptr[i1].pv[k] - ptr[i2].pv[k];
	    dotrv += rhat[k]*vij[k];
	  }
	  
	  randnum = (sep_rand()-0.5)*facchk;
	  for ( int k=0; k<3; k++ ){
	    // Conservative force
	    fC[k] = aij*one_dij*rhat[k];

	    // Dissipative force
	    fD[k] = -gamma*one_dij*one_dij*dotrv*rhat[k];

	    // Random force
	    fR[k] = sigma*one_dij*rhat[k]*isqrtdt*(randnum); 

	    // Summing up
	    ptr[i1].f[k] += fC[k] + fD[k] + fR[k];
	    ptr[i2].f[k] -= fC[k] + fD[k] + fR[k];
	  }

	  // potential energy (conservative force)
	  retval->epot += 0.5*aij*one_dij*one_dij;

	  // pressure
	  for (k=0; k<3; k++){
	    for ( int kk=0; kk<3; kk++ ){
	      retval->pot_P[k][kk] += fC[k]*r[kk];

	      retval->pot_P_conservative[k][kk] += fC[k]*r[kk];
	      retval->pot_P_random[k][kk]       += fR[k]*r[kk];
	      retval->pot_P_dissipative[k][kk]  += fD[k]*r[kk];
	    }
	  }
	
	  if ( sys->molptr->flag_Fij == 1 ){
	    int moli_i2 = ptr[i2].molindex;
	    if ( moli_i1 != -1 && moli_i2 != -1 ){
	      for ( int k=0; k<3;k++ ){
		sys->molptr->Fij[moli_i1][moli_i2][k] += fC[k];
		sys->molptr->Fij[moli_i2][moli_i1][k] -= fC[k];	
#ifdef COUPLE_TENSOR
		sys->molptr->Fiajb[i1][moli_i2][k] += fC[k];
		sys->molptr->Fiajb[i2][moli_i1][k] -= fC[k];
#endif
	      }
	    }	
	  }
	}
	
	//  if r^2 < cutoff
	
      }    //  if atom type is considered
    }
  }        //  neighbor list loops

}





void _sep_make_neighblist_from_llist(seppart *ptr,  int nneighb, 
				     int *list, sepsys sys){
  double r[3], r2, cf2;
  int j1, j2, m1, m1X, m1Y, m1Z, m2, m2X, m2Y, m2Z,
    n, k, offset, nsubbox3, nsubbox2;
  static const int iofX[] = {0,1,1,0,-1,0,1,1,0,-1,-1,-1, 0, 1}; 
  static const int iofY[] = {0,0,1,1, 1,0,0,1,1, 1, 0,-1,-1,-1};
  static const int iofZ[] = {0,0,0,0, 0,1,1,1,1, 1, 1, 1, 1, 1}; 
  int *index;

  index = sep_vector_int(sys.npart);  
  nsubbox2 = sys.nsubbox[1]*sys.nsubbox[0]; 
  nsubbox3 = sys.nsubbox[2]*nsubbox2;
  cf2 = sep_Sq(sys.cf + sys.skin);

  for (n=0; n<sys.npart; n++)
    for (k=0; k<nneighb; k++)  ptr[n].neighb[k] = -1;
  
  for (m1Z = 0; m1Z < sys.nsubbox[2]; m1Z++){
    for (m1Y = 0; m1Y < sys.nsubbox[1]; m1Y++) {
      for (m1X = 0; m1X < sys.nsubbox[0]; m1X++) {

	m1 = m1Z*nsubbox2 + m1Y*sys.nsubbox[0] + m1X;
	
	for (offset = 0; offset < 14; offset++) {

	  m2X = m1X + iofX[offset];    
	  if (m2X == sys.nsubbox[0] ) m2X = 0;    
	  else if (m2X == -1 ) m2X = sys.nsubbox[0]-1;    
    
	  m2Y = m1Y + iofY[offset];    
	  if (m2Y == sys.nsubbox[1] )  m2Y = 0;    
	  else if (m2Y == -1 ) m2Y = sys.nsubbox[1]-1;    
	  
	  m2Z = m1Z + iofZ[offset];    
	  if ( m2Z == sys.nsubbox[2] ) m2Z = 0;    

	  m2 = m2Z*nsubbox2 + m2Y*sys.nsubbox[0] + m2X;

	  j1 = list[m1];

	  while ( j1 != -1 ) {
	    j2 = list[m2];
	    while ( j2 != -1 ) {
	      if ( m1 != m2 || j2 < j1 ) {
		r2 = 0.0;
		for (k = 0; k<3; k++){
		  r[k] = ptr[j1].x[k] - ptr[j2].x[k];
		  sep_Wrap( r[k], sys.length[k] );
		  r2 += r[k]*r[k];
		}

		if (r2 < cf2) {
		  ptr[j1].neighb[index[j1]] = j2;
		  index[j1]++;

		  if ( index[j1] > nneighb - 1 ){
		    sep_error("%s at %d: Too many neighbours\n",
			      __func__, __LINE__);
		  }
		  /*
		  if ( sys.omp_flag == 1 ){
		    ptr[j2].neighb[index[j2]] = j1;
		    index[j2]++;
		  }
		  */
		}
	      }
	      j2 = list[j2+nsubbox3];
	    }
	    j1 = list[j1+nsubbox3];  
	  }

	} } } }

  free(index);
}

