#include "sepcoulomb.h"



void sep_coulomb_sf(seppart *ptr, double cf, sepsys *sys, 
		    sepret *retval, const unsigned opt){

  if ( sys->neighb_update == 0 ) 
    sep_coulomb_sf_brute(ptr, cf, sys, retval, opt);    
  else {
    if ( sys->omp_flag ){
      sep_coulomb_sf_neighb_omp(ptr, cf, sys, retval);
    }
    else 
      sep_coulomb_sf_neighb(ptr, cf, sys, retval);
  }
		
}

void sep_coulomb_sf_brute(seppart *ptr,  double cf, sepsys *sys, 
			  sepret *retval, const int opt){
  int n,m,k;
  double dr[3], r2, ft, f[3], cf2, icf2, icf; 

  cf2 = cf*cf; 
  icf2 = 1.0/cf2; 
  icf = 1.0/cf;

  for ( n=0; n<sys->npart-1; n++ ){

    if ( fabs(ptr[n].z) < DBL_EPSILON ) continue;
    
    int i = ptr[n].molindex;		
    
    for ( m=n+1; m<sys->npart; m++ ){
      
      if ( opt == SEP_NEIGHB_EXCL_BONDED && sep_bonded(ptr, n, m) == 1 ){
	continue;	
      }
      else if ( opt == SEP_NEIGHB_EXCL_SAME_MOL && 
		ptr[n].molindex == ptr[m].molindex && ptr[n].molindex != -1 )
	continue;
      
      r2 = 0.0;
      for (k=0; k<3; k++){
	dr[k]  = ptr[n].x[k]-ptr[m].x[k];
	sep_Wrap( dr[k], sys->length[k] );
	r2   += dr[k]*dr[k];
      }
      
      if ( r2 < cf2 ){ 
	
	double zizj = ptr[n].z*ptr[m].z;
	double r = sqrt(r2);
	ft = zizj*(1.0/r2 - icf2)/r; 
	
	for ( k=0; k<3; k++ ){	
	  f[k] = ft*dr[k];
	  
	  ptr[n].f[k] += f[k];
	  ptr[m].f[k] -= f[k];
	}
	
	for (k=0; k<3; k++)
	  for ( int kk=0; kk<3; kk++ )
	    retval->pot_P[k][kk] += f[k]*dr[kk];
	
	// Force between molecules
	if ( sys->molptr->flag_Fij == 1 ){
	  
	  int j = ptr[m].molindex;
	  
	  if ( i != -1 && j != -1 ){
	    for ( k=0; k<3;k++ ){
	      // pressure tensor	
	      sys->molptr->Fij[i][j][k] += f[k];
	      sys->molptr->Fij[j][i][k] -= f[k];
	      // Couple tensor
	      //sys->molptr->Fiajb[n][j][k] += f[k];
	      //sys->molptr->Fiajb[m][i][k] -= f[k];	
	    }
	  }	
	  
	}
	
	double ecoul = zizj*(1.0/r + (r-cf)*icf2 - icf);
	retval->epot += ecoul;		
	retval->ecoul += ecoul; 
	
      } 
    }
  }

}
		    
void sep_coulomb_sf_neighb(seppart *ptr, double cf, sepsys *sys, sepret *retval) {
  int i1, i2, n, k;
  double r2, ft, f[3], dr[3];
  
  const double cf2 = cf*cf;
  const double icf2 = 1.0/cf2; 
  const double icf = 1.0/cf;


  for (i1=0; i1<sys->npart; i1++){
    
    if ( fabs(ptr[i1].z) < DBL_EPSILON ) continue;

    int moli_i1	= ptr[i1].molindex;
    n = 0;
    while (1){
      i2 = ptr[i1].neighb[n];
      
      if ( i2 == -1 ) break; 
      
      r2 = 0.0;
      for ( k=0; k<3; k++ ){
	dr[k] = ptr[i1].x[k] - ptr[i2].x[k];
	sep_Wrap( dr[k], sys->length[k] );
	r2 += dr[k]*dr[k];
      }
      
      if (r2 < cf2){ 
	
	double zizj = ptr[i2].z*ptr[i1].z;
	double r = sqrt(r2);
	ft = zizj*(1.0/r2 - icf2)/r; 
				
	for (k=0; k<3; k++){
	  f[k] = ft*dr[k];
					
	  ptr[i1].f[k] += f[k];
	  ptr[i2].f[k] -= f[k];
	}
				
	for (k=0; k<3; k++)
	  for ( int kk=0; kk<3; kk++ )
	    retval->pot_P[k][kk] += f[k]*dr[kk];
										
	if ( sys->molptr->flag_Fij == 1 ){
	  int moli_i2 = ptr[i2].molindex;
	  if ( moli_i1 != -1 && moli_i2 != -1 ){
	    for ( int k=0; k<3;k++ ){
	      sys->molptr->Fij[moli_i1][moli_i2][k] += f[k];
	      sys->molptr->Fij[moli_i2][moli_i1][k] -= f[k];
	      
	    }
	  }	
	}
					
	double ecoul = zizj*(1.0/r + (r-cf)*icf2 - icf);
	retval->epot +=  ecoul;
	retval->ecoul += ecoul;
				
      }
      n++;
    }
  }
  
}


void sep_coulomb_sf_neighb_omp(seppart *ptr, double cf, sepsys *sys, 
			       sepret *retval){
  int i1, i2, n, k, kk;
  double r[3], r2, ft, f[3], ecoul, rij, zizj;
  const double cf2 = cf*cf, icf2=1.0/cf2, icf=1.0/cf;
  size_t lvec = sys->npart*3;

  double *force_array = sep_vector(lvec);
  double *pconf = sep_vector(9);

  ecoul = 0.0;
  
#pragma omp parallel for schedule(dynamic)			\
  private(i1, n, i2, k, kk, r, r2, ft, f)			\
  reduction(+:ecoul, force_array[:lvec]) 
  for (i1=0; i1<sys->npart; i1++){
    
    if ( fabs(ptr[i1].z) < DBL_EPSILON ) continue;
    
    n = 0;
    while (1){
      i2 = ptr[i1].neighb[n];
      if ( i2 == -1 ) break; 
      
      for ( k=0; k<3; k++ ){
	r[k] = ptr[i1].x[k] - ptr[i2].x[k];
	sep_Wrap( r[k], sys->length[k] );
      }
      
      r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
      
      if ( r2 < cf2 ){	
	zizj = ptr[i2].z*ptr[i1].z;
	rij = sqrt(r2);
	ft = zizj*(1.0/r2 - icf2)/rij; 
	
	for ( k=0; k<3; k++ ){
	  f[k] = ft*r[k];

	  force_array[i1*3 + k] += f[k];
	  force_array[i2*3 + k] += -f[k];
	}
	
	ecoul +=  zizj*(1.0/rij + (rij-cf)*icf2 - icf);
	for ( k=0; k<3; k++ )
	  for ( kk=0; kk<3; kk++ )
	    pconf[k*3+kk] += f[k]*r[kk];
      }
      
      n++;
    }
  }

    
  for ( n=0; n<sys->npart; n++ )
    for ( k=0; k<3; k++ )
      ptr[n].f[k] += force_array[n*3 + k];
  
  for ( k=0; k<3; k++ )
    for ( kk=0; kk<3; kk++ ) retval->pot_P[k][kk] += pconf[3*k + kk];

  free(force_array); free(pconf);
  
  retval->epot  += 0.5*ecoul;
  retval->ecoul += 0.5*ecoul;
}
	

void sep_coulomb_wolf(seppart *ptr, double alpha, double cf, sepsys *sys, 
		      sepret *retval, const unsigned opt){
  
  if ( sys->neighb_update == 0 ) 
    sep_coulomb_wolf_brute(ptr, alpha, cf, sys, retval, opt);    
  else {
    //if ( sys->omp_flag==1 )
    //	sep_coloumb_wolf_neighb_omp(ptr, cf, sys, retval);
    //else 
    sep_coulomb_wolf_neighb(ptr, alpha, cf, sys, retval);
  }
}		


void sep_coulomb_wolf_brute(seppart *ptr, double alpha, double rcf, sepsys *sys,
			    sepret *ret, unsigned opt) {
  double rij2, f[3], fij, rij[3], ecoul;
  const double rcf2 = rcf*rcf, 
    fac = 2*alpha/sqrt(SEP_PI), alpha2=alpha*alpha, 
    facs = erfc(alpha*rcf)/(2*rcf) + alpha/sqrt(SEP_PI), 
    c = erfc(alpha*rcf)/rcf2, d = fac*exp(-alpha2*rcf2)/rcf,
    b = erfc(alpha*rcf)/rcf;
  
  ecoul= 0.0;
  
  for ( int n=0; n<sys->npart-1; n++ ){
    
    if ( fabs(ptr[n].z) < DBL_EPSILON ) continue;
    
    // self term 
    int i = ptr[n].molindex;
    if ( i == -1 ) ecoul -= facs*ptr[n].z*ptr[n].z;
    
    for ( int m=n+1; m<sys->npart; m++ ){
      
      if ( opt == SEP_NEIGHB_EXCL_BONDED && sep_bonded(ptr, n, m) == 1 )
	continue;	
      else if ( opt == SEP_NEIGHB_EXCL_SAME_MOL && 
		ptr[n].molindex == ptr[m].molindex && ptr[n].molindex != -1 )
	continue;
      
      rij2 = 0.0;
      for ( int k=0; k<3; k++ ){
        rij[k] = ptr[n].x[k] - ptr[m].x[k];
        sep_Wrap( rij[k], sys->length[k] );
        rij2 += rij[k]*rij[k];
      }
      
      if ( rij2 < rcf2 ){
        double rij1 = sqrt(rij2);
        double irij1 = 1.0/rij1;
        double zizj = ptr[n].z*ptr[m].z;
        
        double a = erfc(alpha*rij1)/rij2;
        double aa = fac*exp(-alpha2*rij2)*irij1;

        fij = zizj*(a + aa - c - d)*irij1;
        for ( int k=0; k<3; k++ ){
          f[k] = fij*rij[k];
          ptr[n].f[k] += f[k];
          ptr[m].f[k] -= f[k];
        }
        
        a = erfc(alpha*rij1)*irij1;
        
        if ( ptr[n].molindex == -1 )    // ions 
          ecoul += zizj*(a - b);
        else 
          ecoul += zizj*(a - b + (c + d)*(rij1 - rcf));  // Molecules 
        
        for ( int k=0; k<3; k++ )
          for ( int kk=0; kk<3; kk++ )
            ret->pot_P[k][kk] += f[k]*rij[kk];
	
	// Force between molecules
	if ( sys->molptr->flag_Fij == 1 ){
	  
	  int j = ptr[m].molindex;
	  
	  if ( i != -1 && j != -1 ){
	    for ( int k=0; k<3;k++ ){
	      // pressure tensor	
	      sys->molptr->Fij[i][j][k] += f[k];
	      sys->molptr->Fij[j][i][k] -= f[k];
	      // Couple tensor
	      sys->molptr->Fiajb[n][j][k] += f[k];
	      sys->molptr->Fiajb[m][i][k] -= f[k];	
	    }
	  }	
	}    
	
      }

    }
  }

  ret->ecoul += ecoul;
  ret->epot  += ecoul;

}
				 
void sep_coulomb_wolf_neighb(seppart *ptr, double alpha, double rcf, 
			     sep3D *sys, sepret *retval) {
  int i, j, n, k;
  double rij2, f[3], fij, rij[3], ecoul;
  const double rcf2 = rcf*rcf, 
    fac = 2*alpha/sqrt(SEP_PI), alpha2=alpha*alpha, 
    facs = erfc(alpha*rcf)/(2*rcf) + alpha/sqrt(SEP_PI), 
    c = erfc(alpha*rcf)/rcf2, d = fac*exp(-alpha2*rcf2)/rcf,
    b = erfc(alpha*rcf)/rcf;

  ecoul= 0.0;
  
  for ( i=0; i<sys->npart; i++ ){

    if ( fabs(ptr[i].z) < DBL_EPSILON ) continue;
    
    // self term 
    int moli_i = ptr[i].molindex;
    if ( moli_i == -1 )
      ecoul -= facs*ptr[i].z*ptr[i].z;

    // pair term
    n = 0;
    while (1){
      j = ptr[i].neighb[n];

      if ( j == -1 ) break;
     
      rij2 = 0.0;
      for ( k=0; k<3; k++ ){
        rij[k] = ptr[i].x[k] - ptr[j].x[k];
	sep_Wrap( rij[k], sys->length[k] );
        rij2 += rij[k]*rij[k];
      }

      if ( rij2 < rcf2 ){
        double rij1 = sqrt(rij2);
        double irij1 = 1.0/rij1;
        double zizj = ptr[i].z*ptr[j].z;
        
        double a = erfc(alpha*rij1)/rij2;
        double aa = fac*exp(-alpha2*rij2)*irij1;

        fij = zizj*(a + aa - c - d)*irij1;
        
        for ( k=0; k<3; k++ ){
          f[k] = fij*rij[k];
          ptr[i].f[k] += f[k];
          ptr[j].f[k] -= f[k];
        }
        
        a = erfc(alpha*rij1)*irij1;
        
        // ions 
        if ( ptr[i].molindex == -1 ) 
          ecoul += zizj*(a - b);
        // molecules
        else 
          ecoul += zizj*(a - b + (c + d)*(rij1 - rcf));  
      
	for (k=0; k<3; k++)
	  for ( int kk=0; kk<3; kk++ )
	    retval->pot_P[k][kk] += f[k]*rij[kk];
	
	if ( sys->molptr->flag_Fij == 1 ){
	  int moli_j = ptr[j].molindex;
	  if ( moli_i != -1 && moli_j != -1 ){
	    for ( int k=0; k<3;k++ ){
	      sys->molptr->Fij[moli_i][moli_j][k] += f[k];
	      sys->molptr->Fij[moli_j][moli_i][k] -= f[k];
	      
	      sys->molptr->Fiajb[i][moli_j][k] += f[k];
	      sys->molptr->Fiajb[j][moli_i][k] -= f[k];	
	    }
	  }	
	}
	
      }
      
      n++;
    }
    
  }

  retval->epot +=  ecoul;
  retval->ecoul += ecoul;

}



void sep_ewald_direct(sepatom *ptr, int nrep, const sepsys sys){
  double rij[3], r2, boxId[3];

  for ( int n=-nrep; n<=nrep; n++ ){
    boxId[0] = n; 
    for ( int m=-nrep; m<=nrep; m++ ){
      boxId[1] = m; 
      for ( int l=-nrep; l<=nrep; l++ ){
	boxId[2] = l;

	for ( int i=0; i<sys.npart-1; i++ ){
	  for ( int j=i+1; j<sys.npart; j++ ){

	    r2 = 0.0;
	    for ( int k=0; k<3; k++ ){
	      rij[k] = ptr[i].x[k] - ptr[j].x[k] + boxId[k]*sys.length[k];
	      r2    += rij[k]*rij[k];
	    }
	    double r = sqrt(r2);
	    double fac = ptr[i].z*ptr[j].z/(r2*r);
	    
	    for ( int k=0; k<3; k++ ){    
	      ptr[i].f[k] += fac*rij[k];
	      ptr[j].f[k] -= fac*rij[k];
	    }
	    
	  }
	}
	
      }
    }
  }
  
}


