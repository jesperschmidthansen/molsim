
/* 
* sepintgr.c - This file is a part of the sep-library 
*
* Copyright (C) 2008 Jesper Schmidt Hansen 
* 
* License: GPL - see COPYING for copying conditions.
* There is ABSOLUTELY NO WARRANTY, not even for MERCHANTIBILITY or
* FITNESS FOR A PARTICULAR PURPOSE.
*
*/

#include "sepintgr.h"
#include "sepinit.h"
#include "sepdef.h"

// Can be used by sep_intr.c functions in general
double sep_periodic(sepatom *atoms, unsigned n, sepsys *sys){
  
  double d2 = 0.0;
  for ( int k=0; k<3; k++ ){
    if ( atoms[n].x[k] > sys->length[k] ){
      atoms[n].x[k] -= sys->length[k];
      atoms[n].cross_neighb[k] ++;
      atoms[n].crossings[k] ++;
    }
    else if ( atoms[n].x[k] < 0.0 ){
      atoms[n].x[k] += sys->length[k];
      atoms[n].cross_neighb[k] --;
      atoms[n].crossings[k] --;
    }
  
    double ri = (atoms[n].x[k] + atoms[n].cross_neighb[k]*sys->length[k]) - atoms[n].xn[k];
    
    d2 += sep_Sq(ri);
  }
    
  return d2;
  
}

void sep_leapfrog(seppart *ptr, sepsys *sys, sepret *retval){

  double v[3], sumekin=0.0, d2;
  int n;
  const double dt = sys->dt;

  //  double max_d2 = 0.0;
	
  for ( n=0; n<sys->npart; n++ ){
    
    d2 = 0.0;
    for (int k=0; k<3; k++){
      ptr[n].a[k] = ptr[n].f[k]/ptr[n].m;
      ptr[n].v[k] += ptr[n].a[k]*dt;
      ptr[n].x[k] += ptr[n].v[k]*dt;
      
      v[k] = ptr[n].v[k] - 0.5*ptr[n].a[k]*dt;
      sumekin += v[k]*v[k]*ptr[n].m;
    }

    d2 = sep_periodic(ptr, n, sys);

    if ( d2 > sys->max_dist2 ) sys->max_dist2 = d2;
      
    for ( int k=0; k<3 ; k++ )
      for ( int kk=0; kk<3; kk++ )
	retval->kin_P[k][kk] += v[k]*v[kk]*ptr[n].m;
  }

  
  if ( sqrt(sys->max_dist2) > sys->skin*0.5 ){
    sys->neighb_flag = 1;
    sys->nupdate_neighb ++;
    
    sep_set_xn(ptr, sys->npart);
   
    for ( int n=0; n<sys->npart; n++ ){
      for ( int k=0; k<3; k++ ){
	ptr[n].cross_neighb[k] = 0;
      }
    }
    
  }

  sys->tnow += sys->dt;
  retval->ekin += 0.5*sumekin;
}


void sep_langevinGJF(sepatom *ptr, double temp0, double alpha, sepsys *sys, sepret *retval){
  
  const double dt = sys->dt;
  const double cc  = exp(-alpha*dt);
  
  double sum_ekin = 0.0;
  for ( unsigned n=0; n<sys->npart; n++ ){
    
    double mass = ptr[n].m;
	double imass = 1.0/mass;
    double imass2 = 0.5*imass;
    double fac = sqrt(temp0*(1.0-cc*cc)); // mass?
    
    double c = alpha*dt*imass2;
    double a = (1.0-c)/(1.0+c);
    double b = 1.0/(1.0+c);
      
    for ( unsigned k=0; k<3; k++ ){
      
      ptr[n].v[k] = a*ptr[n].v[k] + dt*imass2*(a*ptr[n].prevf[k]+ptr[n].f[k]) + b*imass*ptr[n].randn[k];
      ptr[n].a[k] = ptr[n].f[k]*imass;

      ptr[n].prevf[k] = ptr[n].f[k];
      
      sum_ekin += ptr[n].v[k]*ptr[n].v[k]*mass;
      
      ptr[n].randn[k] = fac*sep_randn();
	  
      ptr[n].x[k] += b*dt*ptr[n].v[k] + b*dt*dt*imass2*ptr[n].f[k] + b*dt*imass2*ptr[n].randn[k];

    }
    
    double d2 = sep_periodic(ptr, n, sys);
    if ( d2 > sys->max_dist2 ) sys->max_dist2 = d2;
    
    for ( int k=0; k<3 ; k++ )
      for ( int kk=0; kk<3; kk++ )
		  retval->kin_P[k][kk] += ptr[n].v[k]*ptr[n].v[kk]*mass;
  }

  if ( sqrt(sys->max_dist2) > sys->skin*0.5 ){
    sys->neighb_flag = 1;
    sys->nupdate_neighb ++;
    
    sep_set_xn(ptr, sys->npart);
    
    for ( int n=0; n<sys->npart; n++ ){
      for ( int k=0; k<3; k++ ){
	ptr[n].cross_neighb[k] = 0;
      }
    }
  }
  
  sys->tnow += sys->dt;
  retval->ekin += 0.5*sum_ekin;
}

  
void sep_nosehoover(sepatom *ptr, double temp0, double *alpha,
		    const double tau, sepsys *sys){
  
  double ekin = 0.0;
  for (int n=0; n<sys->npart; n++){
    for ( int k=0; k<3; k++ )
      ekin += sep_Sq(ptr[n].v[k])*ptr[n].m;
  }

  ekin = 0.5*ekin/sys->npart;
  double temp = 0.666667*ekin;

  *alpha = *alpha + sys->dt/(tau*tau)*(temp/temp0 - 1.0); 
  
  for ( int n=0; n<sys->npart; n++ ){
    for ( int k=0; k<3; k++ )
      ptr[n].f[k] -= *alpha*ptr[n].m*ptr[n].v[k];
  }
      
} 

void _sep_nosehoover_type(seppart *ptr, char type, double Td, 
			double *alpha, const double Q, sepsys *sys){ 

  int ntype = 0;  double ekin = 0.0;
  for (int n=0; n<sys->npart; n++){
    if ( ptr[n].type == type ){
      ntype++;
      for ( int k=0; k<3; k++ )
	ekin += sep_Sq(ptr[n].v[k])*ptr[n].m;
    }
  }

  const double g = 3*ntype-3;
  
  double tmp = alpha[0];

  alpha[0] = alpha[1];
  alpha[1] = alpha[2];
  alpha[2] = tmp + 2.0*sys->dt*(ekin - g*Td)/Q;	    

  for ( int n=0; n<sys->npart; n++ ){
    if ( ptr[n].type == type ) {	
      for ( int k=0; k<3; k++ ) { 
	ptr[n].f[k] -= alpha[1]*ptr[n].v[k]*ptr[n].m;
      }
    }
  }
  
}


void _sep_nosehoover(seppart *ptr, double Td, double *alpha, 
		    const double Q, sepsys *sys){ 
 
  double ekin = 0.0;
  for (int n=0; n<sys->npart; n++)
    for ( int k=0; k<3; k++ )
      ekin += sep_Sq(ptr[n].v[k])*ptr[n].m;
  
  const double g =  sys->ndof;
  double tmp = alpha[0];

  alpha[0] = alpha[1];
  alpha[1] = alpha[2];
  alpha[2] = tmp + 2.0*sys->dt*(ekin - g*Td)/Q;	    
  
  for ( int n=0; n<sys->npart; n++ )
    for ( int k=0; k<3; k++ ) 
      ptr[n].f[k] -= alpha[1]*ptr[n].v[k]*ptr[n].m;
  
}



static inline double sep_verlet_wrap(double x, double xp, double lbox){

  if ( x - xp > 0.5*lbox )
    x -= lbox;
  else if ( x - xp < -0.5*lbox )
    x += lbox;

  return x;

}

void sep_fp(seppart *ptr, double temp_desired, sepsys *sys, sepret *retval){
  double sumekin=0.0, d2;
  int n;
  const double dt = sys->dt;
  const double fac = sqrt(1.0/12.0);

  //double max_d2 = 0.0;
	
  for ( n=0; n<sys->npart; n++ ){
    
    d2 = 0.0;
    double im = 1.0/ptr[n].m;
    double fric = temp_desired/ptr[n].ldiff;
    double gaussfac = sqrt(24*temp_desired*fric/dt);

    for (int k=0; k<3; k++){
      
      double a = sep_randn()*fac*gaussfac ;
         
      ptr[n].x[k] += dt*ptr[n].v[k]; 
      ptr[n].v[k] += im*dt*(ptr[n].f[k] - fric*ptr[n].v[k] + a);
        
      if ( ptr[n].x[k] > sys->length[k] ){
	ptr[n].x[k] -= sys->length[k];
	ptr[n].cross_neighb[k] ++;
	ptr[n].crossings[k] ++;
      }
      else if ( ptr[n].x[k] < 0.0 ){
	ptr[n].x[k] += sys->length[k];
	ptr[n].cross_neighb[k] --;
	ptr[n].crossings[k] --;
      }
      
      sumekin += ptr[n].v[k]*ptr[n].v[k]*ptr[n].m;

      double ri = (ptr[n].x[k] + ptr[n].cross_neighb[k]*sys->length[k]) 
	- ptr[n].xn[k];
      d2 += sep_Sq(ri);
    }
    
    if ( d2 > sys->max_dist2 ) sys->max_dist2 = d2;
    
    for ( int k=0; k<3 ; k++ )
      for ( int kk=0; kk<3; kk++ )
	retval->kin_P[k][kk] += ptr[n].v[k]*ptr[n].v[kk]*ptr[n].m;
    
  }
  
  if ( sqrt(sys->max_dist2) > sys->skin*0.5 ){
    sys->neighb_flag = 1;
    sys->nupdate_neighb ++;
    sep_set_xn(ptr, sys->npart);
   
    for ( int n=0; n<sys->npart; n++ )
      for ( int k=0; k<3; k++ ) ptr[n].cross_neighb[k] = 0;
  }

  retval->ekin += 0.5*sumekin;
}


void sep_verlet_dpd(seppart *ptr, double lambda, int stepnow,
		    sepsys *sys, sepret *retval){
  double v[3], sumekin=0.0;
  int n;
  const double dt = sys->dt;

  //double max_d2 = 0.0;
  
  for ( n=0; n< sys->npart; n++ ){
    
    for (int k=0; k<3; k++){
      
      ptr[n].a[k] = ptr[n].f[k]/ptr[n].m;

      // Previous time step
      if ( stepnow > 0 )
        ptr[n].v[k] += 0.5*dt*(ptr[n].a[k] + ptr[n].pa[k]);

      ptr[n].x[k] += ptr[n].v[k]*dt + 0.5*dt*dt*ptr[n].a[k];
      ptr[n].pv[k] = ptr[n].v[k] + lambda*dt*ptr[n].a[k];
	      
      ptr[n].pa[k] = ptr[n].a[k];

      /* KE from updated velocities that correspond in time to 
	 positions that were used in the force calculations */
      v[k] = ptr[n].v[k];

      sumekin += v[k]*v[k]*ptr[n].m;

    } 

    double d2 = sep_periodic(ptr, n, sys);
    if ( d2 > sys->max_dist2 ) sys->max_dist2 = d2;
    
    for ( int k=0; k<3 ; k++ )
      for ( int kk=0; kk<3; kk++ )
        retval->kin_P[k][kk] += v[k]*v[kk]*ptr[n].m;

  }  // loop over particles
  
  if ( sqrt(sys->max_dist2) > sys->skin*0.5 ){
    sys->neighb_flag = 1;
    sep_set_xn(ptr, sys->npart);

    for ( int n=0; n<sys->npart; n++ )
      for ( int k=0; k<3; k++ ) ptr[n].cross_neighb[k] = 0;
  }

  retval->ekin += 0.5*sumekin;
}

/**** SHAKE - maybe out? ****/

void sep_set_shake(sepmol *mols, unsigned nuau, 
		   int nb, double *blength, sepsys sys){

  for ( unsigned i=0; i<sys.molptr->num_mols; i++ ){
    if ( mols[i].nuau == nuau ) {
      mols[i].nbonds = nb;
      
      mols[i].blength = sep_vector(nb);
      for ( int n=0; n<nb; n++ )
	mols[i].blength[n] = blength[n];

      mols[i].shake_flag = 1;
    }
  }

}


void sep_shake(sepatom *atoms, sepmol *mols, sepsys *sys, 
	       double tol, sepret *retval){
  
  double dtsq = sys->dt*sys->dt;
  double tol2 = 2.0*tol;
  double wc = 0.0;
  
  double max_d2 = 0.0;

  // LOOP OVER ALL ATOMS NOT IN A MOLECULE
  for ( int n=0; n<sys->npart; n++ ){

    if ( atoms[n].molindex == -1 ){

      for ( int k=0; k<3; k++ ){
	double rw = sep_verlet_wrap(atoms[n].x[k], atoms[n].xp[k], sys->length[k]);
	double a = atoms[n].f[k]/atoms[n].m;
	atoms[n].x[k] = 2.0*rw - atoms[n].xp[k] + dtsq*a;
      }
      
      double d2 = sep_periodic(atoms, n, sys);

      if ( d2 > max_d2 ) max_d2 = d2;
    }

  }      

  // LOOP OVER ALL NON-SHAKE MOLECULES
  for ( unsigned i = 0; i<sys->molptr->num_mols; i++ ){
   
    if ( mols[i].shake_flag == 1 ) continue;

    const unsigned nuau = mols[i].nuau;

    for ( unsigned a=0; a<nuau; a++ ){
      int ia   = mols[i].index[a];
      for ( int k=0; k<3; k++ ){
	double rw = sep_verlet_wrap(atoms[ia].x[k], 
				    atoms[ia].xp[k], sys->length[k]);
	double a = atoms[ia].f[k]/atoms[ia].m;
	atoms[ia].x[k] = 2.0*rw - atoms[ia].xp[k] + dtsq*a;
      }
      
      double d2 = sep_periodic(atoms, ia, sys);

      if ( d2 > max_d2 ) max_d2 = d2;
    }

  }

  // LOOP OVER ALL SHAKE-MOLECULES
  for ( unsigned i = 0; i<sys->molptr->num_mols; i++ ){
   
    if ( mols[i].shake_flag != 1 ) continue;

    const unsigned nuau = mols[i].nuau;

    double *rxi = sep_vector(nuau);
    double *ryi = sep_vector(nuau);
    double *rzi = sep_vector(nuau);
    double *pxi = sep_vector(nuau);
    double *pyi = sep_vector(nuau);
    double *pzi = sep_vector(nuau);
    
    double *mass = sep_vector(nuau); 

    int *moving = sep_vector_int(nuau);
    int *moved  = sep_vector_int(nuau);

    for ( unsigned a=0; a<nuau; a++ ){

      int ia   = mols[i].index[a];
      mass[a] = atoms[ia].m;
      
      rxi[a] = atoms[ia].x[0];
      double xn = sep_verlet_wrap(rxi[a], atoms[ia].xp[0], sys->length[0]);

      ryi[a] = atoms[ia].x[1];
      double yn = sep_verlet_wrap(ryi[a], atoms[ia].xp[1], sys->length[1]);

      rzi[a] = atoms[ia].x[2];
      double zn = sep_verlet_wrap(rzi[a], atoms[ia].xp[2], sys->length[2]);

      pxi[a] = 2.0*xn - atoms[ia].xp[0] + dtsq*atoms[ia].f[0]/mass[a];
      pyi[a] = 2.0*yn - atoms[ia].xp[1] + dtsq*atoms[ia].f[1]/mass[a];
      pzi[a] = 2.0*zn - atoms[ia].xp[2] + dtsq*atoms[ia].f[2]/mass[a];
      
      moving[a] = SEP_FALSE;
      moved[a] = SEP_TRUE;
    }

    int it = 0;
    int done = SEP_FALSE;

    // BEGIN ITERATIVE LOOP
    while ( done==SEP_FALSE && it<SEP_SHAKE_MAXIT ) {
     
      done = SEP_TRUE;
      
      for ( unsigned a=0; a<mols[i].nbonds; a++ ){

	unsigned b = a + 1;
	if ( b == nuau ) b = 0;

	double rabsq = sep_Sq(mols[i].blength[a]); 

	if ( moved[a]==SEP_TRUE || moved[b]==SEP_TRUE ) {

	  double pxab = pxi[a] - pxi[b];
	  sep_Wrap(pxab, sys->length[0]);

	  double pyab = pyi[a] - pyi[b];
	  sep_Wrap(pyab, sys->length[1]);

	  double pzab = pzi[a] - pzi[b];
	  sep_Wrap(pzab, sys->length[2]);

	  double pabsq = pxab*pxab + pyab*pyab + pzab*pzab;
	  double diffsq = rabsq - pabsq;

	  if ( fabs(diffsq) > rabsq*tol2 ){
	    double rxab = rxi[a] - rxi[b];
	    sep_Wrap(rxab, sys->length[0]);

	    double ryab = ryi[a] - ryi[b];
	    sep_Wrap(ryab, sys->length[1]);

	    double rzab = rzi[a] - rzi[b];
	    sep_Wrap(rzab, sys->length[2]);
	    
            double rpab = rxab*pxab + ryab*pyab + rzab*pzab;
	    
	    if ( rpab < rabsq*1.0e-6 )
	      sep_error("Error %d", __LINE__);

	    double rma = 1.0/mass[a];
	    double rmb = 1.0/mass[b];

	    double gab = diffsq/(2.0*(rma+rmb)*rpab);
	    wc += gab*rabsq;

	    double dx = rxab*gab;
	    double dy = ryab*gab;
	    double dz = rzab*gab;
	    
	    pxi[a] += rma*dx;
	    pyi[a] += rma*dy; 
	    pzi[a] += rma*dz; 

	    pxi[b] -= rmb*dx; 
	    pyi[b] -= rmb*dy; 
	    pzi[b] -= rmb*dz; 

	    moving[a] = SEP_TRUE;
	    moving[b] = SEP_TRUE;
	    
	    done = SEP_FALSE;
	  }
	}
      }

      for ( unsigned a=0; a<nuau; a++ ){
	moved[a] = moving[a];
	moving[a] = SEP_FALSE;
      }
    
      it++;
   
    }
    
    if ( done==SEP_FALSE )
     sep_error("Error %d", __LINE__);
  
    for ( unsigned a=0; a<mols[i].nuau; a++ ){
      int ia = mols[i].index[a];
    
      atoms[ia].v[0] = 0.5*(pxi[a] - atoms[ia].xp[0])/sys->dt;
      atoms[ia].v[1] = 0.5*(pyi[a] - atoms[ia].xp[1])/sys->dt;
      atoms[ia].v[2] = 0.5*(pzi[a] - atoms[ia].xp[2])/sys->dt;
   
      atoms[ia].x[0] = pxi[a]; 
      atoms[ia].x[1] = pyi[a]; 
      atoms[ia].x[2] = pzi[a]; 
   
      atoms[ia].xp[0] = rxi[a]; 
      atoms[ia].xp[1] = ryi[a]; 
      atoms[ia].xp[2] = rzi[a]; 
      
      double d2 = sep_periodic(atoms, ia, sys);
      if ( d2 > max_d2 ) max_d2 = d2;
    
      retval->ekin += 0.5*(sep_Sq(atoms[ia].v[0])+sep_Sq(atoms[ia].v[1])+
			   sep_Sq(atoms[ia].v[2]))*atoms[ia].m;

      
    }
   
    free(rxi); free(ryi); free(rzi);
    free(pxi); free(pyi); free(pzi);
    free(mass);
    free(moving); free(moved);

  }
  

  if ( sqrt(max_d2) > sys->skin*0.5 ){
    sys->neighb_flag = 1;
    sys->nupdate_neighb ++;
    sep_set_xn(atoms, sys->npart);
   
    for ( int n=0; n<sys->npart; n++ ){
      for ( int k=0; k<3; k++ ){
	atoms[n].cross_neighb[k] = 0;
      }
    }

  }

  wc /= dtsq/3.0;
}


void sep_set_leapfrog(sepsys *sys, double dt){

  sys->intgr_type = SEP_LEAPFROG;
  sys->dt = dt;

}

