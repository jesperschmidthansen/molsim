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


#include "sepomp.h"

void sep_omp_make_neighb(seppart *ptr,double cf, sepsys *sys, const unsigned opt){

  if ( cf > sys->cf )
    sep_error("cutoff for an interaction cannot be larger than maximum cutoff");

  if ( sys->neighb_flag == 1 ){        
    
    if ( opt == SEP_ALL )
      sep_neighb(ptr, sys);
    else if ( opt == SEP_NEIGHB_EXCL_BONDED )
      sep_neighb_nonbonded(ptr, sys);
    else if ( opt == SEP_NEIGHB_EXCL_SAME_MOL )
      sep_neighb_excl_same_mol(ptr, sys);
  
    sys->neighb_flag = 0;
  }
  
}

void sep_omp_pairs(double **ftot, const seppart *ptr,
		   const char *types, double cf,
		   double (*fun)(double, char), const sepsys *sys) {
  int i1, i2, n, k;
  double r2, ft, f[3], r[3];
  const double cf2 = cf*cf;

  for (i1=0; i1<sys->npart; i1++){
    
    if ( ptr[i1].type != types[0] && ptr[i1].type != types[1] )
      continue;

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
	  for ( k=0; k<3; k++ ) f[k] = ft*r[k];	  
	 
	  for ( k=0; k<3; k++ ){
	    ftot[i1][k] += f[k];
	    ftot[i2][k] -= f[k];
	  }
	 
	}
      }
      n++;
    }
  }

}


void sep_omp_pairs_lj(double **ftot, const seppart *ptr,
		      const char *types, const double *param, const sepsys *sys) {
  int i1, i2, n, k;
  double r2, ft, f[3], r[3];

  const double cf = param[0];
  const double sigma = param[1];
  const double epsilon = param[2];
 
  const double cf2 = cf*cf;

  
  for (i1=0; i1<sys->npart; i1++){
    
    if ( ptr[i1].type != types[0] && ptr[i1].type != types[1] )
      continue;

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

	  double rri = sigma/r2;
	  double rri3 = rri*rri*rri;

	  ft = 48.0/epsilon *rri3*(rri3 - 0.5)*rri;
	 
	  for ( k=0; k<3; k++ ) f[k] = ft*r[k];	  
	 
	  for ( k=0; k<3; k++ ){
	    ftot[i1][k] += f[k];
	    ftot[i2][k] -= f[k];
	  }
	 
	}
      }
      n++;
    }
  }

}


void sep_omp_coulomb(double **ftot, seppart *ptr, double cf, sepsys *sys) {
  int i1, i2, n;
  double r2, ft, f[3], dr[3];
  
  const double cf2 = cf*cf;
  const double icf2 = 1.0/cf2; 
  
  for (i1=0; i1<sys->npart; i1++){
    n = 0;
    while (1){
      i2 = ptr[i1].neighb[n];
      
      if ( i2 == -1 ) break; 

      r2 = 0.0;
      for ( int k=0; k<3; k++ ){
	dr[k] = ptr[i1].x[k] - ptr[i2].x[k];
	sep_Wrap( dr[k], sys->length[k] );
	r2 += dr[k]*dr[k];
      }
      
      if (r2 < cf2){ 
	
	double zizj = ptr[i2].z*ptr[i1].z;
	double r = sqrt(r2);
	
	ft = zizj*(1.0/r2 - icf2)/r; 
	for ( int k=0; k<3; k++ ) f[k] = ft*dr[k];
	
	for ( int k=0; k<3; k++ ){
	  ftot[i1][k] += f[k];
	  ftot[i2][k] -= f[k];
	}
      } 
      n++;
    }
  }
  
}




void sep_omp_bond(double **ftot, sepatom *aptr, int type, 
		  const double lbond, const double ks, 
		  sepsys *sys){
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
	ftot[a][k] += f[k];
	ftot[b][k] -= f[k];
      }

    }
  }

}

void sep_omp_angle(double **ftot, sepatom *ptr, int type, 
		   const double angle0, const double k, 
		   sepsys *sys){
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
						
	ftot[a][k] += f1;
	ftot[b][k] += (-f1-f2);
	ftot[c][k] += f2;
      }
     
    }
  }

}

// From Rapapport
void sep_omp_torsion(double **ftot, sepatom *ptr, int type, 
		 const double g[6],  sepsys *sys){
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
							
	ftot[a][k] += f1;
	ftot[b][k] += (-(1.0 + cR1)*f1 + cR2*f2);
	ftot[c][k] += (cR1*f1 - (1.0 + cR2)*f2);
	ftot[d][k] += f2;
      }
    }
  }

}



void sep_omp_dpd_pairs(double **f_tot, seppart *ptr, const char *types, 
		       const double cf, const double aij, 
		       const double temp_desired, 
		       const double sigma, sepsys *sys){
  int i1, i2, n, k;
  double r2,  r[3], rhat[3], vij[3], fC[3], fD[3], fR[3], dij, one_dij, randnum;
  const double cf2 = cf*cf, isqrtdt = 1.0/sqrt(sys->dt);
  const double facchk = 2.0*sqrt(3.0);
  const double gamma = sigma*sigma/(2.0*temp_desired);


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
	  one_dij = 1.0-dij;

	  double dotrv = 0.0;
	  for ( int k=0; k<3; k++ ){
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
	    fR[k] = sigma*one_dij*rhat[k]*isqrtdt*randnum; 

	    // Summing up
	    f_tot[i1][k] += fC[k] + fD[k] + fR[k];
	    f_tot[i2][k] -= fC[k] + fD[k] + fR[k];
	  }

	}  //  if r^2 < cutoff
	
      }    //  if atom type is considered
      n++;
    }
  }        //  neighbor list loops
}
