#include "sepewald.h"


double sep_ewald_four(sepatom *atom, int kmax, double lbox, int npart){
  int i, kxi, kyi, kzi, kxn, kyn, kzn, knsq;
  double complex eikx[npart][kmax+1], eiky[npart][2*kmax+1], 
    eikz[npart][2*kmax+1];
  double kx, ky, kz, ksq, preffac, cc, sumeikr, sumU, Ucol;
  const double tpil = 2*SEP_PI/lbox, kappa = KAPPA, //sqrt(1.0/lbox), 
    fac1 = 1./(kappa*kappa*4.), fac2 = 8*SEP_PI/(lbox*lbox*lbox);
 

  // Calculare exp factors and sum Eq. 12.1.18 (Frenkel 2nd edition)
  sumeikr = 0.0;
  Ucol = 0.0;
  for ( kxi=0; kxi<=kmax; kxi++ ){

    kxn = kxi;
    kx = tpil*kxn;
    
    for ( kyi=0; kyi<=2*kmax; kyi++ ){
      
      kyn = kyi-kmax;
      ky = tpil*kyn;
      
      for ( kzi=0; kzi<=2*kmax; kzi++ ){
	
	kzn = kzi-kmax;
	kz = tpil*kzn;

	knsq = kxn*kxn + kyn*kyn + kzn*kzn;
	if ( knsq == 0 ) continue;

	sumU = 0.0;
	for ( i=0; i<npart; i++ ){
	  eikx[i][kxi] = exp(-I*kx*atom[i].x[0]);
	  eiky[i][kyi] = exp(-I*ky*atom[i].x[1]);
	  eikz[i][kzi] = exp(-I*kz*atom[i].x[2]);

	  cc = conj(eikx[i][kxi])*conj(eiky[i][kyi])*conj(eikz[i][kzi]);
	  sumU +=  atom[i].z*cc;  
	  sumeikr +=  atom[i].z*cc; 
	}

	ksq = kx*kx + ky*ky + kz*kz;
	preffac = exp(-ksq*fac1)/ksq; 
	Ucol += preffac*creal(conj(sumU)*sumU);

      }
    }

  }

  // Do the actual force calculation
  for ( i=0; i<npart; i++ ){
    
    for ( kxi=0; kxi<=kmax; kxi++ ){

      kxn = kxi;
      kx = tpil*kxn;
      
      for ( kyi=0; kyi<=2*kmax; kyi++ ){

	kyn = kyi-kmax;
	ky = tpil*kyn;
	
	for ( kzi=0; kzi<=2*kmax; kzi++ ){

	  kzn = kzi-kmax;
	  kz = tpil*kzn;

	  knsq = kxn*kxn + kyn*kyn + kzn*kzn;
	  if ( knsq == 0 ) continue;

	  ksq = kx*kx + ky*ky + kz*kz;
	  preffac = exp(-ksq*fac1)/ksq; 
	  
	  cc = eikx[i][kxi]*eiky[i][kyi]*eikz[i][kzi]*sumeikr;
	  cc = cimag(cc); 

	  atom[i].a[0] += fac2*kx*preffac*atom[i].z*cc/atom[i].m;
	  atom[i].a[1] += fac2*ky*preffac*atom[i].z*cc/atom[i].m;
	  atom[i].a[2] += fac2*kz*preffac*atom[i].z*cc/atom[i].m;

	}

      }
    
    }
  
  }
  
  return 4*SEP_PI*Ucol;

}




double sep_ewald_short(sepatom *ptr, double cf2, double lbox, int npart, 
			  int ndim){
  double Ucol, rv[3], r2, r, ft, f, erfckr, zizj, krsq;
  int n, m, k;
  const double kappa = KAPPA, A = 2*kappa/sqrt(SEP_PI);

  Ucol = 0.0;
  for (n=0; n<npart-1; n++){

    for (m=n+1; m<npart; m++){
    
      r2 = 0.0;
      for (k=0; k<ndim; k++){
        rv[k]  = ptr[n].x[k]-ptr[m].x[k];
	sep_Wrap( rv[k], lbox );
	r2   += rv[k]*rv[k];
      }

      if (r2 < cf2){ 
	r = sqrt(r2);
        krsq = kappa*kappa*r2;
	zizj = ptr[n].z*ptr[m].z;
	erfckr = erfc(kappa*r)/r;
	ft = zizj*(erfckr + A*exp(-krsq))/r2;

	for (k=0; k<ndim; k++){
	  f = ft*rv[k];
	  ptr[n].a[k] += f/ptr[n].m;
	  ptr[m].a[k] -= f/ptr[m].m;
	}

	Ucol += zizj*erfckr;
      }
 
    }
  }

  return Ucol;
}



double sep_ewald_self(sepatom *atom, double lbox, int npart){
  double Ucol, kappa = sqrt(5.0/lbox);
  int i;

  Ucol = 0.0;
  for ( i=0; i<npart; i++ ){
    Ucol -= atom[i].z*atom[i].z;
  }

  return Ucol*kappa/sqrt(SEP_PI);
}


double sep_ewald(sepatom *atom, int kmax, double lbox, int npart){
  double Ucol;

  Ucol = sep_ewald_four(atom, kmax, lbox, npart);
  Ucol = sep_ewald_short(atom, lbox*0.5, lbox, npart, 3);
  Ucol += sep_ewald_self(atom, lbox, npart);

  return Ucol;
}


void sep_3D_ewald_four(sepatom *atom, double kappa, int kmax, sep3D *sys){
  int n, i, kxi, kyi, kzi, kxn, kyn, kzn, knsq;
  double kx, ky, kz, ksq, cc, sumeikr, sumU, Ucol, tpil[3];
  const double fac1 = 1./(kappa*kappa*4.), fac2 = 8*SEP_PI/(sys->volume);
  complex double eikx[sys->npart][kmax+1], 
    eiky[sys->npart][2*kmax+1], 
    eikz[sys->npart][2*kmax+1];
  double pffac[(kmax+1)*(2*kmax + 1)*(2*kmax + 1)];
  
  // Calculate 2*pi/Ly prefactor
  for ( n=0; n<3; n++ ) 
    tpil[n] = 2.0*SEP_PI/sys->length[n];

  sumeikr = 0.0;  Ucol = 0.0; n = 0;

  // Calculare exp factors and sum Eq. 12.1.18 (Frenkel 2nd edition)
  for ( kxi=0; kxi<=kmax; kxi++ ){

    kxn = kxi;
    kx = tpil[0]*kxn;

    for ( kyi=0; kyi<=2*kmax; kyi++ ){
      
      kyn = kyi-kmax;
      ky = tpil[1]*kyn;

      for ( kzi=0; kzi<=2*kmax; kzi++ ){
	
	kzn = kzi-kmax;
	kz =  tpil[2]*kzn;

	knsq = kxn*kxn + kyn*kyn + kzn*kzn;
	if ( knsq == 0 ) continue;

	ksq = kx*kx + ky*ky + kz*kz;
	pffac[n] = exp(-ksq*fac1)/ksq; 

	sumU = 0.0;
	for ( i=0; i<sys->npart; i++ ){
	  eikx[i][kxi] = exp(-I*kx*atom[i].x[0]);
	  eiky[i][kyi] = exp(-I*ky*atom[i].x[1]);
	  eikz[i][kzi] = exp(-I*kz*atom[i].x[2]);

	  cc = conj(eikx[i][kxi])*conj(eiky[i][kyi])*conj(eikz[i][kzi]);
	  sumU +=  atom[i].z*cc;  
	  sumeikr +=  atom[i].z*cc; 
	}

	Ucol += pffac[n]*creal(conj(sumU)*sumU);

	n++;

      }
    
    }

  }
  
  
  // Do the actual force calculation
  for ( i=0; i<sys->npart; i++ ){
    
    n = 0;
    
    for ( kxi=0; kxi<=kmax; kxi++ ){

      kxn = kxi;
      kx =  tpil[0]*kxn;

      for ( kyi=0; kyi<=2*kmax; kyi++ ){

	kyn = kyi-kmax;
	ky = tpil[1]*kyn;

	for ( kzi=0; kzi<=2*kmax; kzi++ ){

	  kzn = kzi-kmax;
	  kz = tpil[2]*kzn;

	  knsq = kxn*kxn + kyn*kyn + kzn*kzn;
	  if ( knsq == 0 ) continue;
	   
	  cc = eikx[i][kxi]*eiky[i][kyi]*eikz[i][kzi]*sumeikr;
	  cc = cimag(cc); 

	  atom[i].a[0] += fac2*kx*pffac[n]*atom[i].z*cc/atom[i].m;
	  atom[i].a[1] += fac2*ky*pffac[n]*atom[i].z*cc/atom[i].m;
	  atom[i].a[2] += fac2*kz*pffac[n]*atom[i].z*cc/atom[i].m;

	  n++;
	}

      }
    
    }
  
  }
  
  sys->retval[0] = 4.0*SEP_PI*Ucol;

}

void sep_3D_ewald_short(sepatom *ptr, double kappa, double cf, sep3D *sys){
  double Ucol, rv[3], r2, r, ft, f, erfckr, zizj, krsq;
  int i1, i2, n, k;
  const double A = 2*kappa/sqrt(SEP_PI), cf2 = cf*cf;


  Ucol = 0.0;
  for (i1 = 0; i1 < sys->npart; i1++){
    n = 0;
    while (1){
      i2 = ptr[i1].neighb[n];
      if (i2 == -1) break;

      r2 = 0.0;
      for (k=0; k<3; k++){
	rv[k] = ptr[i1].x[k] - ptr[i2].x[k];
	if (sys->bound[k] == 'p'){
	  sep_Wrap( rv[k], sys->length[k] );
	}	
	r2 += rv[k]*rv[k];
      }
      
      if (r2 < cf2){ 
	r = sqrt(r2);
        krsq = kappa*kappa*r2;
	zizj = ptr[i1].z*ptr[i2].z;
	erfckr = erfc(kappa*r)/r;
	ft = zizj*(erfckr + A*exp(-krsq))/r2;
	for (k=0; k<3; k++){
	  f = ft*rv[k];
	  ptr[i1].a[k] += f/ptr[i1].m;
	  ptr[i2].a[k] -= f/ptr[i2].m;
	}

	Ucol += zizj*erfckr;
      }
      n++;
    }
  }

  sys->retval[0] = Ucol;
}

void _sep_3D_ewald_short_bond(sepatom *ptr, double kappa, double cf, sep3D *sys){
  double Ucol, rv[3], r2, r, ft, f, erfckr=0, erfkr, zizj, krsq;
  int i1, i2, n, k;
  const double A = 2*kappa/sqrt(SEP_PI), cf2 = cf*cf;


  Ucol = 0.0;
  for (i1 = 0; i1 < sys->npart; i1++){
    n = 0;
    while (1){
      i2 = ptr[i1].neighb[n];
      if (i2 == -1) break;

      r2 = 0.0;
      for (k=0; k<3; k++){
	rv[k] = ptr[i1].x[k] - ptr[i2].x[k];
	if (sys->bound[k] == 'p'){
	  sep_Wrap( rv[k], sys->length[k] );
	}	
	r2 += rv[k]*rv[k];
      }

      r = sqrt(r2);
      krsq = kappa*kappa*r2;
      zizj = ptr[i1].z*ptr[i2].z;
       
      // Not same molecule - short ranged real space contribution
      if ( (ptr[i1].molindex == -1 && r2<cf2) || 
	   (ptr[i1].molindex != ptr[i2].molindex && r2<cf2) ){ 
	erfckr = erfc(kappa*r)/r;
	ft = zizj*(erfckr + A*exp(-krsq))/r2;
	for ( k=0; k<3; k++ ){
	  f = ft*rv[k];
	  ptr[i1].a[k] += f/ptr[i1].m;
	  ptr[i2].a[k] -= f/ptr[i2].m;
	}
	
	Ucol += zizj*erfckr;
      }
      // Same molecule no real space contribution, but remove 
      // the self part from the Fourier part
      else if ( ptr[i1].molindex == ptr[i2].molindex ){
	erfkr = erf(kappa*r)/r;
      	ft = -zizj*(erfkr - A*exp(-krsq))/r2;
	for ( k=0; k<3; k++ ){
	  f = ft*rv[k];
	  ptr[i1].a[k] += f/ptr[i1].m;
	  ptr[i2].a[k] -= f/ptr[i2].m;
	}
	
	Ucol -= zizj*erfckr;
      }

      n++;
    }
  }

  sys->retval[0] = Ucol;
}


 
void sep_3D_ewald_self(sepatom *atom, double kappa, sep3D *sys){
  double Ucol;
  int i;

  Ucol = 0.0;
  for ( i=0; i<sys->npart; i++ ){
    Ucol -= atom[i].z*atom[i].z;
  }

  sys->retval[0] = Ucol*kappa/sqrt(SEP_PI);
}




void sep_3D_ewald(sepatom *atom, double kappa, double cf_short, 
		  int kmax, sep3D *sys){
  double Ucol;
  
  sep_3D_ewald_four(atom, kappa, kmax, sys);
  Ucol = sys->retval[0];
    
  sep_3D_ewald_short(atom, kappa, cf_short, sys);
  //_sep_3D_ewald_short_bond(atom, kappa, cf_short, sys);
  Ucol = sys->retval[0];
  
  //sep_3D_ewald_self(atom, kappa, sys);
  
  sys->retval[0] = Ucol;
  
}

