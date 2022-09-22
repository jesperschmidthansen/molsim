#include "sepsampler.h"

// main sampler structure
void sep_add_sampler(sepsampler *sptr, const char *sampler, 
		     sepsys sys, int lvec, ...){
  va_list args;
	 
  va_start(args, lvec); 
  
  if ( strcmp(sampler, "sacf")==0 ){
    if ( sptr->flag_sacf != 1 ){
      double tsample  = va_arg(args, double);
      sptr->sacf = sep_sacf_init(lvec, tsample, sys.dt);
      sptr->flag_sacf = 1;
    }
  }
  else if ( strcmp(sampler, "vacf")==0 ){
    if ( sptr->flag_vacf != 1 ){
      double tsample  = va_arg(args, double);
      sptr->vacf = sep_vacf_init(lvec, tsample, sys);
      sptr->flag_vacf = 1;
    }
  }
  else if ( strcmp(sampler, "msacf")==0 ){
    if ( sptr->flag_msacf != 1 ){
      
      if ( sptr->molptr == NULL ) 
	sep_error("%s at line %d: molpointer not initialized", 
		  __func__, __LINE__);

      double tsample  = va_arg(args, double);
      sptr->msacf = sep_msacf_init(lvec, tsample, sys.dt);
      sptr->flag_msacf = 1;
    }
  }
  else if ( strcmp(sampler, "gh")==0 ){
    if ( sptr->flag_gh != 1 ){
      double tsample  = va_arg(args, double);
      int nwave = va_arg(args, int);

      sptr->gh = sep_gh_init(lvec, tsample, sys.dt, nwave, sys.length[1]);
      sptr->flag_gh = 1;
      
    }
  }
  else if ( strcmp(sampler, "mgh")==0 ){
    if ( sptr->flag_mgh != 1 ){
      if ( sptr->molptr == NULL ) 
	sep_error("%s at line %d: molpointer not initialized", 
		  __func__, __LINE__);

      double tsample  = va_arg(args, double);
      int nwave = va_arg(args, int);
      int safe = va_arg(args, int);
      
      sptr->mgh = sep_mgh_init(lvec, tsample, sys.dt, nwave,
			       sys.length[1], safe); 
      sptr->flag_mgh = 1;
    }
  }
  else if ( strcmp(sampler, "profs")==0 ){
    if ( sptr->flag_profs != 1 ){
      char type = va_arg(args, int);
      int isample = va_arg(args, int);

      sptr->profs = sep_profs_init(type, lvec, isample);
      sptr->flag_profs = 1;
    }
  }
  else if ( strcmp(sampler, "mprofs")==0 ){
    if ( sptr->flag_mprofs != 1 ){

      char type = va_arg(args, int);
      int isample = va_arg(args, int);
      // ACHTUNG!!!!
      int dir = 2;
      int dirvel = 0;  
      int diramom = 1;

      sptr->mprofs = sep_mprofs_init(type, lvec, isample, dir, dirvel, diramom);
      sptr->flag_mprofs = 1;
    }
  }
  else if ( strcmp(sampler, "mcacf")==0 ){
    if ( sptr->flag_mcacf != 1 ){
      
      if ( sptr->molptr == NULL ) 
	sep_error("%s at line %d: molpointer not initialized", 
		  __func__, __LINE__);

      double tsample  = va_arg(args, double);
      sptr->mcacf = sep_mcacf_init(lvec, tsample, sys.dt);
      sptr->flag_mcacf = 1;
    }
  }
  else if ( strcmp(sampler, "mvacf")==0 ){
    if ( sptr->flag_mvacf != 1 ){
      double tsample  = va_arg(args, double);
      sptr->mvacf = sep_mvacf_init(lvec, tsample, sys);
      sptr->flag_mvacf = 1;
    }
  }
  else if ( strcmp(sampler, "mavacf")==0 ){
    if ( sptr->flag_mavacf != 1 ){
      double tsample  = va_arg(args, double);
      sptr->mavacf = sep_mavacf_init(lvec, tsample, sys);
      sptr->flag_mavacf = 1;
    }
  }
  else if ( strcmp(sampler, "radial")==0 ){
    if ( sptr->flag_radial != 1 ){
      int nsample = va_arg(args, int);
      char *t = va_arg(args, char *);
      sptr->radial = sep_radial_init(lvec, nsample, t);
      sptr->flag_radial = 1;
    }
  }
  else if ( strcmp(sampler, "msd")==0 ){
    if ( sptr->flag_msd != 1 ){
      double tsample  = va_arg(args, double);
      int nwave = va_arg(args, int);
      char type = va_arg(args, int);
      sptr->msd = sep_msd_init(lvec, tsample, nwave, type, sys);
      sptr->flag_msd = 1;
    }
  }
  else if ( strcmp(sampler, "mmsd")==0 ){
    if ( sptr->flag_mmsd != 1 ){
      double tsample  = va_arg(args, double);
      int nwave = va_arg(args, int);
      char type = va_arg(args, int);
      sptr->mmsd = sep_mol_msd_init(lvec, tsample, nwave, type, sys);
      sptr->flag_mmsd = 1;
    }
  }
  else {
    sep_error("%s and at line %d: Sampler %s is not recognized", 
	      __func__, __LINE__, sampler);
  }
  
  va_end(args);
}

void sep_add_mol_sampler(sepsampler *sptr, sepmol *mols){
  
  sptr->molptr = mols;

}

void sep_close_sampler(sepsampler *ptr){
  
  if ( ptr->flag_sacf == 1 )    sep_sacf_close(ptr->sacf);
  if ( ptr->flag_vacf == 1 )    sep_vacf_close(ptr->vacf);
  if ( ptr->flag_msacf == 1 )   sep_msacf_close(ptr->msacf);
  if ( ptr->flag_gh == 1 )      sep_gh_close(ptr->gh);
  if ( ptr->flag_mgh == 1 )     sep_mgh_close(ptr->mgh);
  if ( ptr->flag_profs == 1 )   sep_profs_close(ptr->profs);
  if ( ptr->flag_mprofs == 1 )  sep_mprofs_close(ptr->mprofs);
  if ( ptr->flag_mcacf == 1 )   sep_mcacf_close(ptr->mcacf);
  if ( ptr->flag_mvacf == 1 )   sep_mvacf_close(ptr->mvacf);
  if ( ptr->flag_mavacf == 1 )  sep_mavacf_close(ptr->mavacf);
  if ( ptr->flag_radial == 1 )  sep_radial_close(ptr->radial);
  if ( ptr->flag_msd == 1 )     sep_msd_close(ptr->msd);
  if ( ptr->flag_mmsd == 1 )    sep_mol_msd_close(ptr->mmsd);

}

void sep_sample(seppart *pptr, sepsampler *sptr, sepret *ret, sepsys sys, 
		unsigned n){
  
  if ( sptr->flag_sacf==1 && n%sptr->sacf->isample == 0 )
    sep_sacf_sample(sptr->sacf, ret, sys);
  
  if ( sptr->flag_vacf==1 && n%sptr->vacf->isample == 0 )
    sep_vacf_sample(pptr, sptr->vacf, sys);

  if ( sptr->flag_msacf==1 && n%sptr->msacf->isample == 0 )
    sep_msacf_sample(sptr->msacf, pptr, sptr->molptr, ret, sys);

  if ( sptr->flag_gh==1 && n%sptr->gh->isample == 0 )
    sep_gh_sampler(pptr, sptr->gh, sys);

  if ( sptr->flag_mgh==1 && n%sptr->mgh->isample == 0 )
    sep_mgh_sampler(pptr, sptr->molptr, sptr->mgh, sys);

  if ( sptr->flag_profs==1 && n%sptr->profs->isample == 0 )
    sep_profs_sampler(pptr, sptr->profs, sys);

  if ( sptr->flag_mprofs==1 && n%sptr->mprofs->isample == 0 )
    sep_mprofs_sampler(pptr, sptr->molptr, sptr->mprofs, sys);

  if ( sptr->flag_mcacf==1 && n%sptr->mcacf->isample == 0 )
    sep_mcacf_sample(sptr->mcacf, pptr, sptr->molptr, ret, sys);

  if ( sptr->flag_mvacf==1 && n%sptr->mvacf->isample == 0 )
    sep_mvacf_sample(pptr, sptr->mvacf, sptr->molptr, sys);  
  
  if ( sptr->flag_mavacf==1 && n%sptr->mavacf->isample == 0 )
    sep_mavacf_sample(pptr, sptr->mavacf, sptr->molptr, sys);  

  if ( sptr->flag_radial==1 && n%sptr->radial->isample == 0 )
    sep_radial_sample(sptr->radial, pptr, sys);

  if ( sptr->flag_msd==1 ){

    sep_msd_crossing(pptr, sptr->msd, sys);

    if ( sptr->msd->logflag==false && n%sptr->msd->isample == 0 ){
      sep_msd_sample(pptr, sptr->msd, sys);
    }
    else if (  sptr->msd->logflag==true && sptr->msd_counter%sptr->msd->logcounter==0 ){
      sep_msd_sample(pptr, sptr->msd, sys);
      if ( sptr->msd->i == 0 ) {
	sptr->msd->logcounter = 1;
      }
      else {
	sptr->msd->logcounter = 2*sptr->msd->logcounter;
      }
    }

    sptr->msd_counter++;
       
  }
  
  if ( sptr->flag_mmsd==1 && n%sptr->mmsd->isample == 0 )
    sep_mol_msd_sample(pptr, sptr->molptr, sptr->mmsd, sys);

}

sepsampler sep_init_sampler(void){
  
  sepsampler a;

  a.molptr = NULL;

  a.flag_sacf = 0;  a.sacf = NULL;
  a.flag_msacf = 0;  a.msacf = NULL;

  a.flag_vacf = 0;  a.vacf = NULL;
  a.flag_mvacf = 0;  a.mvacf = NULL;
  a.flag_mavacf = 0;  a.mavacf = NULL;

  a.flag_gh = 0; a.gh = NULL;
  a.flag_mgh = 0; a.mgh = NULL;

  a.flag_profs = 0; a.profs = NULL;
  a.flag_mprofs = 0; a.mprofs = NULL;

  a.flag_mcacf = 0; a.mcacf = NULL;

  a.flag_radial = 0; a.radial = NULL;
  
  a.flag_msd = 0; a.msd = NULL;
  a.flag_mmsd = 0; a.mmsd = NULL;

  a.msd_counter = 1;
  return a;
}

/************************************************
 *
 * Function definitions for the individual samplers
 *
 *************************************************/

// Stress ACF
sepsacf *sep_sacf_init(int lvec, double tsample, double dt){
  
  sepsacf *ptr = malloc(sizeof(sepsacf));
  if ( ptr==NULL )
    sep_error("%s at line %d: Couldn't allocate memory", __func__, __LINE__);
  
  ptr->dt   = dt;
  ptr->lvec = lvec;
  ptr->i    = 0; 
  
  ptr->nsample = 0;

  ptr->dtsample = tsample/lvec;
  ptr->isample = (int)(ptr->dtsample/dt);
  if ( ptr->isample < 1 ) 
    sep_error("%s at line %d: isample is too small - CHECK lvec argument",
	      __func__, __LINE__);
  
  ptr->sacf = sep_vector(lvec);
  ptr->stress = sep_matrix(lvec, 3);

  return ptr;
}


void sep_sacf_sample(sepsacf *sptr, sepret *ret, sepsys sys){
  unsigned int k, n, nn;
  

  sep_pressure_tensor(ret, &sys);

  int index = sptr->i;

  sptr->stress[index][0] = -ret->P[0][1];
  sptr->stress[index][1] = -ret->P[0][2];
  sptr->stress[index][2] = -ret->P[1][2];

  (sptr->i)++;

  if ( sptr->i == sptr->lvec ){
    double *parray = sep_vector(sptr->lvec);

#pragma omp parallel for schedule(dynamic) 	\
  private(k, n, nn)				\
  reduction(+:parray[:sptr->lvec])
    for ( k=0; k<3; k++ ){
      for ( n=0; n<sptr->lvec; n++ ){
        for ( nn=0; nn<sptr->lvec-n; nn++ ){
          parray[n] += sptr->stress[nn][k]*sptr->stress[n+nn][k];
        }
      }
    }
    
    for ( n=0; n<sptr->lvec; n++ )  sptr->sacf[n] += parray[n];
    free(parray);

    (sptr->nsample)++;

    FILE *fout1 = fopen("sacf.dat", "w");
    for ( unsigned n=0; n<sptr->lvec; n++ ){
      double t   = n*sptr->dtsample;
      double fac = sys.volume/(3*(sptr->lvec-n)*sptr->nsample);
      fprintf(fout1, "%f %f\n", t, sptr->sacf[n]*fac);
    }

    fclose(fout1); 
    
    sptr->i = 0;
  }
}

void sep_sacf_close(sepsacf *ptr){
  
  free(ptr->sacf);
  sep_free_matrix(ptr->stress, ptr->lvec);
  
  free(ptr);

}

// Atomic radial distriution
sepradial *sep_radial_init(int lvec, int sampleinterval, char types[]){

  sepradial *ptr = malloc(sizeof(sepradial));
  if ( ptr==NULL )
    sep_error("%s at line %d: Couldn't allocate memory", __func__, __LINE__);
  
  ptr->lvec = lvec;
  ptr->isample = sampleinterval;

  // To be changed for multiple neighbour types
  ptr->ntypes = (int)strlen(types);
  for ( int n=0; n<ptr->ntypes; n++ )
    ptr->types[n] = types[n];
    
  ptr->ncomb=0;
  for ( int n=1; n<=ptr->ntypes; n++ ) (ptr->ncomb) += n;

  ptr->hist = sep_matrix_int(lvec, ptr->ncomb);

  FILE *fout=  fopen("radial_info.dat", "w");
  if ( fout == NULL )
    sep_error("%s: Couldn't open file", __func__);

  fprintf(fout, "Pairs in radial.dat columns are\n");
  for ( int nt_1 = 0; nt_1<ptr->ntypes; nt_1++ ){
    for  ( int nt_2 = nt_1; nt_2<ptr->ntypes; nt_2++ ){
      fprintf(fout, "%c%c  ", ptr->types[nt_1],  ptr->types[nt_2]);
    }
  }
     
  fclose(fout);
  
  return ptr;
}

void sep_radial_sample(sepradial *sptr,  seppart *atom, sepsys sys){

  const int npart = sys.npart;
  const double lbox = sys.length[0];
  
  const double dg = 0.5*lbox/sptr->lvec;

  double rv[3];
  
  for ( int i=0; i<npart-1; i++ ){
    for ( int j=i+1; j<npart; j++ ){

      double r2 = 0.0;
      for ( int k=0; k<3; k++){
	rv[k]  = atom[i].x[k]-atom[j].x[k];
	sep_Wrap( rv[k], lbox );
	r2   += rv[k]*rv[k];
      }
	
      double r = sqrt(r2);
      int index = (int)(r/dg);

      if ( index < sptr->lvec ){      
	int counter = 0;
	for ( int nt_1 = 0; nt_1<sptr->ntypes; nt_1++ ){
	  for  ( int nt_2 = nt_1; nt_2<sptr->ntypes; nt_2++ ){
	    if ( (atom[i].type == sptr->types[nt_1] &&
		  atom[j].type == sptr->types[nt_2])
		 ||
		 (atom[i].type == sptr->types[nt_2] &&
		  atom[j].type == sptr->types[nt_1]) ){
	      sptr->hist[index][counter] += 1;
	    }
	    counter++;
	  }
	}
      }
      
    }
  }
  
  sptr->nsample++;

  FILE *fout = fopen("radial.dat", "w");
  if ( fout == NULL )
    sep_error("%s: Couldn't open file", __func__);

  for ( int i=0; i<sptr->lvec; i++ ){
    double vi = pow(i*dg, 3.0);
    double vii = pow((i+1)*dg, 3.0);
    double r = (i+0.5)*dg;
    
    fprintf(fout, "%f ", r);
    
    for ( int n=0; n<sptr->ncomb; n++ ){
      double fac = 1.0/((vii-vi)*(sptr->nsample));
      double g = (double)sptr->hist[i][n]*fac;
    
      fprintf(fout, "%f ", g);
    }
    fprintf(fout, "\n");
  }
  
  fclose(fout);
  
}

void sep_radial_close(sepradial *ptr){

  sep_free_matrix_int(ptr->hist, ptr->lvec);
  
}

sepmsd *sep_msd_init(int lvec, double tsample, int nk, char type, sepsys sys){

  const int npart = sys.npart;
    
  sepmsd *ptr = malloc(sizeof(sepmsd));

  if ( ptr==NULL )
    sep_error("%s at line %d: Couldn't allocate memory", __func__, __LINE__);

  ptr->dt   = sys.dt;
  ptr->i    = 0; 
  ptr->nk   = nk;
  ptr->npart = sys.npart;
  ptr->nsample = 0;
  ptr->type = type;

  
  if (lvec > 0 ){
    ptr->lvec = lvec;
    ptr->dtsample = tsample/lvec;
    ptr->isample = (int)(ptr->dtsample/sys.dt);
    ptr->logflag = false;
    ptr->time  = sep_vector(lvec);
    for ( int n=0; n<lvec; n++ ){
      ptr->time[n] = ptr->dtsample*(n+1);
    }
  }
  else if ( lvec == 0 ){
    double t = sys.dt; lvec = 1;
    while ( t<tsample ){
      lvec++;
      t = 2*t;
    }
    ptr->lvec = lvec;
    ptr->logflag = true;
    ptr->time  = sep_vector(lvec);
    int n=1; int c=0;
    while (c<lvec){
      ptr->time[c] = sys.dt*n;
      n = 2*n; c++;
    }
    ptr->logcounter = 1;
  }
  
  ptr->msd   = sep_vector(lvec);
  ptr->Fs    = sep_complex_matrix(lvec, nk);
  ptr->msdsq = sep_vector(lvec);  

  ptr->crossings = sep_matrix_int(npart, 3);
  ptr->prev_pos = sep_matrix(npart, 3);
  ptr->pos0 = sep_matrix(npart, 3);
  ptr->k  = sep_vector(nk);

  FILE *fout = fopen("msd-k.dat", "w");
  if ( fout==NULL )
    sep_error("%s at line %d: Couldn't open file", __func__, __LINE__);
  
  for ( int n=1; n<=nk; n++ ){
    ptr->k[n-1] = 2*SEP_PI/sys.length[0]*n;
    fprintf(fout, "%f\n", ptr->k[n-1]);
  }
  fclose(fout);
  
  return ptr;
}

void sep_msd_crossing(seppart *atom, sepmsd *sptr, sepsys sys){

  for ( int n=0; n<sys.npart; n++ ){
    for ( int k=0; k<3; k++ ){
      
      if ( sptr->prev_pos[n][k] - atom[n].x[k]  > 0.5*sys.length[k] )
	sptr->crossings[n][k] ++;
      else if ( sptr->prev_pos[n][k] - atom[n].x[k]  < -0.5*sys.length[k] )
	sptr->crossings[n][k] --;

      sptr->prev_pos[n][k] = atom[n].x[k];
    }
  
  }

}

void sep_msd_sample(seppart *atom, sepmsd *sptr, sepsys sys){

  int index = sptr->i;
  if ( index == 0 ){
    for ( int n=0; n<sys.npart; n++ ){
      for ( int k=0; k<3; k++ ) {
	sptr->prev_pos[n][k] = sptr->pos0[n][k] = atom[n].x[k];
	sptr->crossings[n][k] = 0;
      }
    }
  }

  double sd = 0, qd=0.0, dr[3]={0.0};
  for ( int n=0; n<sys.npart; n++ ){
    if ( atom[n].type == sptr->type ){
      for ( int k=0; k<3; k++ ){

	dr[k] = atom[n].x[k] + sptr->crossings[n][k]*sys.length[k] 
	  - sptr->pos0[n][k];
      }

      double a = sep_dot(dr, dr, 3);
      
      sd += a;
      qd += a*a;
    }
  }
  
  for ( int i=0; i<sptr->nk; i++ ){
    for ( int n=0; n<sys.npart; n++ ){
      if ( atom[n].type == sptr->type ){
	dr[0] = atom[n].x[0] + sptr->crossings[n][0]*sys.length[0] 
	  - sptr->pos0[n][0];
	sptr->Fs[index][i] += cexp(I*sptr->k[i]*dr[0]);     
      }
    }
  }

  sptr->msd[index]   += sd;
  sptr->msdsq[index] += qd;
  index ++;
   
  if ( index == sptr->lvec ){
    sptr->nsample++;
    const int ntype = sep_count_type(atom, sptr->type, sptr->npart);

    FILE *fout = fopen("msd.dat", "w");
    if ( fout==NULL )
      sep_error("%s at line %d: Couldn't open file", __func__, __LINE__);
	   
    for ( int n=0; n<sptr->lvec; n++ ){      
      fprintf(fout, "%f %f \n",
	      sptr->time[n], sptr->msd[n]/(ntype*sptr->nsample));
    }
    fclose(fout);
    
    fout = fopen("msd-gaussparam.dat", "w");
    if ( fout==NULL )
      sep_error("%s at line %d: Couldn't open file", __func__, __LINE__);
  
    for ( int n=0; n<sptr->lvec; n++ ){
      double a = sptr->msdsq[n]/(ntype*sptr->nsample);
      double b = sep_Sq(sptr->msd[n]/(ntype*sptr->nsample));
      
      fprintf(fout, "%f %f \n", sptr->time[n], 3.0*a/(5.0*b)-1.);
    }
    fclose(fout);
    

    fout = fopen("msd-incoherent.dat", "w");
    if ( fout==NULL )
      sep_error("%s at line %d: Couldn't open file", __func__, __LINE__);

    for ( int n=0; n<sptr->lvec; n++ ){      
      fprintf(fout, "%f ", sptr->time[n]);
      for ( int i=0; i<sptr->nk; i++ )
	fprintf(fout, "%f ", creal(sptr->Fs[n][i])/(sptr->nsample*ntype));
      fprintf(fout, "\n");
    }
    fclose(fout);
    index = 0;
  }

  sptr->i = index;
}


void sep_msd_close(sepmsd *ptr){

  free(ptr->msd);
  free(ptr->k);
  free(ptr->msdsq);
  free(ptr->time);
   
  sep_free_complex_matrix(ptr->Fs, ptr->lvec);
  sep_free_matrix(ptr->prev_pos, ptr->npart);
  sep_free_matrix(ptr->pos0, ptr->npart);
  sep_free_matrix_int(ptr->crossings, ptr->npart);
  
  free(ptr);

}

// Velocity ACF
sepvacf *sep_vacf_init(int lvec, double tsample, sepsys sys){
  
  sepvacf *ptr = malloc(sizeof(sepvacf));
  if ( ptr==NULL )
    sep_error("%s at line %d: Couldn't allocate memory", __func__, __LINE__);
  
  ptr->dt   = sys.dt;
  ptr->lvec = lvec;
  ptr->i    = 0; 
  
  ptr->nsample = 0;
  
  ptr->dtsample = tsample/lvec;
  ptr->isample = (int)(ptr->dtsample/sys.dt);
  if ( ptr->isample < 1 ) 
    sep_error("%s at line %d: isample is too small - CHECK lvec argument",
	      __func__, __LINE__);

  ptr->vacf = sep_vector(lvec);
  ptr->vels = sep_matrix(lvec, sys.npart);
  
  return ptr;
}

 
void sep_vacf_sample(seppart *ptr, sepvacf *vptr, sepsys sys){
  
  int index = vptr->i;
  for ( int i=0; i<sys.npart; i++ )
    vptr->vels[index][i] = ptr[i].v[0];
  
  (vptr->i)++;
  
  if ( vptr->i == vptr->lvec ){
    
    for ( int i=0; i<sys.npart; i++ )
      for ( unsigned k=0; k<vptr->lvec; k++ )
        for ( unsigned kk=0; kk<vptr->lvec-k; kk++ )
          vptr->vacf[k] += vptr->vels[kk][i]*vptr->vels[k+kk][i];   
  
    (vptr->nsample)++;

    FILE *fout = fopen("vacf.dat", "w");
    for ( unsigned k=0; k<vptr->lvec; k++ ){
      double t   = k*vptr->dtsample;
      double fac = 1.0/((vptr->lvec-k)*sys.npart*vptr->nsample);
      fprintf(fout, "%f %f\n", t, vptr->vacf[k]*fac);
    }
    fclose(fout);

    vptr->i=0;
    
  }

}

void sep_vacf_close(sepvacf *ptr){

  free(ptr->vacf);
  sep_free_matrix(ptr->vels, ptr->lvec);
  
  free(ptr);

}


// Molecular stress ACF
sepmsacf *sep_msacf_init(int lvec, double tsample, double dt){
  
  sepmsacf *ptr = malloc(sizeof(sepmsacf));
  if ( ptr==NULL )
    sep_error("%s at line %d: Couldn't allocate memory", __func__, __LINE__);
  
  ptr->dt   = dt;
  ptr->lvec = lvec;
  ptr->i    = 0; 
  
  ptr->nsample = 0;

  ptr->dtsample = tsample/lvec;
  ptr->isample = (int)(ptr->dtsample/dt);
  if ( ptr->isample < 1 ) 
    sep_error("%s at line %d: isample is too small - CHECK lvec argument",
	      __func__, __LINE__);


  ptr->sacf = sep_matrix(lvec, 2);
  ptr->sstress = sep_matrix(lvec, 3);
  ptr->astress = sep_matrix(lvec, 3);

  return ptr;
}


void sep_msacf_sample(sepmsacf *sptr, sepatom *atoms, sepmol *mols, 
		      sepret *ret, sepsys sys){
  unsigned k, n, nn;
  
  sep_mol_pressure_tensor(atoms, mols, ret, &sys);

  int index = sptr->i;

  sptr->sstress[index][0] = 0.5*(ret->P_mol[0][1]+ret->P_mol[1][0]);
  sptr->sstress[index][1] = 0.5*(ret->P_mol[0][2]+ret->P_mol[2][0]);
  sptr->sstress[index][2] = 0.5*(ret->P_mol[1][2]+ret->P_mol[2][1]);

  sptr->astress[index][0] = 0.5*(ret->P_mol[0][1]-ret->P_mol[1][0]);
  sptr->astress[index][1] = 0.5*(ret->P_mol[0][2]-ret->P_mol[2][0]);
  sptr->astress[index][2] = 0.5*(ret->P_mol[1][2]-ret->P_mol[2][1]);

  (sptr->i)++;

  if ( sptr->i == sptr->lvec ){
    double *parray_sym = sep_vector(sptr->lvec);
    double *parray_asym = sep_vector(sptr->lvec);

#pragma omp parallel for schedule(dynamic) 	\
  private(k, n, nn)				\
  reduction(+:parray_sym[:sptr->lvec],parray_asym[:sptr->lvec])
    for ( k=0; k<3; k++ ){
      for ( n=0; n<sptr->lvec; n++ ){
        for ( nn=0; nn<sptr->lvec-n; nn++ ){
          parray_sym[n] += sptr->sstress[nn][k]*sptr->sstress[n+nn][k];
	  parray_asym[n] += sptr->astress[nn][k]*sptr->astress[n+nn][k];
        }
      }
    }

    for ( unsigned n=0; n<sptr->lvec; n++ ){ 
      sptr->sacf[n][0] += parray_sym[n];
      sptr->sacf[n][1] += parray_asym[n];
    }
    free(parray_sym); free(parray_asym);
    
    (sptr->nsample)++;

    FILE *fout1 = fopen("msacf.dat", "w");
   
    for ( unsigned n=0; n<sptr->lvec; n++ ){
      double t   = n*sptr->dtsample;
      double fac = sys.volume/(3*(sptr->lvec-n)*sptr->nsample);
      fprintf(fout1, "%f %f %f\n", t, sptr->sacf[n][0]*fac, 
	      sptr->sacf[n][1]*fac);
    }

    fclose(fout1); 
    
    sptr->i = 0;
    
  }

}

void sep_msacf_close(sepmsacf *ptr){
  
 
  sep_free_matrix(ptr->sacf, ptr->lvec);
  sep_free_matrix(ptr->sstress, ptr->lvec);
  sep_free_matrix(ptr->astress, ptr->lvec);
  
  free(ptr);

}


// Gen. hydro sampler (atomic)
sepgh *sep_gh_init(int lvec, double tsample, double dt,
		   int nwave, double Ldir){
  
  sepgh *ptr = malloc(sizeof(sepgh));
  if ( ptr==NULL )
    sep_error("%s at line %d: Couldn't allocate memory", __func__, __LINE__);
  
  ptr->dt   = dt;
  ptr->lvec = lvec;
  ptr->i    = 0; 
  
  ptr->nsample = 0;

  ptr->ncalls = 0;
  ptr->avekin = 0.0;
  
  ptr->dtsample = tsample/lvec;
  ptr->isample = (int)(ptr->dtsample/dt);
  if ( ptr->isample < 1 ) 
    sep_error("%s at line %d: isample is too small - CHECK lvec argument",
	      __func__, __LINE__);

  ptr->kdir = 1; // Wavevector k=(0,k,0)
  ptr->tdir = 0; // Tranverse direction is x 

  ptr->nwave = nwave;
  ptr->k = sep_vector(nwave);
  FILE *fout = fopen("gh-wavevector.dat", "w");
  if ( fout == NULL )
    sep_error("%s at %d: Couldn't open file k.dat");

  for ( int n=1; n<=nwave; n++ ){
    ptr->k[n-1] = 2*SEP_PI*n/Ldir;
    fprintf(fout, "%f\n", ptr->k[n-1]);
  }

  fclose(fout);

  ptr->fk_tv  = sep_complex_matrix(lvec, nwave);
  ptr->fkm_tv = sep_complex_matrix(lvec, nwave);

  ptr->fk_lv  = sep_complex_matrix(lvec, nwave);
  ptr->fkm_lv = sep_complex_matrix(lvec, nwave);

  ptr->fk_rho  = sep_complex_matrix(lvec, nwave);
  ptr->fkm_rho = sep_complex_matrix(lvec, nwave);
  
  ptr->fk_e  = sep_complex_matrix(lvec, nwave);
  ptr->fkm_e = sep_complex_matrix(lvec, nwave);
  
  ptr->fk_X  = sep_complex_matrix(lvec, nwave);
  ptr->fkm_X = sep_complex_matrix(lvec, nwave);

  ptr->c_tv     = sep_complex_matrix(lvec, nwave);
  ptr->c_lv     = sep_complex_matrix(lvec, nwave);
  ptr->c_scatt  = sep_complex_matrix(lvec, nwave);
  ptr->c_e      = sep_complex_matrix(lvec, nwave);
  ptr->c_re     = sep_complex_matrix(lvec, nwave);
  ptr->c_rj     = sep_complex_matrix(lvec, nwave);
  ptr->c_ej     = sep_complex_matrix(lvec, nwave);
  ptr->c_er     = sep_complex_matrix(lvec, nwave);
  ptr->c_jr     = sep_complex_matrix(lvec, nwave);
  ptr->c_je     = sep_complex_matrix(lvec, nwave);
  
  ptr->c_X      = sep_complex_matrix(lvec, nwave);

  return ptr;
}


void sep_gh_close(sepgh *ptr){
  int lvec = ptr->lvec;

  sep_free_complex_matrix(ptr->fk_tv, lvec);
  sep_free_complex_matrix(ptr->fkm_tv, lvec);
  sep_free_complex_matrix(ptr->fk_lv, lvec);
  sep_free_complex_matrix(ptr->fkm_lv, lvec);
  sep_free_complex_matrix(ptr->fk_rho, lvec);
  sep_free_complex_matrix(ptr->fkm_rho, lvec);
  sep_free_complex_matrix(ptr->fk_e, lvec);
  sep_free_complex_matrix(ptr->fkm_e, lvec);

 
  sep_free_complex_matrix(ptr->fk_X, lvec);
  sep_free_complex_matrix(ptr->fkm_X, lvec);

  sep_free_complex_matrix(ptr->c_lv, lvec);
  sep_free_complex_matrix(ptr->c_scatt, lvec);
  sep_free_complex_matrix(ptr->c_tv, lvec);
  sep_free_complex_matrix(ptr->c_e, lvec);
  sep_free_complex_matrix(ptr->c_re, lvec);
  sep_free_complex_matrix(ptr->c_rj, lvec);
  sep_free_complex_matrix(ptr->c_ej, lvec);
  sep_free_complex_matrix(ptr->c_er, lvec);
  sep_free_complex_matrix(ptr->c_jr, lvec);
  sep_free_complex_matrix(ptr->c_je, lvec);

  sep_free_complex_matrix(ptr->c_X, lvec);
   
 }


void sep_gh_sampler(sepatom *ptr, sepgh *gh, sepsys sys){

  int index = gh->i;

  sep_eval_xtrue(ptr, &sys);
  unsigned kdir = gh->kdir;
  unsigned tdir = gh->tdir;

  // Calculate the average kinetic energy
  double sumv2 = 0.0;
  for ( int m=0; m<sys.npart; m++ ){
    for ( int kk=0; kk<3; kk++)
      sumv2 += 0.5*ptr[m].m*sep_Sq(ptr[m].v[kk]);
  }

  (gh->ncalls)++ ;
  gh->avekin += sumv2/(sys.npart*gh->ncalls);
  
  for ( unsigned n=0; n<gh->nwave; n++ ){
    
    gh->fk_tv[index][n]    = gh->fkm_tv[index][n]    = 0.0 + 0.0*I;
    gh->fk_lv[index][n]    = gh->fkm_lv[index][n]    = 0.0 + 0.0*I;
    gh->fk_rho[index][n]   = gh->fkm_rho[index][n]   = 0.0 + 0.0*I;
    gh->fk_e[index][n]     = gh->fkm_e[index][n]     = 0.0 + 0.0*I;
    
    gh->fk_X[index][n] = gh->fkm_X[index][n] = 0.0 + 0.0*I;
    
    for ( int m=0; m<sys.npart; m++ ){
      double mass = ptr[m].m;

      complex double kfac   = cexp(I*gh->k[n]*ptr[m].xtrue[kdir]);
      complex double kfac_m = cexp(-I*gh->k[n]*ptr[m].xtrue[kdir]);

      gh->fk_rho[index][n]   += mass*kfac;
      gh->fkm_rho[index][n]  += mass*kfac_m;
     
      gh->fk_tv[index][n]   +=  mass*ptr[m].v[tdir]*kfac;
      gh->fkm_tv[index][n]  +=  mass*ptr[m].v[tdir]*kfac_m;
         
      gh->fk_lv[index][n]   +=  mass*ptr[m].v[kdir]*kfac;
      gh->fkm_lv[index][n]  +=  mass*ptr[m].v[kdir]*kfac_m;
      
      double ekin = 0.0;
      for ( int kk=0; kk<3; kk++) ekin += 0.5*mass*sep_Sq(ptr[m].v[kk]); 
      gh->fk_e[index][n]   +=  (ekin - gh->avekin)*kfac;
      gh->fkm_e[index][n]  +=  (ekin - gh->avekin)*kfac_m;

      // Auxillary quantity
      gh->fk_X[index][n]   +=  mass*(ptr[m].a[1]-I*gh->k[n]*sep_Sq(ptr[m].v[1]))*kfac;
      gh->fkm_X[index][n]  +=  mass*(ptr[m].a[1]-I*gh->k[n]*sep_Sq(ptr[m].v[1]))*kfac_m;

    }
  }

  (gh->i)++;
  
  if ( gh->i == gh->lvec ){

    for ( unsigned k=0; k<gh->nwave; k++ ){
      for ( unsigned n=0; n<gh->lvec; n++ ){
	for ( unsigned nn=0; nn<gh->lvec-n; nn++ ){
	  // Auto-correlations
	  gh->c_tv[n][k]    += gh->fkm_tv[nn][k]*gh->fk_tv[n+nn][k];
	  gh->c_lv[n][k]    += gh->fkm_lv[nn][k]*gh->fk_lv[n+nn][k];
	  gh->c_scatt[n][k] += gh->fk_rho[nn][k]*gh->fkm_rho[n+nn][k];
	  gh->c_e[n][k]     += gh->fk_e[nn][k]*gh->fkm_e[n+nn][k];

	  // Cross-correlations
	  gh->c_re[n][k]     += gh->fk_rho[nn][k]*gh->fkm_e[n+nn][k];
	  gh->c_rj[n][k]     += gh->fk_rho[nn][k]*gh->fkm_lv[n+nn][k];
	  gh->c_ej[n][k]     += gh->fk_e[nn][k]*gh->fkm_lv[n+nn][k];
	  gh->c_er[n][k]     += gh->fk_rho[nn+n][k]*gh->fkm_e[nn][k];
	  gh->c_jr[n][k]     += gh->fk_rho[nn+n][k]*gh->fkm_lv[nn][k];
	  gh->c_je[n][k]     += gh->fk_e[nn+n][k]*gh->fkm_lv[nn][k];

	  // Auxillary-correlation function
	  gh->c_X[n][k] += gh->fk_X[nn][k]*gh->fkm_lv[n+nn][k];
	}
      }
    }

    (gh->nsample)++;

    FILE *fout_scatt = fopen("gh-rho-acf.dat", "w");
    FILE *fout_tv = fopen("gh-trans-momentum-acf.dat", "w");
    FILE *fout_lv = fopen("gh-long-momentum-acf.dat", "w");
    FILE *fout_e = fopen("gh-energy-acf.dat", "w");

    FILE *fout_de = fopen("gh-rho-energy-ccf.dat", "w");
    FILE *fout_dlv = fopen("gh-rho-long-momentum-ccf.dat", "w");
    FILE *fout_elv = fopen("gh-energy-long-momentum-ccf.dat", "w");

    FILE *fout_ed = fopen("gh-energy-rho-ccf.dat", "w");
    FILE *fout_lvd = fopen("gh-long-momentum-rho-ccf.dat", "w");
    FILE *fout_lve = fopen("gh-long-momentum-energy-ccf.dat", "w");
   
    FILE *fout_X = fopen("gh-X-acf.dat", "w");
    
    if ( fout_tv==NULL || fout_lv==NULL || fout_scatt==NULL || fout_e==NULL ||
	 fout_de==NULL || fout_dlv==NULL || fout_elv==NULL ||
	 fout_ed==NULL || fout_lvd==NULL || fout_lve==NULL || fout_X==NULL )
      sep_error("%s at line %d: Error opening files", __func__, __LINE__);
    
    for ( unsigned n=0; n<gh->lvec; n++ ){
      double t   = n*gh->dtsample;
      
      fprintf(fout_tv, "%f ", t);
      fprintf(fout_lv, "%f ", t);
      fprintf(fout_scatt, "%f ", t);
      fprintf(fout_e, "%f ", t);

      fprintf(fout_de, "%f ", t);
      fprintf(fout_dlv, "%f ", t);
      fprintf(fout_elv, "%f ", t);

      fprintf(fout_ed, "%f ", t);
      fprintf(fout_lvd, "%f ", t);
      fprintf(fout_lve, "%f ", t);

      fprintf(fout_X, "%f ", t);
         
      double fac = 1.0/(gh->nsample*sys.volume*(gh->lvec-n));
      
      for ( unsigned k=0; k<gh->nwave; k++ ){
	fprintf(fout_tv, "%f %f ", creal(gh->c_tv[n][k]*fac),
		cimag(gh->c_tv[n][k]*fac));
	fprintf(fout_lv, "%f %f ", creal(gh->c_lv[n][k]*fac),
		cimag(gh->c_lv[n][k]*fac));
	fprintf(fout_scatt, "%f %f ", creal(gh->c_scatt[n][k]*fac),
		cimag(gh->c_scatt[n][k]*fac)) ;
	fprintf(fout_e, "%f %f ", creal(gh->c_e[n][k]*fac),
		cimag(gh->c_e[n][k]*fac));

	fprintf(fout_de, "%f %f ", creal(gh->c_re[n][k]*fac), 
		cimag(gh->c_re[n][k]*fac));
	fprintf(fout_dlv, "%f %f ", creal(gh->c_rj[n][k]*fac), 
		cimag(gh->c_rj[n][k]*fac));
	fprintf(fout_elv, "%f %f ", creal(gh->c_ej[n][k]*fac), 
		cimag(gh->c_ej[n][k]*fac));

	fprintf(fout_ed, "%f %f ", creal(gh->c_er[n][k]*fac), 
		cimag(gh->c_er[n][k]*fac));
	fprintf(fout_lvd, "%f %f ", creal(gh->c_jr[n][k]*fac), 
		cimag(gh->c_jr[n][k]*fac));
	fprintf(fout_lve, "%f %f ", creal(gh->c_je[n][k]*fac), 
		cimag(gh->c_je[n][k]*fac));

	fprintf(fout_X, "%f %f ", creal(gh->c_X[n][k]*fac),
		cimag(gh->c_X[n][k]*fac));
      }
      
      fprintf(fout_tv, "\n"); 
      fprintf(fout_lv, "\n"); 
      fprintf(fout_scatt, "\n");
      fprintf(fout_e, "\n");

      fprintf(fout_de, "\n");
      fprintf(fout_dlv, "\n");
      fprintf(fout_elv, "\n");

      fprintf(fout_ed, "\n");
      fprintf(fout_lvd, "\n");
      fprintf(fout_lve, "\n");
      
      
      fprintf(fout_X, "\n");

    }
	
    fclose(fout_tv); fclose(fout_lv); fclose(fout_scatt); fclose(fout_e);

    fclose(fout_de);    fclose(fout_dlv);    fclose(fout_elv);
    fclose(fout_ed);    fclose(fout_lvd);    fclose(fout_lve);

    fclose(fout_X);
      
    gh->i = 0;
  }

  
}


// Molecular gen. hydrodynamics sampler
sepmgh *sep_mgh_init(int lvec, double tsample, double dt, int nwave, 
		     double Ly, int safe){

  
  sepmgh *ptr = malloc(sizeof(sepmgh));
  if ( ptr==NULL )
    sep_error("%s at line %d: Couldn't allocate memory", __func__, __LINE__);
  
  ptr->dt   = dt;
  ptr->lvec = lvec;
  ptr->i    = 0; 
  ptr->ncalls = 0;
  ptr->nsample = 0;

  ptr->dtsample = tsample/lvec;
  ptr->isample = (int)(ptr->dtsample/dt);
  if ( ptr->isample < 1 ) 
    sep_error("%s at line %d: isample is too small - CHECK lvec argument",
	      __func__, __LINE__);

  ptr->nwave = nwave;
  ptr->k = sep_vector(nwave);
  FILE *fout = fopen("mgh-wavevector.dat", "w");
  if ( fout == NULL )
    sep_error("%s at %d: Couldn't open file k.dat");

  for ( int n=1; n<=nwave; n++ ){
    ptr->k[n-1] = 2*SEP_PI*n/Ly;
    fprintf(fout, "%f\n", ptr->k[n-1]);
  }

  fclose(fout);

  if ( safe==true )
    ptr->safe = true;
  else {
    sep_warning("ACHTUNG - unsafe mode for mgh sampler. Assuming uniaxial single component system");
    ptr->safe = false;
  }
 
  ptr->fk_tv  = sep_complex_matrix(lvec, nwave);
  ptr->fkm_tv = sep_complex_matrix(lvec, nwave);

  ptr->fk_lv  = sep_complex_matrix(lvec, nwave);
  ptr->fkm_lv = sep_complex_matrix(lvec, nwave);

  ptr->fk_rho  = sep_complex_matrix(lvec, nwave);
  ptr->fkm_rho = sep_complex_matrix(lvec, nwave);

  ptr->fk_tav  = sep_complex_matrix(lvec, nwave);
  ptr->fkm_tav = sep_complex_matrix(lvec, nwave);

  ptr->fk_lav  = sep_complex_matrix(lvec, nwave);
  ptr->fkm_lav = sep_complex_matrix(lvec, nwave);

  ptr->fk_e    = sep_complex_matrix(lvec, nwave);
  ptr->fkm_e   = sep_complex_matrix(lvec, nwave);

  ptr->fk_vav  =  sep_complex_matrix(lvec, nwave);
  ptr->fkm_vav =  sep_complex_matrix(lvec, nwave);

  ptr->fk_dip  =  sep_complex_matrix(lvec, nwave);
  ptr->fkm_dip =  sep_complex_matrix(lvec, nwave);

  ptr->c_tv     = sep_complex_matrix(lvec, nwave);
  ptr->c_lv     = sep_complex_matrix(lvec, nwave);
  ptr->c_scatt  = sep_complex_matrix(lvec, nwave);
  ptr->c_tav    = sep_complex_matrix(lvec, nwave);
  ptr->c_lav    = sep_complex_matrix(lvec, nwave);
  ptr->c_e      = sep_complex_matrix(lvec, nwave);
  ptr->c_vav    = sep_complex_matrix(lvec, nwave);
  ptr->c_dip    = sep_complex_matrix(lvec, nwave);

  
  // Template
  ptr->fk_X  = sep_complex_matrix(lvec, nwave);
  ptr->fkm_X = sep_complex_matrix(lvec, nwave);

  ptr->c_X   = sep_complex_matrix(lvec, nwave);

  return ptr;
}

void sep_mgh_close(sepmgh *ptr){
  int lvec = ptr->lvec;

  sep_free_complex_matrix(ptr->fk_tv, lvec);
  sep_free_complex_matrix(ptr->fkm_tv, lvec);
  sep_free_complex_matrix(ptr->fk_lv, lvec);
  sep_free_complex_matrix(ptr->fkm_lv, lvec);
  sep_free_complex_matrix(ptr->fk_rho, lvec);
  sep_free_complex_matrix(ptr->fkm_rho, lvec);
  sep_free_complex_matrix(ptr->fk_tav, lvec);
  sep_free_complex_matrix(ptr->fkm_tav, lvec);
  sep_free_complex_matrix(ptr->fk_lav, lvec);
  sep_free_complex_matrix(ptr->fkm_lav, lvec);
  sep_free_complex_matrix(ptr->fk_e, lvec);
  sep_free_complex_matrix(ptr->fkm_e, lvec);
  sep_free_complex_matrix(ptr->fk_vav, lvec);
  sep_free_complex_matrix(ptr->fkm_vav, lvec);
  sep_free_complex_matrix(ptr->fk_dip, lvec);
  sep_free_complex_matrix(ptr->fkm_dip, lvec);

  
  sep_free_complex_matrix(ptr->c_lv, lvec);
  sep_free_complex_matrix(ptr->c_scatt, lvec);
  sep_free_complex_matrix(ptr->c_tv, lvec);
  sep_free_complex_matrix(ptr->c_tav, lvec);
  sep_free_complex_matrix(ptr->c_lav, lvec);
  sep_free_complex_matrix(ptr->c_e, lvec);
  sep_free_complex_matrix(ptr->c_vav, lvec);
  sep_free_complex_matrix(ptr->c_dip, lvec);
  
  // The template
  sep_free_complex_matrix(ptr->fk_X, lvec);
  sep_free_complex_matrix(ptr->fkm_X, lvec);

  sep_free_complex_matrix(ptr->c_X, lvec);
  
}

void sep_mgh_sampler(sepatom *ptr, sepmol *mols, sepmgh *mgh, sepsys sys){

  const int num_mols = sys.molptr->num_mols;
  int index = mgh->i;

  // Evaluate the center-of-mass properties of each molecule
  sep_mol_cm(ptr, mols, &sys);
  sep_mol_velcm(ptr, mols, &sys);
  sep_mol_eval_xtrue(ptr, mols, sys);  

  sep_mol_spin(ptr, mols, &sys, mgh->safe); // ACHTUNG HERE!!!!
  sep_mol_dipoles(ptr, mols, &sys);

  double tot_mass = 0.0;
  double sumv2 = 0.0;
  for ( int m=0; m<num_mols; m++ ){
    tot_mass += mols[m].m;
    for ( int kk=0; kk<3; kk++)
      sumv2 += 0.5*mols[m].m*sep_Sq(mols[m].v[kk]);
  }

  (mgh->ncalls)++ ;
  mgh->avekin += sumv2/(num_mols*mgh->ncalls);


  for ( unsigned n=0; n<mgh->nwave; n++ ){ 
    
    mgh->fk_tv[index][n]  = mgh->fkm_tv[index][n]  = 0.0 + 0.0*I;
    mgh->fk_lv[index][n]  = mgh->fkm_lv[index][n]  = 0.0 + 0.0*I;
    mgh->fk_rho[index][n] = mgh->fkm_rho[index][n] = 0.0 + 0.0*I;
    mgh->fk_tav[index][n] = mgh->fkm_tav[index][n] = 0.0 + 0.0*I;
    mgh->fk_lav[index][n] = mgh->fkm_lav[index][n] = 0.0 + 0.0*I;
    mgh->fk_vav[index][n] = mgh->fkm_vav[index][n] = 0.0 + 0.0*I;
    mgh->fk_dip[index][n] = mgh->fkm_dip[index][n] = 0.0 + 0.0*I;

    mgh->fk_X[index][n] = mgh->fkm_X[index][n] = 0.0 + 0.0*I;
  
    for ( unsigned m=0; m<sys.molptr->num_mols; m++ ){
      const double mass = mols[m].m;

      complex double kfac   = cexp(I*mgh->k[n]*mols[m].x[1]);
      complex double kfac_m = cexp(-I*mgh->k[n]*mols[m].x[1]);
  
      mgh->fk_rho[index][n]   +=  mass*kfac;
      mgh->fkm_rho[index][n]  +=  mass*kfac_m;

      mgh->fk_tv[index][n]   +=  mass*mols[m].v[0]*kfac;
      mgh->fkm_tv[index][n]  +=  mass*mols[m].v[0]*kfac_m;

      mgh->fk_lv[index][n]   +=  mass*mols[m].v[1]*kfac;
      mgh->fkm_lv[index][n]  +=  mass*mols[m].v[1]*kfac_m;

      double ekin = 0.0;
      for ( int kk=0; kk<3; kk++) ekin += 0.5*mass*sep_Sq(mols[m].v[kk]);
      mgh->fk_e[index][n]  += (ekin - mgh->avekin)*kfac;
      mgh->fkm_e[index][n] += (ekin - mgh->avekin)*kfac_m;


      mgh->fk_tav[index][n]   +=  mols[m].s[2]*kfac;
      mgh->fkm_tav[index][n]  +=  mols[m].s[2]*kfac_m;
     
      mgh->fk_lav[index][n]   +=  mols[m].s[1]*kfac;
      mgh->fkm_lav[index][n]  +=  mols[m].s[1]*kfac_m;

      mgh->fk_vav[index][n]   +=  mols[m].s[2]*kfac;
      mgh->fkm_vav[index][n]  +=  mass*mols[m].v[0]*kfac_m;

      mgh->fk_dip[index][n]   +=  mols[m].pel[1]*kfac;
      mgh->fkm_dip[index][n]  +=  mols[m].pel[1]*kfac_m;

      if ( mgh->safe == false ){
	mgh->fk_X[index][n]  +=  mols[m].w[0]*kfac;
	mgh->fkm_X[index][n] +=  mols[m].w[0]*kfac_m;
      }
      
    }
  }

  (mgh->i)++;

  if ( mgh->i == mgh->lvec ){
    for ( unsigned k=0; k<mgh->nwave; k++ ){
      for ( unsigned n=0; n<mgh->lvec; n++ ){
	for ( unsigned nn=0; nn<mgh->lvec-n; nn++ ){
	  
	  mgh->c_tv[n][k]    += mgh->fkm_tv[nn][k]*mgh->fk_tv[n+nn][k];
	  mgh->c_lv[n][k]    += mgh->fkm_lv[nn][k]*mgh->fk_lv[n+nn][k];
	  mgh->c_scatt[n][k] += mgh->fk_rho[nn][k]*mgh->fkm_rho[n+nn][k];
	  mgh->c_e[n][k]     += mgh->fk_e[nn][k]*mgh->fkm_e[n+nn][k];

	  mgh->c_tav[n][k]   += mgh->fkm_tav[nn][k]*mgh->fk_tav[n+nn][k];
	  mgh->c_lav[n][k]   += mgh->fkm_lav[nn][k]*mgh->fk_lav[n+nn][k];
	 
	  mgh->c_vav[n][k]   += mgh->fkm_vav[nn][k]*mgh->fk_vav[n+nn][k];

	  mgh->c_dip[n][k]   += mgh->fkm_dip[nn][k]*mgh->fk_dip[n+nn][k];

	  if ( mgh->safe == false )
	    mgh->c_X[n][k]   += mgh->fkm_X[nn][k]*mgh->fk_X[n+nn][k];
	}
      }
    }

    (mgh->nsample)++;

    FILE *fout_tv = fopen("mgh-trans-momentum-acf.dat", "w");
    FILE *fout_lv = fopen("mgh-long-momentum-acf.dat", "w");
    FILE *fout_scatt = fopen("mgh-rho-acf.dat", "w");
    FILE *fout_tav = fopen("mgh-trans-angmomentum-acf.dat", "w");
    FILE *fout_lav = fopen("mgh-long-angmomentum-acf.dat", "w");
    FILE *fout_e = fopen("mgh-energy-acf.dat", "w");
    FILE *fout_vav = fopen("mgh-momentum-angmomentum-ccf.dat", "w");
    FILE *fout_dip = fopen("mgh-dipole-acf.dat", "w");
    FILE *fout_X = fopen("mgh-X-cf.dat", "w");

    if ( fout_tv==NULL || fout_lv==NULL || fout_scatt==NULL ||
	 fout_tav==NULL || fout_lav==NULL || fout_e==NULL || fout_vav==NULL
	 || fout_X==NULL || fout_dip==NULL )
        sep_error("%s at line %d: Error opening files", __func__, __LINE__);

    
    for ( unsigned n=0; n<mgh->lvec; n++ ){
      double t   = n*mgh->dtsample;
      
      fprintf(fout_tv, "%f ", t);
      fprintf(fout_lv, "%f ", t);
      fprintf(fout_scatt, "%f ", t);
      fprintf(fout_tav, "%f ", t);
      fprintf(fout_lav, "%f ", t);
      fprintf(fout_e, "%f ", t);
      fprintf(fout_vav, "%f ", t);
      fprintf(fout_dip, "%f ", t);
      if ( mgh->safe == false )
	fprintf(fout_X, "%f ", t);
 
      double fac = 1.0/(mgh->nsample*sys.volume*(mgh->lvec-n));
      
      for ( unsigned k=0; k<mgh->nwave; k++ ){
	fprintf(fout_tv, "%f %f ", creal(mgh->c_tv[n][k]*fac),
		cimag(mgh->c_tv[n][k]*fac));
	fprintf(fout_lv, "%f %f ", creal(mgh->c_lv[n][k]*fac),
		cimag(mgh->c_lv[n][k]*fac));
	fprintf(fout_scatt, "%f %f ", creal(mgh->c_scatt[n][k]*fac),
		cimag(mgh->c_scatt[n][k]*fac));
	fprintf(fout_e, "%f %f ", creal(mgh->c_e[n][k]*fac),
		cimag(mgh->c_e[n][k]*fac));
	fprintf(fout_tav, "%f %f ", creal(mgh->c_tav[n][k]*fac),
		cimag(mgh->c_tav[n][k]*fac));
	fprintf(fout_lav, "%f %f ", creal(mgh->c_lav[n][k]*fac),
		cimag(mgh->c_lav[n][k]*fac));
	fprintf(fout_vav, "%f %f ",
		creal(mgh->c_vav[n][k]*fac), cimag(mgh->c_vav[n][k]*fac));
	fprintf(fout_dip, "%f %f ",
		creal(mgh->c_dip[n][k]*fac), cimag(mgh->c_dip[n][k]*fac));

	if ( mgh->safe == false )
	  fprintf(fout_X, "%f %f ",
		  creal(mgh->c_X[n][k]*fac), cimag(mgh->c_X[n][k]*fac));

      }

      fprintf(fout_tv, "\n");
      fprintf(fout_lv, "\n");
      fprintf(fout_scatt, "\n");
      fprintf(fout_tav, "\n");
      fprintf(fout_lav, "\n");
      fprintf(fout_vav, "\n");
      fprintf(fout_dip, "\n");
      fprintf(fout_e, "\n");

      if ( mgh->safe == false )
	fprintf(fout_X, "\n");
      
    }
	
    fclose(fout_tv); fclose(fout_lv); fclose(fout_scatt);
    fclose(fout_tav); fclose(fout_lav); fclose(fout_vav);
    fclose(fout_dip); fclose(fout_e); fclose(fout_X);
    
    mgh->i = 0;
  }

}



// Hydro profiler (atomic) 
sepprofs *sep_profs_init(char type, int lvec, int isample){
  
  sepprofs *ptr = malloc(sizeof(sepprofs));
  if ( ptr==NULL )
    sep_error("%s at line %d: Couldn't allocate memory", __func__, __LINE__);
  
  ptr->lvec = lvec;
  ptr->nsample = 0;
  ptr->isample = isample;
  ptr->type = type;

  // default - should be changed to sampler argument/option
  // dir: Spatial directio
  // dirvel: vector component evaluated in velocity
  ptr->dir = 2;
  ptr->dirvel = 0;

  ptr->momc  = sep_vector(lvec);
  ptr->svel  = sep_vector(lvec);
  ptr->dens  = sep_vector(lvec);
  ptr->temp  = sep_vector(lvec);

  return ptr;
}

void sep_profs_close(sepprofs *ptr){
 

  free(ptr->momc);  free(ptr->svel);
  free(ptr->dens);  free(ptr->temp);

}



void sep_profs_sampler(seppart *pptr, sepprofs *ptr, sepsys sys){
  double dl, *j, *rho, *sumv2, dV;
  int n, i, k=0, *numb;
  FILE *fout;

  const int dir = ptr->dir;
  const int dirvel = ptr->dirvel;

  dl = sys.length[dir]/ptr->lvec;

  if ( dir == 0 )
    dV = dl*sys.length[1]*sys.length[2];
  else if ( dir == 1 )
    dV = sys.length[0]*dl*sys.length[2];
  else 
    dV = sys.length[0]*sys.length[1]*dl;
  	     	          
  j = sep_vector(ptr->lvec);
  rho = sep_vector(ptr->lvec);
  sumv2 = sep_vector(ptr->lvec);
  numb = sep_vector_int(ptr->lvec);
  
  for ( n=0; n<sys.npart; n++ ){
    if ( pptr[n].type == ptr->type ) {

      i = pptr[n].x[dir]/dl;            
      j[i] += pptr[n].m*pptr[n].v[dirvel];
      rho[i] += pptr[n].m;
      
      for ( k=0; k<3; k++ ){			
	if ( k != dirvel ) 
	  sumv2[i] += pptr[n].m*sep_Sq(pptr[n].v[k]);			
      }

      numb[i]++;
    }
  }
  
  (ptr->nsample)++;
  double idV = 1.0/dV;

  for ( unsigned n=0; n<ptr->lvec; n++ ){
    ptr->momc[n] += j[n]*idV;
    ptr->dens[n] += rho[n]*idV;
    if ( numb[n] > 0 )
      ptr->temp[n] += sumv2[n]; 
  }
  
  free(j); free(rho); free(sumv2); free(numb);
  
  
  if ( ptr->nsample%100 == 0 ){

    fout = fopen("profs.dat", "w");
    if ( fout == NULL )
      sep_error("%s at %d: Couldn't open file", __func__, __LINE__);
    
    for ( unsigned n=0; n<ptr->lvec; n++ ){

      double insample = 1.0/ptr->nsample;
      double temp_fac = 0;
      if ( ptr->dens[n] > 0.0 ){
	ptr->svel[n] = ptr->momc[n]/ptr->dens[n];	
	temp_fac = 1.0/(2.0*dV*ptr->dens[n]);
      }     
      
      double z = (n + 0.5)*dl;
      fprintf(fout, "%f %f %f %f %f\n",  z, ptr->momc[n]*insample, 
	      ptr->dens[n]*insample, ptr->temp[n]*temp_fac, ptr->svel[n]); 
	
    }
    fclose(fout);
  }  

}				


// Hydro profiler (molecular) 
sepmprofs *sep_mprofs_init(char type, int lvec, int isample, int dir, int dirvel, 
			   int diramom){
  sepmprofs *ptr = malloc(sizeof(sepmprofs));
  if ( ptr==NULL )
    sep_error("%s at line %d: Couldn't allocate memory", __func__, __LINE__);

  ptr->lvec    = lvec;
  ptr->nsample = 0;
  ptr->isample = isample;
  ptr->type    = type;
  ptr->dir     = dir;
  ptr->dirvel  = dirvel;
  ptr->diramom = diramom;

  ptr->momc  = sep_vector(lvec);
  ptr->dens  = sep_vector(lvec);
  ptr->amom  = sep_matrix(lvec, 3);
  ptr->temp  = sep_vector(lvec);
  ptr->numb  = sep_vector(lvec);

  ptr->inertia = sep_tensor(lvec,3,3);

  return ptr;
}

void sep_mprofs_close(sepmprofs *ptr){
 
  
  free(ptr->momc);  
  sep_free_matrix(ptr->amom, ptr->lvec);
  free(ptr->dens);  
  free(ptr->temp);
  free(ptr->numb);

  sep_free_tensor(ptr->inertia, ptr->lvec, 3);
  
}



void sep_mprofs_sampler(seppart *pptr, sepmol *mols, 
			sepmprofs *ptr, sepsys sys){
  double dl, dV, w[3]={0.0};
  unsigned int n, i, k=0;
  FILE *fout;

  sep_mol_spin(pptr, mols, &sys, 0); //cm is also evaluated here
  sep_mol_velcm(pptr, mols, &sys);

  const int dir     = ptr->dir;
  const int dirvel  = ptr->dirvel;
  const int diramom = ptr->diramom;

  dl = sys.length[dir]/ptr->lvec;

  if ( dir == 0 )
    dV = dl*sys.length[1]*sys.length[2];
  else if ( dir == 1 )
    dV = sys.length[0]*dl*sys.length[2];
  else 
    dV = sys.length[0]*sys.length[1]*dl;
  
  for ( n=0; n<sys.molptr->num_mols; n++ ){
    if ( mols[n].type == ptr->type ) {
      
      i = mols[n].x[dir]/dl;            
      
      if ( i >= ptr->lvec ) 
	sep_error("%s at %d: Array index larger than vector length", 
		  __func__, __LINE__);

      ptr->momc[i] += mols[n].m*mols[n].v[dirvel];
      ptr->dens[i] += mols[n].m;
      ptr->numb[i] ++;
      for ( k=0; k<3; k++ ) ptr->amom[i][k] += mols[n].s[k];
      
      for ( k=0; k<3; k++ ){			
	if ( k != (unsigned)dirvel ) 
	  ptr->temp[i] += mols[n].m*sep_Sq(mols[n].v[k]); 		
	for ( int kk=0; kk<3; kk++ )
	  ptr->inertia[i][k][kk] += mols[n].inertia[k][kk];
      }
      
    }
  }
  
  (ptr->nsample)++;

  if ( ptr->nsample%100 == 0 ){
    
    fout = fopen("mprofs.dat", "w");
    if ( fout == NULL )
      sep_error("%s at %d: Couldn't open file", __func__, __LINE__);

    double fac = 1/(ptr->nsample*dV);    
    for ( unsigned n=0; n<ptr->lvec; n++ ){

      double **inertia = sep_matrix(3,3);

      double x = (n + 0.5)*dl;

      if ( ptr->numb[n] > 0 ){
	
	double a = ptr->dens[n]*fac;          // Mass density 
	double b = ptr->momc[n]*fac;          // Momentum current density   
	
	double av_numb = ptr->numb[n]/ptr->nsample; 
	double c = ptr->temp[n]/(2*ptr->nsample*av_numb); //Temperature

	double spin[3];
	for (  int k=0; k<3; k++ ){
	  spin[k] = ptr->amom[n][k]*fac;
	  for ( int kk=0; kk<3; kk++ )
	    inertia[k][kk] = ptr->inertia[n][k][kk]*fac; //(a*ptr->nsample);
	}

	sep_solvelineq1(w, spin, inertia, 3);
	fprintf(fout, "%e %e %e %e %e ",  x, a, b/a, c, w[diramom]);
	fprintf(fout, "%e %e %e %e %e %e %e %e %e %e %e %e\n", 
		spin[0], spin[1], spin[2], 
		inertia[0][0], inertia[0][1], inertia[0][2], 
		inertia[1][0], inertia[1][1], inertia[1][2],
		inertia[2][0], inertia[2][1], inertia[2][2]);
	
      }
      else {
	fprintf(fout, "%e %e %e %e %e %e ",  x, 0.0, 0.0, 0.0, 0.0, 0.0);
	fprintf(fout, "%e %e % e%e %e %e %e %e %e %e %e %e\n", 
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
      }
      sep_free_matrix(inertia, 3);

      
    }
    fclose(fout);
  }  
  
}				


// Molecular couple tensor ACFs
sepmcacf *sep_mcacf_init(int lvec, double tsample, double dt){
  
  sepmcacf *ptr = malloc(sizeof(sepmcacf));
  if ( ptr==NULL )
    sep_error("%s at line %d: Couldn't allocate memory", __func__, __LINE__);
  
  ptr->dt   = dt;
  ptr->lvec = lvec;
  ptr->i    = 0; 
  
  ptr->nsample = 0;

  ptr->dtsample = tsample/lvec;
  ptr->isample = (int)(ptr->dtsample/dt);
  if ( ptr->isample < 1 ) 
    sep_error("%s at line %d: isample is too small - CHECK lvec argument",
	      __func__, __LINE__);


  ptr->cacf = sep_matrix(lvec, 2);
  ptr->scouple = sep_matrix(lvec, 3);
  ptr->acouple = sep_matrix(lvec, 3);

  return ptr;
}


void sep_mcacf_sample(sepmcacf *sptr, sepatom *atoms, sepmol *mols, 
		      sepret *ret, sepsys sys){

  sep_eval_mol_couple_tensor(atoms, mols, ret, &sys);

  int index = sptr->i;

  sptr->scouple[index][0] = 0.5*(ret->T_mol[0][1]+ret->T_mol[1][0]);
  sptr->scouple[index][1] = 0.5*(ret->T_mol[0][2]+ret->T_mol[2][0]);
  sptr->scouple[index][2] = 0.5*(ret->T_mol[1][2]+ret->T_mol[2][1]);

  sptr->acouple[index][0] = 0.5*(ret->T_mol[0][1]-ret->T_mol[1][0]);
  sptr->acouple[index][1] = 0.5*(ret->T_mol[0][2]-ret->T_mol[2][0]);
  sptr->acouple[index][2] = 0.5*(ret->T_mol[1][2]-ret->T_mol[2][1]);

  (sptr->i)++;

  if ( sptr->i == sptr->lvec ){
     
    for ( int k=0; k<3; k++ ){
      for ( unsigned n=0; n<sptr->lvec; n++ ){
        for ( unsigned nn=0; nn<sptr->lvec-n; nn++ ){
          sptr->cacf[n][0] += sptr->scouple[nn][k]*sptr->scouple[n+nn][k];
	  sptr->cacf[n][1] += sptr->acouple[nn][k]*sptr->acouple[n+nn][k];
        }
      }
    }

    (sptr->nsample)++;

    FILE *fout1 = fopen("mcacf.dat", "w");
   
    for ( unsigned n=0; n<sptr->lvec; n++ ){
      double t   = (n+0.5)*sptr->dtsample;
      double fac = sys.volume/(3*(sptr->lvec-n)*sptr->nsample);
      fprintf(fout1, "%f %f %f\n", t, sptr->cacf[n][0]*fac, sptr->cacf[n][1]*fac);
    }

    fclose(fout1); 
    
    sptr->i = 0;
    
  }

}

void sep_mcacf_close(sepmcacf *ptr){
  
 
  sep_free_matrix(ptr->cacf, ptr->lvec);
  sep_free_matrix(ptr->scouple, ptr->lvec);
  sep_free_matrix(ptr->acouple, ptr->lvec);
  
  free(ptr);

}

// Molecular velocity ACF
sepmvacf *sep_mvacf_init(int lvec, double tsample, sepsys sys){
  
  sepmvacf *ptr = malloc(sizeof(sepmvacf));
  if ( ptr==NULL )
    sep_error("%s at line %d: Couldn't allocate memory", __func__, __LINE__);
  
  ptr->dt   = sys.dt;
  ptr->lvec = lvec;
  ptr->i    = 0; 

  ptr->nsample = 0;
  
  ptr->dtsample = tsample/lvec;
  ptr->isample = (int)(ptr->dtsample/sys.dt);
  if ( ptr->isample < 1 ) 
    sep_error("%s at line %d: isample is too small - CHECK lvec argument",
	      __func__, __LINE__);

  ptr->vacf = sep_vector(lvec);
  ptr->vels = sep_matrix(lvec, sys.molptr->num_mols);
  
  return ptr;
}

 
void sep_mvacf_sample(seppart *ptr, sepmvacf *vptr, sepmol *mols, sepsys sys){
    
  sep_mol_velcm(ptr, mols, &sys);

  int index = vptr->i;

  const int num_mols = sys.molptr->num_mols;

  for ( int i=0; i<num_mols; i++ )
    vptr->vels[index][i] = mols[i].v[0];
  
  (vptr->i)++;
  
  if ( vptr->i == vptr->lvec ){
    
    for ( int i=0; i<num_mols; i++ )
      for ( unsigned k=0; k<vptr->lvec; k++ )
        for ( unsigned kk=0; kk<vptr->lvec-k; kk++ )
          vptr->vacf[k] += vptr->vels[kk][i]*vptr->vels[k+kk][i];   
  
    (vptr->nsample)++;

    FILE *fout = fopen("mvacf.dat", "w");
    for ( unsigned k=0; k<vptr->lvec; k++ ){
      double t   = k*vptr->dtsample;
      double fac = 1.0/((vptr->lvec-k)*num_mols*vptr->nsample);
      fprintf(fout, "%f %f\n", t, vptr->vacf[k]*fac);
    }
    fclose(fout);

    vptr->i=0;
    
  }

}

void sep_mvacf_close(sepmvacf *ptr){

  free(ptr->vacf);
  sep_free_matrix(ptr->vels, ptr->lvec);
 
  free(ptr);

}

// Molecular angular velocity ACF
sepmavacf *sep_mavacf_init(int lvec, double tsample, sepsys sys){
  
  sepmavacf *ptr = malloc(sizeof(sepmavacf));
  if ( ptr==NULL )
    sep_error("%s at line %d: Couldn't allocate memory", __func__, __LINE__);
  
  ptr->dt   = sys.dt;
  ptr->lvec = lvec;
  ptr->i    = 0; 
  
  ptr->nsample = 0;
  
  ptr->dtsample = tsample/lvec;
  ptr->isample = (int)(ptr->dtsample/sys.dt);
  if ( ptr->isample < 1 ) 
    sep_error("%s at line %d: isample is too small - CHECK lvec argument",
	      __func__, __LINE__);

  ptr->avacf = sep_vector(lvec);
  ptr->avels = sep_matrix(lvec, sys.molptr->num_mols);
  
  return ptr;
}

 
void sep_mavacf_sample(seppart *ptr, sepmavacf *vptr, sepmol *mols, sepsys sys){
    
  sep_mol_spin(ptr, mols, &sys, 1);


  // To be removed - from here
  sep_mol_ete(ptr, mols, 'A', 0, mols[0].nuau-1, sys);
  double q = fabs(ptr[0].z);
  // to here

  int index = vptr->i;

  // To be removed - from here
  const int num_mols = sys.molptr->num_mols;

  double p[3];
  for ( int i=0; i<num_mols; i++ ){
    double mass = mols[i].m;
    for ( int k=0; k<3; k++ )
      p[k] = mols[i].ete[k]*q*mass;

    //vptr->avels[index][i] = mols[i].w[1]*p[2]-mols[i].w[2]*p[1];
    vptr->avels[index][i] = p[1];
  }
  // To here
	
  //  const int num_mols = sys.molptr->num_mols;
  //  for ( int i=0; i<num_mols; i++ )
    //    vptr->avels[index][i] = mols[i].w[0];
  
  (vptr->i)++;
  
  if ( vptr->i == vptr->lvec ){
    
    for ( int i=0; i<num_mols; i++ )
      for ( unsigned k=0; k<vptr->lvec; k++ )
        for ( unsigned kk=0; kk<vptr->lvec-k; kk++ )
          vptr->avacf[k] += vptr->avels[kk][i]*vptr->avels[k+kk][i];   
  
    (vptr->nsample)++;

    FILE *fout = fopen("mavacf.dat", "w");
    for ( unsigned k=0; k<vptr->lvec; k++ ){
      double t   = (k+0.5)*vptr->dtsample;
      double fac = 1.0/((vptr->lvec-k)*sys.molptr->num_mols*vptr->nsample);
      fprintf(fout, "%f %f\n", t, vptr->avacf[k]*fac);
    }
    fclose(fout);

    vptr->i=0;
    
  }

}

void sep_mavacf_close(sepmavacf *ptr){

  free(ptr->avacf);
  sep_free_matrix(ptr->avels, ptr->lvec);
 
  free(ptr);

}


sepmmsd *sep_mol_msd_init(int lvec, double tsample, int nk, char type, sepsys sys){

  int nmol = sys.molptr->num_mols;
  sepmmsd *ptr = malloc(sizeof(sepmmsd));

  if ( ptr==NULL )
    sep_error("%s at line %d: Couldn't allocate memory", __func__, __LINE__);
  
  ptr->msd = sep_vector(lvec);
  ptr->Fs  = sep_complex_matrix(lvec, nk);

  ptr->crossings = sep_matrix_int(nmol, 3);
  ptr->prev_pos = sep_matrix(nmol, 3);
  ptr->pos0 = sep_matrix(nmol, 3);
  ptr->k  = sep_vector(nk);

  FILE *fout = fopen("mmsd-wavevector.dat", "w");
  if ( fout==NULL )
    sep_error("%s at line %d: Couldn't open file", __func__, __LINE__);
  
  for ( int n=1; n<=nk; n++ ){
    ptr->k[n-1] = 2*SEP_PI/sys.length[0]*n;
    fprintf(fout, "%f\n", ptr->k[n-1]);
  }
  fclose(fout);
  
  ptr->dt   = sys.dt;
  ptr->lvec = lvec;
  ptr->i    = 0; 
  ptr->nk   = nk;
  ptr->nmol = nmol;
  ptr->nsample = 0;
  ptr->type = type;
  
  ptr->dtsample = tsample/lvec;
  ptr->isample = (int)(ptr->dtsample/sys.dt);
 
  return ptr;
}

void sep_mol_msd_sample(seppart *atom, sepmol *mols, sepmmsd *sptr, sepsys sys){

  sep_mol_cm(atom, mols, &sys);
	 
  int index = sptr->i;
  if ( index == 0 ){
    for ( int n=0; n<sptr->nmol; n++ ){
      for ( int k=0; k<3; k++ ) {
	sptr->prev_pos[n][k] = sptr->pos0[n][k] = mols[n].x[k];
	sptr->crossings[n][k] = 0;
      }
    }
  }

  double sd = 0, dr[3]={0.0};
  for ( int n=0; n<sptr->nmol; n++ ){
    if ( mols[n].type == sptr->type ){
      for ( int k=0; k<3; k++ ){
	if ( sptr->prev_pos[n][k] - mols[n].x[k]  > 0.5*sys.length[k] )
	  sptr->crossings[n][k] ++;
	else if ( sptr->prev_pos[n][k] - mols[n].x[k]  < -0.5*sys.length[k] )
	  sptr->crossings[n][k] --;
	
	dr[k] = mols[n].x[k] + sptr->crossings[n][k]*sys.length[k] 
	  - sptr->pos0[n][k];
	
	sptr->prev_pos[n][k] = mols[n].x[k];
      }
      sd += sep_dot(dr, dr, 3);
    }
  }

  sptr->msd[index] += sd;

  
  for ( int i=0; i<sptr->nk; i++ ){
    for ( int n=0; n<sptr->nmol; n++ ){
      if ( mols[n].type == sptr->type ){
	dr[0] = mols[n].x[0] + sptr->crossings[n][0]*sys.length[0] 
	  - sptr->pos0[n][0];
	sptr->Fs[index][i] += cexp(I*sptr->k[i]*dr[0]);     
      }
    }
  }


  index ++;
  if ( index == sptr->lvec ){
    sptr->nsample++;

    int ntype = 0;
    for ( int n=0; n<sptr->nmol; n++ )
      if ( mols[n].type == sptr->type) ntype ++;

    FILE *fout = fopen("mmsd.dat", "w");
    if ( fout==NULL )
      sep_error("%s at line %d: Couldn't open file", __func__, __LINE__);
  
    for ( int n=0; n<sptr->lvec; n++ ){      
      fprintf(fout, "%f %f \n", n*sptr->dtsample, 
	      sptr->msd[n]/(ntype*sptr->nsample));
    }
    fclose(fout);

    fout = fopen("mmsd-incoherent.dat", "w");
    if ( fout==NULL )
      sep_error("%s at line %d: Couldn't open file", __func__, __LINE__);

    for ( int n=0; n<sptr->lvec; n++ ){      
      fprintf(fout, "%f ", n*sptr->dtsample);
      for ( int i=0; i<sptr->nk; i++ )
	fprintf(fout, "%f ", creal(sptr->Fs[n][i])/(sptr->nsample*ntype));
      fprintf(fout, "\n");
    }
    fclose(fout);
    index = 0;
  }

  sptr->i = index;
}


void sep_mol_msd_close(sepmmsd *ptr){

  free(ptr->msd);
  free(ptr->k);
  
  sep_free_complex_matrix(ptr->Fs, ptr->lvec);
  sep_free_matrix(ptr->prev_pos, ptr->nmol);
  sep_free_matrix(ptr->pos0, ptr->nmol);
  sep_free_matrix_int(ptr->crossings, ptr->nmol);
  
  free(ptr);

}
