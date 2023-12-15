#include "sepcudasampler.h"

// General purpose functions
double** sep_cuda_matrix(size_t nrow, size_t ncol){
  double **ptr;
  size_t n, m;
 
  ptr = (double **)malloc(nrow*sizeof(double *));

  for (n=0; n<nrow; n++)
    ptr[n] = (double *)malloc(ncol*sizeof(double));

  for (n=0; n<nrow; n++)
    for (m=0; m<ncol; m++)
      ptr[n][m] = 0.0;

  return ptr;
}


void sep_cuda_free_matrix(double **ptr, size_t nrow){
  size_t n;

  for (n=0; n<nrow; n++)
    free(ptr[n]);

  free(ptr);
}

// The gen. hydrodynamic sampler - atomic
sepcugh* sep_cuda_sample_gh_init(sepcusys *sysptr, int lvec, unsigned nk, double dtsample){
	
	sepcugh *sptr = (sepcugh *)malloc(sizeof(sepcugh));
	
	sptr->dacf = (double **)sep_cuda_matrix(lvec, nk);
	sptr->tmacf = (double **)sep_cuda_matrix(lvec, nk);
	sptr->stress = (double **)sep_cuda_matrix(lvec, nk);
	
	sptr->wavevector = (double *)malloc(nk*sizeof(double));
	
	sptr->mcoskrArray = (double **)sep_cuda_matrix(lvec, nk);
	sptr->msinkrArray = (double **)sep_cuda_matrix(lvec, nk);

	sptr->vcoskrArray = (double **)sep_cuda_matrix(lvec, nk);
	sptr->vsinkrArray = (double **)sep_cuda_matrix(lvec, nk);
	
	sptr->stressa = (double **)sep_cuda_matrix(lvec, nk);
	sptr->stressb = (double **)sep_cuda_matrix(lvec, nk);
	
	sptr->nwaves = nk; sptr->lvec=lvec; sptr->dtsample = dtsample;
	
	sptr->index = 0; sptr->nsample = 0;
	
	FILE *fout = fopen("gh-wavevectors.dat", "w");
	if ( fout == NULL ) sep_cuda_file_error();
	
	for ( unsigned n=0; n<nk; n++ ){
		sptr->wavevector[n] = 2*SEP_CUDA_PI*(n+1)/sysptr->lbox.y;
		fprintf(fout, "%f\n", sptr->wavevector[n]);
	}
	fclose(fout);
	
	return sptr;
}	


void sep_cuda_sample_gh_free(sepcugh *ptr){
	
	sep_cuda_free_matrix(ptr->dacf, ptr->lvec);
	sep_cuda_free_matrix(ptr->tmacf, ptr->lvec);
	sep_cuda_free_matrix(ptr->stress, ptr->lvec);

	free(ptr->wavevector);
	
	sep_cuda_free_matrix(ptr->mcoskrArray, ptr->lvec);
	sep_cuda_free_matrix(ptr->msinkrArray, ptr->lvec);
	
	sep_cuda_free_matrix(ptr->vcoskrArray, ptr->lvec);
	sep_cuda_free_matrix(ptr->vsinkrArray, ptr->lvec);
	
	sep_cuda_free_matrix(ptr->stressa, ptr->lvec);
	sep_cuda_free_matrix(ptr->stressb, ptr->lvec);
	
	free(ptr);
}


void sep_cuda_sample_gh(sepcugh *sampleptr, sepcupart *pptr, sepcusys *sptr){
	
	sep_cuda_copy(pptr, 'x', 'h');
	sep_cuda_copy(pptr, 'v', 'h');
	sep_cuda_copy(pptr, 'f', 'h');
	
	unsigned index = sampleptr->index;
	
	for ( unsigned k=0; k<sampleptr->nwaves; k++ ){
	  
		double mcoskr = 0.0;	double msinkr = 0.0;
		double vcoskr = 0.0;	double vsinkr = 0.0;

		double stressa = 0.0; double stressb = 0.0;
		
		for ( unsigned n=0; n<sptr->npart; n++ ){
			double kr = sampleptr->wavevector[k]*pptr->hx[n].y;
			double mass = pptr->hx[n].w; double fx = pptr->hf[n].x;
			double velx = pptr->hv[n].x; double vely = pptr->hv[n].y;
			
			double ckr = cos(kr); double skr = sin(kr);
			mcoskr += mass*ckr; msinkr += mass*skr;
			vcoskr += mass*velx*ckr; vsinkr += mass*velx*skr;
			
			stressa += fx/sampleptr->wavevector[k]*ckr - mass*velx*vely*skr;
			stressb += fx/sampleptr->wavevector[k]*skr + mass*velx*vely*ckr;			
		}
		
		sampleptr->mcoskrArray[index][k] = mcoskr;
		sampleptr->msinkrArray[index][k] = msinkr;
		
		sampleptr->vcoskrArray[index][k] = vcoskr;
		sampleptr->vsinkrArray[index][k] = vsinkr;
		
		sampleptr->stressa[index][k] = stressa;
		sampleptr->stressb[index][k] = stressb;
	
	}
	
	(sampleptr->index)++;
	if ( sampleptr->index == sampleptr->lvec){
	
	   for ( unsigned k=0; k<sampleptr->nwaves; k++ ){
			
			for ( unsigned n=0; n<sampleptr->lvec; n++ ){
				for ( unsigned nn=0; nn<sampleptr->lvec-n; nn++ ){
					// God I miss C99!				
					double costerm = (sampleptr->mcoskrArray[nn][k])*(sampleptr->mcoskrArray[nn+n][k]);
					double sinterm = (sampleptr->msinkrArray[nn][k])*(sampleptr->msinkrArray[nn+n][k]);
					
					sampleptr->dacf[n][k]  += costerm + sinterm;
					
					costerm = (sampleptr->vcoskrArray[nn][k])*(sampleptr->vcoskrArray[nn+n][k]);
					sinterm = (sampleptr->vsinkrArray[nn][k])*(sampleptr->vsinkrArray[nn+n][k]);
					
					sampleptr->tmacf[n][k] += costerm + sinterm;
					
					double asqr = (sampleptr->stressa[nn][k])*(sampleptr->stressa[nn+n][k]);
					double bsqr = (sampleptr->stressb[nn][k])*(sampleptr->stressb[nn+n][k]);
					
					sampleptr->stress[n][k] += asqr + bsqr;
				}
			}
		}
		(sampleptr->nsample)++;

		FILE *fout_tmacf = fopen("gh-tmacf.dat", "w");
		FILE *fout_dacf = fopen("gh-dacf.dat", "w");
		FILE *fout_stress = fopen("gh-stress.dat", "w");

		if ( fout_dacf == NULL || fout_tmacf == NULL || fout_stress == NULL ){
			fprintf(stderr, "Couldn't open file(s)\n");
		}

		double volume = sptr->lbox.x*sptr->lbox.y*sptr->lbox.z;

		for ( unsigned n=0; n<sampleptr->lvec; n++ ){
			double fac = 1.0/(sampleptr->nsample*volume*(sampleptr->lvec-n));
			double t   = n*sampleptr->dtsample;
      
			fprintf(fout_dacf, "%f ", t); fprintf(fout_tmacf, "%f ", t); fprintf(fout_stress, "%f ", t); 
		 
			 for ( unsigned k=0; k<sampleptr->nwaves; k++ ) {
				 fprintf(fout_dacf, "%f ", sampleptr->dacf[n][k]*fac);
				 fprintf(fout_tmacf, "%f ", sampleptr->tmacf[n][k]*fac);
				 fprintf(fout_stress, "%f ", sampleptr->stress[n][k]*fac);
			 }
			 fprintf(fout_dacf, "\n"); fprintf(fout_tmacf, "\n");fprintf(fout_stress, "\n");
		}
	
		fclose(fout_dacf); 
		fclose(fout_tmacf);	
		fclose(fout_stress);
	
		sampleptr->index = 0;
	}
}


// The gen. hydrodynamic sampler - molecular UNDER CONSTRUCTION 
/*sepcumgh* sep_cuda_sample_mgh_init(sepcusys *sysptr, int lvec[2], unsigned nk, double dtsample){
	
	sepcumgh *sptr = (sepcumgh *)malloc(sizeof(sepcumgh));
	
	sptr->wavevector = (double *)malloc(nk*sizeof(double));

	sptr->stress = (double **)sep_cuda_matrix(lvec[0], nk);
	sptr->stressax = (double **)sep_cuda_matrix(lvec[0], nk);
	sptr->stressbx = (double **)sep_cuda_matrix(lvec[0], nk);
	sptr->stressay = (double **)sep_cuda_matrix(lvec[0], nk);
	sptr->stressby = (double **)sep_cuda_matrix(lvec[0], nk);
	sptr->stressaz = (double **)sep_cuda_matrix(lvec[0], nk);
	sptr->stressbz = (double **)sep_cuda_matrix(lvec[0], nk);

	sptr->stresslvec = lvec[0]; 
	sptr->stressindex = 0; sptr->stressnsample = 0; 

	sptr->dipole = (double **)sep_cuda_matrix(lvec[1], nk);
	sptr->dipoleax = (double **)sep_cuda_matrix(lvec[1], nk);
	sptr->dipolebx = (double **)sep_cuda_matrix(lvec[1], nk);
	sptr->dipoleay = (double **)sep_cuda_matrix(lvec[1], nk);
	sptr->dipoleby = (double **)sep_cuda_matrix(lvec[1], nk);
	sptr->dipoleaz = (double **)sep_cuda_matrix(lvec[1], nk);
	sptr->dipolebz = (double **)sep_cuda_matrix(lvec[1], nk);
	
	sptr->dipolelvec = lvec[1];
	sptr->dipoleindex = 0; sptr->dipolensample = 0;

	sptr->nwaves = nk; sptr->dtsample = dtsample;
	
	FILE *fout = fopen("mgh-wavevectors.dat", "w");
	if ( fout == NULL ) sep_cuda_file_error();
	
	for ( unsigned n=0; n<nk; n++ ){
		sptr->wavevector[n] = 2*SEP_CUDA_PI*(n+1)/sysptr->lbox.y;
		fprintf(fout, "%f\n", sptr->wavevector[n]);
	}
	fclose(fout);
	
	return sptr;
}	


void sep_cuda_sample_mgh_free(sepcumgh *ptr){
	
	free(ptr->wavevector);
	
	sep_cuda_free_matrix(ptr->stress, ptr->stresslvec);
	sep_cuda_free_matrix(ptr->stressax, ptr->stresslvec);
	sep_cuda_free_matrix(ptr->stressbx, ptr->stresslvec);
	sep_cuda_free_matrix(ptr->stressay, ptr->stresslvec);
	sep_cuda_free_matrix(ptr->stressby, ptr->stresslvec);
	sep_cuda_free_matrix(ptr->stressaz, ptr->stresslvec);
	sep_cuda_free_matrix(ptr->stressbz, ptr->stresslvec);
	
	sep_cuda_free_matrix(ptr->dipole, ptr->dipolelvec);
	sep_cuda_free_matrix(ptr->dipoleax, ptr->dipolelvec);
	sep_cuda_free_matrix(ptr->dipolebx, ptr->dipolelvec);
	sep_cuda_free_matrix(ptr->dipoleay, ptr->dipolelvec);
	sep_cuda_free_matrix(ptr->dipoleby, ptr->dipolelvec);
	sep_cuda_free_matrix(ptr->dipoleaz, ptr->dipolelvec);
	sep_cuda_free_matrix(ptr->dipolebz, ptr->dipolelvec);
	
	free(ptr);
}


void sep_cuda_print_current_corr(sepcusys *sptr, sepcumgh *sampler, const char quantity, const char *filename){
	double **corr, **a, **b, **ax, **bx, **ay, **by, **az, **bz;
   	unsigned lvec; int nsample; 	
	
	if ( quantity=='S' ){
		corr = sampler->stress; a = sampler->stressa; b = sampler->stressb;
		lvec = sampler->stresslvec; nsample = sampler->stressnsample;
	}
	else if ( quantity=='D' ){
		corr = sampler->dipole; 
		ax = sampler->dipoleax; bx = sampler->dipolebx;
		ay = sampler->dipoleay; by = sampler->dipoleby;
		az = sampler->dipoleaz; bz = sampler->dipolebz;

		lvec = sampler->dipolelvec; nsample = sampler->dipolensample;
	}
	else {
		fprintf(stderr, "Invalid argument\n");
		return;
	}

	for ( unsigned k=0; k<sampler->nwaves; k++ ){
		
		for ( unsigned n=0; n<lvec; n++ ){
			for ( unsigned nn=0; nn<lvec-n; nn++ ){
				if ( quantity=='S' ){
					double asqr = a[nn][k]*a[nn+n][k];
					double bsqr = b[nn][k]*b[nn+n][k];
					
					corr[n][k] += asqr + bsqr;
				}
				else if ( quantity=='D' ){
					double asqr = ax[nn][k]*ax[nn+n][k] + ay[nn][k]*ay[nn+n][k] + az[nn][k]*az[nn+n][k];
					double bsqr = bx[nn][k]*bx[nn+n][k] + by[nn][k]*by[nn+n][k] + bz[nn][k]*bz[nn+n][k];
					
					corr[n][k] += (asqr + bsqr)/3.0;
				}
			}	
		}
	}

	FILE *fout = fopen(filename, "w");
	if ( fout == NULL ) fprintf(stderr, "Couldn't open file(s)\n");
		
	double volume = sptr->lbox.x*sptr->lbox.y*sptr->lbox.z;
		
	for ( unsigned n=0; n<lvec; n++ ){
		double fac = 1.0/(nsample*volume*(lvec-n));
		double t   = n*sampler->dtsample;
	
		fprintf(fout, "%f ", t); 
 
		for ( unsigned k=0; k<sampler->nwaves; k++ ) {
			fprintf(fout, "%f ", corr[n][k]*fac); 
		}
		fprintf(fout, "\n") ;
	}
	
	fclose(fout); 
	
}

void sep_cuda_sample_mgh(sepcumgh *sampleptr, sepcupart *pptr, sepcusys *sptr, sepcumol *mptr){

	if ( !pptr->sptr->molprop ) {
		fprintf(stderr, "Mol. properties flag not set to 'on' - stress correlator not calculated\n");
		return;
	}

	if ( !sptr->cmflag ) 	
		sep_cuda_mol_calc_cmprop(pptr, mptr); // Calculations done and saved on host
	
	// Forces on molecular  - wavevector depedent stress and mech. properties
	sep_cuda_mol_calc_forceonmol(pptr, mptr);
	sep_cuda_copy(pptr, 'M', 'h');

	// Dipole moments - dielectric properties
	sep_cuda_mol_calc_dipoles(pptr, mptr);

	// Just to ease the symbolism
	unsigned idxS = sampleptr->stressindex;
	unsigned idxD = sampleptr->dipoleindex;
	for ( unsigned k=0; k<sampleptr->nwaves; k++ ){
	  
		double stressa[3] = {0.0}; double stressb[3] = {0.0};
		double dipolea[3] = {0.0}; double dipoleb[3] = {0.0};

		for ( unsigned m=0; m<mptr->nmols; m++ ){
			
			double kr3[3];
			kr3[0] = sampleptr->wavevector[k]*mptr->hx[m].x;
			kr3[1] = sampleptr->wavevector[k]*mptr->hx[m].y;
			kr3[2] = sampleptr->wavevector[k]*mptr->hx[m].z;
		
			double ckr3[3]; double skr3[3]; 
		   for ( int i=0; i<3; i+ ){ ckr3[i] = cos(kr3[i]); skr3[i] = sin(kr3[i]) };

			double mass = mptr->masses[m]; double fx = mptr->hf[m].x;	
			double velx = mptr->hv[m].x; double vely = mptr->hv[m].y;
			
			stressa += fx/sampleptr->wavevector[k]*ckr - mass*velx*vely*skr;
			stressb += fx/sampleptr->wavevector[k]*skr + mass*velx*vely*ckr;		

			double kr3[3];
			if ( k>0 ) {
				kr3[0] = sampleptr->wavevector[k-1]*mptr->hx[m].x;
				kr3[1] = sampleptr->wavevector[k-1]*mptr->hx[m].y;
				kr3[2] = sampleptr->wavevector[k-1]*mptr->hx[m].z;
			}

	
				dipolea[0] += mptr->hpel[m].x*cos(kr3[0]); 
			dipoleb[0] += mptr->hpel[m].x*sin(kr3[0]); 
			dipolea[1] += mptr->hpel[m].y*cos(kr3[1]); 
			dipoleb[1] += mptr->hpel[m].y*sin(kr3[1]); 
			dipolea[2] += mptr->hpel[m].z*cos(kr3[2]); 
			dipoleb[2] += mptr->hpel[m].z*sin(kr3[2]); 
	
		}
				
		sampleptr->stressa[idxS][k] = stressa;	sampleptr->stressb[idxS][k] = stressb;
		
		sampleptr->dipoleax[idxD][k] = dipolea[0];	sampleptr->dipolebx[idxD][k] = dipoleb[0];
		sampleptr->dipoleay[idxD][k] = dipolea[1];	sampleptr->dipoleby[idxD][k] = dipoleb[1];
		sampleptr->dipoleaz[idxD][k] = dipolea[2];	sampleptr->dipolebz[idxD][k] = dipoleb[2];
	}
	
	idxS++; idxD++;

	if ( idxS == sampleptr->stresslvec ){
		(sampleptr->stressnsample)++;
		idxS = 0;
		sep_cuda_print_current_corr(sptr, sampleptr, 'S', "mgh-stress.dat");
	}	
	
	if ( idxD == sampleptr->dipolelvec ){
		(sampleptr->dipolensample)++;
		idxD = 0;
		sep_cuda_print_current_corr(sptr, sampleptr, 'D', "mgh-dipole.dat");
	}	

	sampleptr->stressindex=idxS; 
	sampleptr->dipoleindex=idxD;

}
*/
// Wavevector dipole/polarization correlator  
sepcusampler_dipole* sep_cuda_sample_dipole_init(sepcusys *sptr, int lvec, unsigned nk, double dtsample){
	
	sepcusampler_dipole *sampleptr = (sepcusampler_dipole *)malloc(sizeof(sepcusampler_dipole));
	if ( sampleptr==NULL ) sep_cuda_mem_error();

	sampleptr->lvec = lvec; 
	sampleptr->index = 0; 
	sampleptr->nsample = 0; 

	sampleptr->nwaves = nk; 
	sampleptr->dtsample = dtsample;
	
	sampleptr->corr = (double **)sep_cuda_matrix(lvec, nk);
	sampleptr->dipoleax = (double **)sep_cuda_matrix(lvec, nk);
	sampleptr->dipolebx = (double **)sep_cuda_matrix(lvec, nk);
	sampleptr->dipoleay = (double **)sep_cuda_matrix(lvec, nk);
	sampleptr->dipoleby = (double **)sep_cuda_matrix(lvec, nk);
	sampleptr->dipoleaz = (double **)sep_cuda_matrix(lvec, nk);
	sampleptr->dipolebz = (double **)sep_cuda_matrix(lvec, nk);
	
	if ( sampleptr->corr == NULL || 
			sampleptr->dipoleax == NULL || sampleptr->dipolebx == NULL || 
			sampleptr->dipoleay == NULL || sampleptr->dipoleby == NULL ||
			sampleptr->dipoleaz == NULL || sampleptr->dipolebz == NULL )
		sep_cuda_mem_error();


	sampleptr->wavevector = (double *)malloc(nk*sizeof(double));
	if ( sampleptr->wavevector == NULL ) sep_cuda_mem_error();

	FILE *fout = fopen("wavevectors-polcor.dat", "w");
	if ( fout == NULL ) sep_cuda_file_error();
	
	// We include the zero vector
	for ( unsigned n=0; n<nk; n++ ){
		sampleptr->wavevector[n] = 2*SEP_CUDA_PI*n/sptr->lbox.y;
		fprintf(fout, "%f\n", sampleptr->wavevector[n]);
	}
	fclose(fout);
	
	return sampleptr;
}	


void sep_cuda_sample_dipole_free(sepcusampler_dipole *ptr){
	
	free(ptr->wavevector);

	sep_cuda_free_matrix(ptr->corr, ptr->lvec);
	sep_cuda_free_matrix(ptr->dipoleax, ptr->lvec);
	sep_cuda_free_matrix(ptr->dipolebx, ptr->lvec);
	sep_cuda_free_matrix(ptr->dipoleay, ptr->lvec);
	sep_cuda_free_matrix(ptr->dipoleby, ptr->lvec);
	sep_cuda_free_matrix(ptr->dipoleaz, ptr->lvec);
	sep_cuda_free_matrix(ptr->dipolebz, ptr->lvec);
	
	free(ptr);
}



void sep_cuda_sample_dipole(sepcusampler_dipole *sampleptr, sepcupart *pptr, sepcusys *sptr, sepcumol *mptr){

	if ( !sptr->cmflag ) 	
		sep_cuda_mol_calc_cmprop(pptr, mptr); // Calculations done and saved on host
	
	// Dipole moments - dielectric properties
	sep_cuda_mol_calc_dipoles(pptr, mptr);

	for ( unsigned k=0; k<sampleptr->nwaves; k++ ){
	  
		double dipolea[3] = {0.0}; double dipoleb[3] = {0.0};

		for ( unsigned m=0; m<mptr->nmols; m++ ){
			
			double kr3[3];
			kr3[0] = sampleptr->wavevector[k]*mptr->hx[m].x;
			kr3[1] = sampleptr->wavevector[k]*mptr->hx[m].y;
			kr3[2] = sampleptr->wavevector[k]*mptr->hx[m].z;
	
			dipolea[0] += mptr->hpel[m].x*cos(kr3[0]); 
			dipoleb[0] += mptr->hpel[m].x*sin(kr3[0]); 

			dipolea[1] += mptr->hpel[m].y*cos(kr3[1]); 
			dipoleb[1] += mptr->hpel[m].y*sin(kr3[1]); 
			
			dipolea[2] += mptr->hpel[m].z*cos(kr3[2]); 
			dipoleb[2] += mptr->hpel[m].z*sin(kr3[2]); 
	
		}
				
		unsigned idx = sampleptr->index;	
		sampleptr->dipoleax[idx][k] = dipolea[0];	sampleptr->dipolebx[idx][k] = dipoleb[0];
		sampleptr->dipoleay[idx][k] = dipolea[1];	sampleptr->dipoleby[idx][k] = dipoleb[1];
		sampleptr->dipoleaz[idx][k] = dipolea[2];	sampleptr->dipolebz[idx][k] = dipoleb[2];
	}
	
	(sampleptr->index)++;
	
	if ( sampleptr->index == sampleptr->lvec ){
		(sampleptr->nsample)++;
		sampleptr->index = 0;

		for ( unsigned k=0; k<sampleptr->nwaves; k++ ){
			
			for ( unsigned n=0; n<sampleptr->lvec; n++ ){
				for ( unsigned nn=0; nn<sampleptr->lvec-n; nn++ ){
					double asqr = sampleptr->dipoleax[nn][k]*sampleptr->dipoleax[nn+n][k] + 
						sampleptr->dipoleay[nn][k]*sampleptr->dipoleay[nn+n][k] + 
						sampleptr->dipoleaz[nn][k]*sampleptr->dipoleaz[nn+n][k];
					double bsqr =  sampleptr->dipolebx[nn][k]*sampleptr->dipolebx[nn+n][k] + 
						sampleptr->dipoleby[nn][k]*sampleptr->dipoleby[nn+n][k] + 
						sampleptr->dipolebz[nn][k]*sampleptr->dipolebz[nn+n][k];

					sampleptr->corr[n][k] += (asqr + bsqr)/3.0;
				}
			}	
		}

		FILE *fout = fopen("polcor.dat", "w");
		if ( fout == NULL ) sep_cuda_file_error();
		
		double volume = sptr->lbox.x*sptr->lbox.y*sptr->lbox.z;
		
		for ( unsigned n=0; n<sampleptr->lvec; n++ ){
			double fac = 1.0/(sampleptr->nsample*volume*(sampleptr->lvec-n));
			double t   = n*sampleptr->dtsample;
	
			fprintf(fout, "%f ", t); 
 
			for ( unsigned k=0; k<sampleptr->nwaves; k++ ) {	
				fprintf(fout, "%f ", sampleptr->corr[n][k]*fac); 
			}
			fprintf(fout, "\n") ;
		}
	
		fclose(fout);
	}	

}

sepcusampler_stress* sep_cuda_sample_stress_init(sepcusys *sysptr, int lvec, unsigned nk, double dtsample){
	
	sepcusampler_stress *sptr = (sepcusampler_stress *)malloc(sizeof(sepcusampler_stress));
	if ( sptr==NULL ) sep_cuda_mem_error();

	sptr->lvec = lvec; 
	sptr->index = 0; 
	sptr->nsample = 0; 

	sptr->nwaves = nk; 
	sptr->dtsample = dtsample;

	sptr->corr = (double **)sep_cuda_matrix(lvec, nk);
	sptr->stressa = (double **)sep_cuda_matrix(lvec, nk);
	sptr->stressb = (double **)sep_cuda_matrix(lvec, nk);

	if ( sptr->corr == NULL || sptr->stressa == NULL || sptr->stressb == NULL )
		sep_cuda_mem_error();

	sptr->wavevector = (double *)malloc(nk*sizeof(double));
	if ( sptr->wavevector == NULL ) sep_cuda_mem_error();

	FILE *fout = fopen("wavevectors-stresscor.dat", "w");
	if ( fout == NULL ) sep_cuda_file_error();
	
	for ( unsigned n=0; n<nk; n++ ){
		sptr->wavevector[n] = 2*SEP_CUDA_PI*(n+1)/sysptr->lbox.y;
		fprintf(fout, "%f\n", sptr->wavevector[n]);
	}
	fclose(fout);
	
	return sptr;
}	


void sep_cuda_sample_stress_free(sepcusampler_stress *ptr){
	
	free(ptr->wavevector);
	
	sep_cuda_free_matrix(ptr->corr, ptr->lvec);
	sep_cuda_free_matrix(ptr->stressa, ptr->lvec);
	sep_cuda_free_matrix(ptr->stressb, ptr->lvec);

	free(ptr);
}

void sep_cuda_sample_stress(sepcusampler_stress *sampleptr, sepcupart *pptr, sepcusys *sptr, sepcumol *mptr){

	if ( !pptr->sptr->molprop ) {
		fprintf(stderr, "Mol. properties flag not set to 'on' - stress correlator not calculated\n");
		return;
	}

	if ( !sptr->cmflag ) 	
		sep_cuda_mol_calc_cmprop(pptr, mptr); // Calculations done and saved on host
	
	sep_cuda_mol_calc_forceonmol(pptr, mptr);
	sep_cuda_copy(pptr, 'M', 'h');

	unsigned idx = sampleptr->index;
	

	for ( unsigned k=0; k<sampleptr->nwaves; k++ ){

		double stressa = 0.0; double stressb = 0.0;
                
		for ( unsigned m=0; m<mptr->nmols; m++ ){
			double kr = sampleptr->wavevector[k]*mptr->hx[m].y;

			double mass = mptr->masses[m]; double fx = mptr->hf[m].x;
			double velx = mptr->hv[m].x; double vely = mptr->hv[m].y;
							
			double ckr = cos(kr); double skr = sin(kr);
						
			stressa += fx/sampleptr->wavevector[k]*ckr - mass*velx*vely*skr;
			stressb += fx/sampleptr->wavevector[k]*skr + mass*velx*vely*ckr;                       
		}

		sampleptr->stressa[idx][k] = stressa;
	   	sampleptr->stressb[idx][k] = stressb;
	}
	
	(sampleptr->index)++;
	
	if ( sampleptr->index == sampleptr->lvec ){
		(sampleptr->nsample)++;
		sampleptr->index = 0;
	
		for ( unsigned k=0; k<sampleptr->nwaves; k++ ){
			
			for ( unsigned n=0; n<sampleptr->lvec; n++ ){
				for ( unsigned nn=0; nn<sampleptr->lvec-n; nn++ ){
					double asqr = (sampleptr->stressa[nn][k])*(sampleptr->stressa[nn+n][k]);
					double bsqr = (sampleptr->stressb[nn][k])*(sampleptr->stressb[nn+n][k]);

					sampleptr->corr[n][k] += asqr + bsqr;
				}
			}	
		}

		FILE *fout = fopen("stresscor.dat", "w");
		if ( fout == NULL ) sep_cuda_file_error();
		
		double volume = sptr->lbox.x*sptr->lbox.y*sptr->lbox.z;
		
		for ( unsigned n=0; n<sampleptr->lvec; n++ ){
			double fac = 1.0/(sampleptr->nsample*volume*(sampleptr->lvec-n));
			double t   = n*sampleptr->dtsample;
	
			fprintf(fout, "%f ", t); 
 
			for ( unsigned k=0; k<sampleptr->nwaves; k++ ) {	
				fprintf(fout, "%f ", sampleptr->corr[n][k]*fac); 
			}
			fprintf(fout, "\n") ;
		}
	
		fclose(fout);
	}	

}


