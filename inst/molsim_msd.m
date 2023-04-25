


function [msd scatt] = molsim_msd(ptype='A', numbwave=10, maxconfidx=10000)

	
	if ( !exist("molsim_readxyz") )
		error("molsim_msd dependent on molsim_evcorr and molsim_readxyz");
	endif
	
	indx = 0; 
	while ( 1 )
		
		filename = sprintf("molsim-%05d.xyz", indx);
		
		if ( !exist(filename) || indx > maxconfidx ) 
			break; 
		endif
		
		[npart, parttypes, pos, vel, mass, charge, lbox] = molsim_readxyz(filename, 'sep'); 
		
		if ( indx==0 )
			pidx = find( parttypes==ptype ); nptypes = length(pidx);
			prevPos = pos0 = pos(pidx,:);
			crossings = zeros(nptypes, 3);
			
			waves = 2*pi.*[1:1:numbwave]./lbox(1);
			
			printf("Found %d particles of type %c\n", nptypes, ptype);
			fflush(stdout);
		endif
		
		dr = pos(pidx,:) - prevPos;
		
		for k=1:3
			indxWrap = find( dr(:,k) > 0.5*lbox(k) );
			crossings(indxWrap,k) += 1;
			
			indxWrap = find( dr(:,k) < -0.5*lbox(k) );
			crossings(indxWrap,k) += -1;
		endfor
		prevPos = pos(pidx,:);
		
		dr = pos(pidx,:) - crossings.*lbox - pos0;
				
		sd = sum(dot(dr', dr'))./nptypes;
		
		## MSD
		indx++; msd(indx) = sd;
		
		## SCATT
		for nk=1:numbwave
			sumsc = 0.0;
			for np=1:nptypes
				sumsc = sumsc + real(exp(-I.*waves(nk)*dr(np,1)));
			end
			scatt(indx, nk) = sumsc/nptypes;
		end
		
		printf("\rProcessed %d  ", indx);fflush(stdout);
	endwhile 
	
	printf("\n");
endfunction 



#{
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
#}
