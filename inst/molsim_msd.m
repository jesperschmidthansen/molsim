


function [msd scatt _time] = molsim_msd(ptype='A', numbwave=10, haveCrossing=true, maxconfidx=10000)

	
	if ( !exist("molsim_readxyz") )
		error("molsim_msd dependent on molsim_evcorr and molsim_readxyz");
	endif
	
	indx = 0; _time(1)=0; 
	while ( 1 )
		
		filename = sprintf("molsim-%05d.xyz", indx);
		
		if ( !exist(filename) || indx > maxconfidx ) 
			break; 
		endif
		
		[npart, parttypes, pos, vel, mass, charge, lbox] = molsim_readxyz(filename, 'sep'); 
		
		if ( indx==0 )
			pidx = find( parttypes==ptype ); 
			nptypes = length(pidx);
			prevPos = pos0 = pos(pidx,:);
			
			if ( !haveCrossing )
				crossings = zeros(nptypes, 3);
			endif
			
			waves = 2*pi.*[1:1:numbwave]./lbox(2);
			
			printf("Found %d particles of type %c\n", nptypes, ptype);
			fflush(stdout);
		endif
		
		dr = pos(pidx,:) - prevPos;
		
		if ( !haveCrossing )  
			for k=1:3
				indxWrap = find( dr(:,k) > 0.5*lbox(k) );
				crossings(indxWrap,k) += 1;
				
				indxWrap = find( dr(:,k) < -0.5*lbox(k) );
				crossings(indxWrap,k) += -1;
			endfor
		else
			datacross = load(sprintf("crossings-%05d.dat", indx));
			_time(indx+1) = datacross(1,1);
			crossings = datacross(pidx+1,:);
		end

		prevPos = pos(pidx,:);
						
		dr = pos(pidx,:) + crossings.*lbox - pos0;
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

