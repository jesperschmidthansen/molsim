##
## usage:  [msd scatt time] = molsim_msd(ptype, numbwave, haveCrossing (bool), maximum config file index)
##
## Calculates the mean square distance from configuration files molsim_%05d.xyz and crossing files
## crossings-%05d.dat. If the later does not exist an attempt to estimate the particles' box crossings.
## 


function [msd scatt _time] = molsim_msd(ptype='A', numbwave=10, haveCrossing=true, maxconfidx=10000)

	
	if ( !exist("molsim_readxyz") )
		error("molsim_msd dependent on molsim_evcorr and molsim_readxyz");
	endif
	
	if ( nargin > 4 )
		error("Usage: molsim_msd(part. type, no. wavevectors, have crossing files, max. configurations)");
	endif
	
	_time=[];
	
	indx = 0;  
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
		else 
			if ( !haveCrossing )  
				## Attempts to find the true distance from prev. positions
				dr = pos(pidx,:) - prevPos;
				for k=1:3
					indxWrap = find( dr(:,k) > 0.5*lbox(k) );
					crossings(indxWrap,k) += 1;
					
					indxWrap = find( dr(:,k) < -0.5*lbox(k) );
					crossings(indxWrap,k) += -1;
				endfor
				
				prevPos = pos(pidx,:);
				
			else
				filenameCross = sprintf("crossings-%05d.dat", indx);
				if  ( !exist(filenameCross) )
					error("Cannot find crossing files\n");
				endif
				
				datacross = load(filenameCross);
				_time(indx) = datacross(1,1);
				crossings = datacross(pidx+1,:);
			endif
	
			dr = pos(pidx,:) + crossings.*lbox - pos0;
			sd = sum(dot(dr', dr'))./nptypes;
					
			## MSD
			msd(indx) = sd;
			
			## SCATT
			for nk=1:numbwave
				sumsc = 0.0;
				for np=1:nptypes
					sumsc = sumsc + real(exp(-I.*waves(nk)*dr(np,1)));
				end
				scatt(indx, nk) = sumsc/nptypes;
			end
			
		endif
		
		indx++;
		
		printf("\rProcessed %d  ", indx); fflush(stdout);
	endwhile 
	
	printf("\n");
endfunction 

