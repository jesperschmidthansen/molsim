##
## Usage: css = molsim_calccf(a, b, nblocks)
##  
## Calculates (directly) the correlation function between data arrays a and b
##  
## Arguments: 
##   - a and b; data arrays
##   - nblocks; Number of blocks the arrays are divided into (default 1)
##
## Output
##	 - The correlations function
##
	
function css = molsim_calccf(data1, data2, nblocks=1)

	if ( !exist("evcorr.oct") )
		error("molsim_calccf needs evcorr.oct");
	endif
	
	if ( nargin > 3 )
		error("molsim_calccf called with too many arguments");
	endif
	
	ldata = rows(data1); 
	nseries = columns(data2);
	lblock = floor(ldata/nblocks);

	cor = zeros(lblock,1);
	
	for n=1:nblocks
		for m=1:nseries
			a = data1((n-1)*lblock + 1:n*lblock, m);
			b = data2((n-1)*lblock + 1:n*lblock, m);
			cor = cor + evcorr(a,b);
		endfor
	endfor
	
	css = cor./(nseries*nblocks);
	
endfunction	
