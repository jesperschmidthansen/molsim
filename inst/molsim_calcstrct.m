##
## Usage: s = molsim_calcstrct(g, r, k, rho0)
##  
## Calculates the structure factor from the radial distribution function
##  
## Arguments: 
##   - g Radial distribution function
##   - r Radial coordinate
##   - k Wave vector
##   - rho0 System density (default 1)
##
## Output
##	 - The structure factor
##
	
function s = molsim_calcstrct(g, r, k, rho0=1)
 
	if nargin < 3 || nargin > 4 
		error("Input not correct");
	end

	h = g-1;
	s = zeros(1, length(k));
	for n=1:length(k)
		integrand =  r.^2.*h.*sin(r.*k(n))./(r.*k(n)) ;
		s(n) = 1 + 4*pi*rho0*trapz(r,integrand);
	end

end 

