##
## usage:[pz, pdens, pvel, ptemp] = molsim_pprofs(z, dens, vel, temp, threshold)
##
## Pretty profiles. Cut and shifts the data point
## 
## Arguments: 
##  - z; position data points
##  - dens; density data points
##  - vel; velocity data points
##  - temp; temperature data points
##  - threshold; Cut density threshold value (Default: 1e-2) 
##
## Output:
##  The corresponding data, but cut and shifted 
##

function [pz, pdens, pvel, ptemp] = molsim_pprofs(z, dens, vel, temp, threshold=1.0e-2)

	if ( nargin < 4 || nargin > 5 )
		error("Number of arguments is 4 or 5");
	end

	for n=1:length(z)
		if dens(n)>threshold
			idx0 = n;
			break;
		end
	end
	
	for n=length(z):-1:1
		if dens(n)>threshold
			idx1 = n;
			break;
		end
	end

	pz = z(idx0:idx1) - z(idx0); pdens = dens(idx0:idx1); 
	pvel = vel(idx0:idx1); ptemp = temp(idx0:idx1);

end
