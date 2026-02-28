#
# Usage: dr = ms_wrap(dr, lbox)
#
function dr = ms_wrap(dr, lbox)
	
	for k=1:3
		if dr(k) > 0.5*lbox(k)
			dr(k) = dr(k) - lbox(k);
		elseif dr(k) < -0.5*lbox(k)
			dr(k) = dr(k) + lbox(k);
		end
	end

end
