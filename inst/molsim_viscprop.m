
function retval = molsim_viscprop(_time, stress, nblocks=10, volume=1, temp=1, ... 
									_verbose=true, omega=logspace(-2,1))

	if ( !exist("molsim_calccf") || !exist(fltrans) || !exist(hann) )
		error(" The functions molsim_calccf, fltrans, and hann are needed, but at least one not found - bailing out");
	endif
	
	if ( columns(stress) != 3 )
		error("The stress must be on the form nx3");
	endif
	
	if ( rows(_time) != rows(stress) )
		error("Number of rows in time and stress arrays must be the same");
	endif
	
	for n=1:3
		css(:,n) = molsim_calccf(stress(:,n), stress(:,n), nblocks);
	endfor
	
	## The actual stress correlation function averaged over 
	## the three off-diagonal elements	
	retval.css_av = mean(css')';

	## The sample time window
	nr_0 = rows(css); retval.time = linspace(0, _time(2)*nr_0, nr_0);
	
	## Attempt to estimate a convergence from the stress autocorrelation
	## function tail (defined as last 1/3 of data set).
	nr_1 = int64(nr_0*2.0/3.0);
	retval.mean_tail = mean(retval.css_av(nr_1:end));
	retval.std_tail = std(retval.css_av(nr_1:end));
	
	## Actual viscosity
	retval.integrale_css = cumtrapz(retval.css_av);
	retval.integrale_css_hann = cumtrapz(hann(retval.css_av));
	
	retval.eta0_estimate = retval.integrale_css_hann(end)*volume/temp;
	
	## Complex viscosity and modulus
	retval.etaw = fltrans(_time, hann(retval.css_av), omega);
	retval.modulus = I*omega.*retval.etaw;
		
	## Some verbose plots; running integrals, etc
	if _verbose
		subplot(3,1,1);
		semilogx(retval.time(2:end), retval.css_av(2:end)); 
		grid on
		
		subplot(3,1,2);
		plot(retval.time(nr_1:end), retval.css_av(nr_1:end));
		grid on
		
		subplot(3,1,3);
		plot(retval.integrale_css); hold on;
		plot(retval.integrale_css_hann); hold off
		
		printf("Convergence %.2e \pm %.2e \n", retval.mean_tail, retval.std_tail);
		fflush(stdout);
	endif
	
	
endfunction
