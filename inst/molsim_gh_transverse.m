
##
## usage: retval = molsim_gh_transverse(sample_time, sample_tvacf, wavectors, temperature, density, fitstrct)
## Input: 
## sample_time: Time span of the correlation function
## sample_tmacf: The tmacf - rows -> time, columns -> the different wavevectors
## temperature: System temperature
## density: System density
## fitstrct: Structure with elements
##		verbose: true or false
##		Cjj_fun: fitting for tmacf - values 'classic' or 'maxwell' (classical NS or Maxwell type) 
##		Cjj_param0: Init params for Cjj fitting - scalar for 'classic' 2-vector for 'maxwell'
##      etak_param0: Init params for etak fitting - 3-vector
##      omega: Frequencies in the Fourier-Laplace representation
##

function retval = molsim_gh_transverse(sample_time, sample_tmacf, wavectors, temperature=1.0, density=1.0, fitstrct)

	function _Cuu = _tvacf_0(_t, _param)
				
		_Cuu = temp0/rho0*exp(-k^2*param(1).*_t./rho0); 
	
	endfunction
	
	function _Cuu = _tvacf_1(_t, _param)
				
		_tauM = _param(1); _cT = _param(2); 
		_Cuu = temp0/rho0.*exp(-_t./_tauM).*cos(k*_cT.*_t); 
	
	endfunction
	
	
	function _eta = _lorentz(_k, _param)
		
		_eta0 = _param(1); _alpha = _param(2); _beta = _param(3);
		
		_eta = _eta0./(1+_alpha.* _k.^_beta);
	
	endfunction

	temp0 = temperature; rho0=density; omega = fitstrct.omega;
	
	## For leasqr call
	global verbose; verbose=[0 0]; 
	
	for n=1:length(wavectors)
	
		k = wavectors(n);
		Cjj = sample_tmacf(:,n);
		
		if strcmp(fitstrct.Cjj_fun, 'classic')
			[y_fit(:,n) param0_fit(n,:)] = leasqr(sample_time, Cjj, fitstrct.Cjj_param0, @_tvacf_0);
		elseif strcmp(fitstrct.Cjj_fun, 'maxwell')
			[y_fit(:,n) param0_fit(n,:)] = leasqr(sample_time, Cjj, fitstrct.Cjj_param0, @_tvacf_1);
		else
			error("Incorrect Cuu fitting function");
		endif 
		
		param0 = param0_fit(n,:);
	
		Cjj_w = fltrans(sample_time, hann(Cjj), fitstrct.omega);
		C0(n) = Cjj_w(1);
		
		Cjj_w = fltrans(sample_time, density.^2.*hann(Cjj), fitstrct.omega);
		eta_kw(:,n) = density./(k^2).*(Cjj(1).*density^2 - I.*omega.*Cjj_w)./Cjj_w;
		
		if fitstrct.verbose
			plot(sample_time, Cjj, 'o', sample_time, y_fit(:,n), '--');
			pause(1.0);
		endif
	endfor
	
	[etak, param1] = leasqr(wavectors, eta_kw(1,:)', fitstrct.etak_param0, @_lorentz);
	
	retval.Cjj_fit = y_fit; retval.Cjj_params = param0_fit;
	
	retval.etak_data = eta_kw(1,:);
	retval.etak_fit = etak; retval.etak_params = param1;
	retval.eta_kw = eta_kw;
	
endfunction 
