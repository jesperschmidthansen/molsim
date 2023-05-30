

function retval = molsim_gh_transverse(sample_time, sample_tvacf, wavectors, temperature, density, ...
									 param0 = [1,1], param1 = [1, 1, 1], omega=linspace(0, 1,10))

	function _Cuu = _tvacf_0(_t, _param)
				
		_eta0 = _param(1); 
		_Cuu = temp0/rho0*exp(-k^2*_eta0.*_t./rho0); 
	
	endfunction
	
	function _Cuu = _tvacf_1(_t, _param)
				
		_tauM = _param(1); _cT = _param(2); 
		_Cuu = temp0/rho0.*exp(-_t./_tauM).*cos(k*_cT.*_t); 
	
	endfunction
	
	
	function _eta = _lorentz(_k, _param)
		
		_eta0 = _param(1); _alpha = _param(2); _beta = _param(3);
		
		_eta = _eta0./(_alpha + _k.^_beta);
	
	endfunction
	
	temp0 = temperature; rho0=density;
	
	global verbose; verbose(1)=verbose(2)=0;
	
	for n=1:length(wavectors)
	
		k = wavectors(n);
		Cuu = sample_tvacf(:,n);
		
		[y_fit(:,n) param0_fit(n,:)] = leasqr(sample_time, Cuu, param0, @_tvacf_1);

		param0 = param0_fit(n,:);
	
		C0(n) = real(fltrans(sample_time, hann(Cuu), omega))(1);
		
	endfor
	
	etak_data = temperature./(wavectors.^2.*C0');
	[etak, param1] = leasqr(wavectors, etak_data', param1, @_lorentz);
	
	retval.Cuu_fit = y_fit; retval.Cuu_params = param0_fit;
	
	retval.etak_data = etak_data;
	retval.etak_fit = etak; retval.etak_params = param1;
	
endfunction 
