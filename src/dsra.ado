*! version 1.0.0  08may2020
/*-----------------------------------------------------------------------------
	double selection for regression adjustment
-----------------------------------------------------------------------------*/
program dsra, eclass
	version 16.0	
	
	/*----- preserve START -----*/
	preserve

	/*----- parse syntax ----*/
	tempvar touse1 touse2
	_dsra_parse_syntax, touse1(`touse1') touse2(`touse2') : `0'

	/*----- compute -----*/
	Compute, depvar(`r(depvar)')		///
		x1_vars(`r(x1_vars)')		///
		x2_var(`r(x2_vars)')		///
		x3i_vars(`r(x3i_vars)')		///
		x3e_vars(`r(x3e_vars)')		///
		touse1(`r(touse1)')		///
		touse2(`r(touse2)')		///
		lassoopts(`r(lassoopts)')	///
		adjust(`r(adjust)')
	
	/*----- preserve END -----*/
	restore

	/*--- Display -----*/
	_dsra_display
end

					//----------------------------//
					//  compute
					//----------------------------//
program Compute
	syntax, depvar(string)		///
		x1_vars(string)		///
		x2_var(string)		///
		x3i_vars(string)	///
		x3e_vars(string)	///
		touse1(string)		///
		touse2(string)		///
		adjust(passthru)	///
		[lassoopts(string)]
	
	/*----- 1. In sample 2, select x3_sel ------------------------------*/
	SelectX3, x2_var(`x2_var') 	///
		x3i_vars(`x3i_vars') 	///
		x3e_vars(`x3e_vars')	///
		touse2(`touse2')	///
		lassoopts(`lassoopts')	
	local x3_sel `r(x3_sel)'
	local x3_nosel `r(x3_nosel)'

	/*----- 2. In sample 2, compute gamma_sel and its variance ---------*/
	ComputeGamma, x2_var(`x2_var')	///
		x3_sel(`x3_sel')	///
		x3_nosel(`x3_nosel')	///
		touse2(`touse2')	
	tempname vg
	mat `vg' = e(V)

	/*----- 3. In sample 1, predict x2 ---------------------------------*/	
	tempvar x2hat
	PredictX2, x2_var(`x2_var')	///
		touse1(`touse1')	///
		x2hat(`x2hat')		///
		x3_sel(`x3_sel')	///
		x3e_vars(`x3e_vars')
	local x3_sel `r(x3_sel)'
	local x3e_nosel `r(x3e_nosel)'

	/*----- 4. In sample 1, run double selection lasso for main model --*/
	ComputeBeta, depvar(`depvar')	///
		x1_vars(`x1_vars')	///
		x2hat(`x2hat')		///
		x3i_vars(`x3i_vars')	///
		x3_sel(`x3_sel')	///
		x3e_nosel(`x3e_nosel')	///
		touse1(`touse1')	///
		lassoopts(`lassoopts')
	local controls_sel `r(controls_sel)'
	local b2 = `r(b2)'
	local xa `r(xa)'
	tempname vb
	mat `vb' = e(V)

	/*------ 5. Adjust variance ----------------------------------------*/
	AdjustVar, vg(`vg') 		///
		vb(`vb') 		///
		b2(`b2') 		///
		x3_sel(`x3_sel') 	///
		xa(`xa')		///
		touse1(`touse1')	///
		touse2(`touse2')	///
		`adjust'

	/*------ 6. post result --------------------------------------------*/
	local x3e_sel `controls_sel' 
	PostResult, 				///
		depvar(`depvar')		///
		x1_vars(`x1_vars')		///
		x2_var(`x2_var')		///
		x3i_vars(`x3i_vars')		///
		x3e_vars(`x3e_vars')		///
		touse1(`touse1') 		///
		touse2(`touse2') 		///
		x3e_sel(`x3e_sel')		///
		`adjust'
end

					//----------------------------//
					//  post coeff
					//----------------------------//
program PostResult, eclass
	syntax , depvar(string)		///
		x1_vars(string)		///
		x2_var(string)		///
		x3i_vars(string)	///
		x3e_vars(string)	///
		touse1(string)		///
		touse2(string)		///
		adjust(string)		///
		[x3e_sel(string)]

	local vars `x1_vars' `x2_var' `x3i_vars'
	PostCoef, vars(`vars')

	sum `touse1' if `touse1', meanonly
	eret scalar N_1 = r(sum)

	sum `touse2' if `touse2', meanonly
	eret scalar N_2 = r(sum)

	eret local vcetype Robust
	eret local vce robust

	eret local depvar `depvar'
	eret local controls `x3e_vars'
	eret local controls_sel `x3e_sel'

	eret scalar k_controls = `:list sizeof x3e_vars'
	eret scalar k_controls_sel = `:list sizeof x3e_sel'

	eret local title "Double-selection on two samples"

	if (`"`adjust'"' == "noadjust") {
		eret local cmd dsra_unadjust
	}
	else {
		eret loca cmd dsra
	}

end

					//----------------------------//
					//  get b and V
					//----------------------------//
program PostCoef, eclass
	syntax , vars(string)	

	tempname b V
	mat `b' = e(b)
	mat `V' = e(V)

	local k : list sizeof vars
	mat `b' = `b'[1, 1..`k']
	mat `V' = `V'[1..`k', 1..`k']

	mat colname `b' = `vars'
	mat colname `V' = `vars'
	mat rowname `V' = `vars'

	eret post `b' `V', buildfv
end

					//----------------------------//
					//  select x3
					//----------------------------//
program SelectX3, rclass
	syntax, x2_var(string)		///
		x3i_vars(string)	///
		x3e_vars(string)	///
		touse2(string)		///
		[lassoopts(string) ]

	/*----- 1. In sample 2, select x3_sel -------------*/
	di as txt "Step 1: In sample 2, get x3_sel"

	// using sample 2, run lasso x2 on x3
	di as txt "Estimating lasso for {it:`x2_var'}"
	local x3 `x3i_vars' `x3e_vars'
	qui lasso linear `x2_var' `x3' if `touse2', `lassoopts'

	// get x3_sel and x3_nosel
	local allvars `e(allvars)'
	local x3_sel `e(allvars_sel)'
	local x3_nosel : list allvars - x3_sel

	di

	ret local x3_sel `x3_sel'
	ret local x3_nosel `x3_nosel'
end
					//----------------------------//
					// compute gamma sel 
					//----------------------------//
program ComputeGamma
	syntax, x2_var(string)		///
		x3_sel(string)		///
		x3_nosel(string)	///
		touse2(string)	

	/*----- 2. In sample 2, compute gamma_sel and its variance ----*/
	di as txt "Step 2: In sample 2, get consistent estimate of gamma"

	// using sample 2, run dsregress of dsregress on x3_sel
	di as txt "dsregress {it:`x2_var'}"
	qui dsregress `x2_var' `x3_sel' if `touse2', controls(`x3_nosel')

	di
end
					//----------------------------//
					//  Predict X2
					//----------------------------//
program PredictX2, rclass
	syntax , x2_var(string)	///
		touse1(string)	///
		x2hat(string)	///
		x3_sel(string)	///
		x3e_vars(string)

	/*----- 3. In sample 1, predict x2 ---------------------------*/	

	di as txt "Step 3: In sample 1, predict x2"

	// using sample 1, get x3_sel * gamma_sel
	di as txt "predict {it:`x2_var'}"
	qui predict double `x2hat' if `touse1', xb

	// get x3e not selected
	local x3_sel : list uniq x3_sel
	local x3e_nosel : list x3e_vars - x3_sel

	di

	ret local x3_sel `x3_sel'
	ret local x3e_nosel `x3e_nosel'
end
					//----------------------------//
					//   compute beta
					//----------------------------//
program ComputeBeta, rclass
	syntax,	depvar(string)		///
		x1_vars(string)		///
		x2hat(string)		///
		x3i_vars(string)	///
		x3_sel(string)		///
		x3e_nosel(string)	///
		touse1(string)		///
		lassoopts(string)

	/*----- 4. In sample 1, run double selection lasso for main model --*/

	di "Step 4: In sample 1, run double selection for the main model"

	// selected the controls for y, x1, x3i, and x3_sel

	local vars `depvar' `x1_vars' `x3i_vars' `x3_sel'
	local vars : list uniq vars

	foreach var in `vars' {
		di as txt "Estimating lasso for {it:`var'}"
		tempvar tmp
		qui gen double `tmp' = `var' if `touse1'
		qui lasso linear `tmp' `x3e_nosel' if `touse1', `lassoopts'
		local controls_sel `controls_sel' `e(allvars_sel)'
	}
	local controls_sel : list uniq controls_sel

	// OLS y on varsofinterest and selected controls
	di as txt "Estimating regress for {it:`depvar'} on selected vars"
	qui regress `depvar' `x1_vars' `x2hat' `x3i_vars' `controls_sel' ///
		if `touse1', vce(robust)

	local b2 = _b[`x2hat']
	local xa `x1_vars' `x2hat' `x3i_vars' `controls_sel'
	
	ret local controls_sel `controls_sel'
	ret local b2 `b2'
	ret local xa `xa'
end
					//----------------------------//
					//  adjust variance
					//----------------------------//
program AdjustVar, eclass
	syntax, vg(string)	///
		vb(string)	///
		b2(string)	///
		x3_sel(string)	///
		xa(string)	///
		touse1(string)	///
		touse2(string)	///
		adjust(string)
	
	if (`"`adjust'"' == "noadjust") {
		exit
		// NotReached
	}
	
	di
	di as txt "Step 5: adjust variance"
	mata: adjust_var(`"`vg'"', `"`vb'"', `b2', `"`x3_sel'"', `"`xa'"', ///
		`"`touse1'"', `"`touse2'"')
end

/*----------------------------------------------------------------------------*/ 
// 		Mata functions
/*----------------------------------------------------------------------------*/ 
mata:
mata set matastrict on

void adjust_var(			///
	string scalar	_vg,		///
	string scalar	_vb,		///
	real scalar	_b2,		///
	string scalar	_x3_sel,	///
	string scalar	_xa,		///
	string scalar	_touse1,	///
	string scalar	_touse2)
{
	real scalar	m, n, k
	real matrix	XA, Phi_XA_XA, X3s, Phi_XA_X3s, V1, Vg, V
	string scalar	Vs

	n = sum(st_data(., _touse1, _touse1))
	m = sum(st_data(., _touse2, _touse2))
	Xa = (st_data(., _xa, _touse1), J(n, 1, 1))
	X3s = st_data(., _x3_sel, _touse1)

	k = n/m
	Phi_XA_XA = invsym(cross(Xa, Xa)/n)
	Phi_XA_X3s = cross(Xa, X3s)/n

	V1 = st_matrix(_vb)

	Vg = st_matrix(_vg)
	Psi2 = Phi_XA_X3s*(m*Vg*_b2^2)*(Phi_XA_X3s')
	V2 = k* Phi_XA_XA*Psi2*Phi_XA_XA/n

	V = V1 + V2
	Vs = st_tempname()
	st_matrix(Vs, V)

	stata(sprintf("eret repost V=%s", Vs))
}


end
