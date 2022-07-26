*!version 1.0.0  16jun2020
program dsprobit
	
	version 16.0

	if (replay()) {
		Display `0'
	}
	else {
		Estimate `0'
	}
end
					//----------------------------//
					// display
					//----------------------------//
program Display
	syntax [, *]

	local title as txt "`e(title)'"

	local col = 39
						//  nobs
	local nobs _col(`col') as txt "Number of obs" _col(67) "="	///
		_col(69) as res %10.0fc e(N)

						//  number of controls
	local k_controls _col(`col') as txt "Number of controls"	///
		_col(67) "=" _col(69) as res %10.0fc e(k_controls)

						//  number of selected controls
	local k_controls_sel _col(`col') as txt	///
		"Number of selected controls" ///
		_col(67) "=" _col(69) as res %10.0fc e(k_controls_sel)


	di 
	di `title' `nobs'
	di `k_controls'
	di `k_controls_sel'
	_coef_table, `options'
end
					//----------------------------//
					//  estimate
					//----------------------------//
program Estimate
	syntax anything(name=eq) [if] [in] , controls(passthru) [nooptimal]
	
	marksample touse
					// parse syntax
	_dslasso_parse `eq', `controls' model(probit)
	local depvar `s(depvar)'
	local dvars `s(Dvars)'
	local xvars `s(Xvars)'
	local dvars_omit `s(dvars_omit)'
	local dvars_full `s(dvars_full)'
					//  markout
	markout `touse' `depvar' `dvars' `xvars'
	
					//  compute 
	Compute, depvar(`depvar')	///
		dvars(`dvars')		///
		xvars(`xvars')		///
		touse(`touse')		///
		optimal(`optimal')

  	_dslasso_stripe_b , 			///
		dvars_omit(`dvars_omit') 	///
		dvars_full(`dvars_full' _cons) 	///
                lasso_noneed(0)

					//  display
	Display
end

					//----------------------------//
					//  compute
					//----------------------------//
program Compute, eclass
	syntax, depvar(string)	///
		dvars(string)	///
		xvars(string)	///
		touse(string)	///
		[optimal(string)]
	
	/*---- step 1: post-lasso-probit on d and x ----*/
	di 
	di as txt "Step 1: post-lasso-probit"
	fvexpand `xvars'
	local xvars `r(varlist)'
					// lasso 
	qui lasso probit `depvar' `dvars' `xvars' if `touse', sel(plugin)
	local x1_sel `e(allvars_sel)'
	local x1_sel : list x1_sel - dvars


					// post-lasso
	qui probit `depvar' `dvars' `x1_sel' if `touse'

					// get weight
	tempvar xb	
	predict double `xb', xb

	tempvar wvar
	qui gen double `wvar' = normalden(`xb')

	tempvar pr 

	if (`"`optimal'"' != "nooptimal") {
		qui gen double `pr' = normal(`xb')
		qui replace `pr' = `pr'*(1 - `pr')
	}
	else {
		qui gen double `pr' = 1
	}

	tempvar fvar
	qui gen double `fvar' = `wvar'^2/`pr'

	/*---- step 2: lasso ols of d on x with iweight fvar ----*/
	di as txt "Step 2: ols lasso"
	qui foreach var of local dvars {
		tempvar tmp
		qui gen double `tmp' = `var'
		lasso linear `tmp' `xvars' [iweight = `fvar'] if `touse', ///
			sel(plugin) depname(`var')
		local x2_sel `x2_sel' `e(allvars_sel)'
	}
	local x2_sel : list uniq x2_sel


	/*---- step 3: weighted probit with union------*/
	di as txt "Step 3: reweighted probit"
	local x_sel : list x1_sel | x2_sel
	tempvar fvar2
	gen double `fvar2' = `wvar'/`pr'
	qui probit `depvar' `dvars' `x_sel' [iweight=`fvar2'] if `touse',  ///
		vce(robust)
	local nobs = e(N)

	tempvar esample
	qui gen byte `esample' = e(sample)
	
	/*----- step 4: post e(b) and e(V) ----*/
	fvexpand `dvars'
	local myvars `r(varlist)'
	local k : list sizeof myvars
	tempname b V
	mata: get_b(`"`depvar'"', `"`myvars'"', `k', "`b'", "`V'")

	eret post `b' `V', esample(`esample')

	eret local title Double selection probit
	eret local controls `xvars'
	eret local varsofinterest `dvars'
	eret local controls_sel `x_sel'
	eret local depvar `depvar'
	eret local vce robust
	eret local vcetype Robust
	eret scalar k_controls = `: list sizeof xvars'
	eret scalar k_controls_sel = `: list sizeof x_sel'
	eret scalar N = `nobs'
	eret local cmd dsprobit
end

mata:
mata set matastrict on

void get_b(			///
	string scalar depvar,	///
	string scalar dvars,	///
	real scalar _k,		///
	string scalar _myb,	///
	string scalar _myV)
{
	real matrix	b, V, p
	string matrix	bs

	b = st_matrix("e(b)")
	V = st_matrix("e(V)")

	p = (1.._k), length(b)
	
	b = b[1, p]
	V = V[p, p]

	bs = J(length(p),2, "")
	bs[., 1] = J(length(p), 1, depvar)
	bs[., 2] = (tokens(dvars), "_cons")'

	st_matrix(_myb, b)
	st_matrix(_myV, V)
	st_matrixcolstripe(_myb, bs)
	st_matrixcolstripe(_myV, bs)
	st_matrixrowstripe(_myV, bs)
}


end
