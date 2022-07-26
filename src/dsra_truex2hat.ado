*! version 1.0.0  06apr2020
/*-----------------------------------------------------------------------------
	infeasible double selection lasso for matching based estimation using
	true x2hat
-----------------------------------------------------------------------------*/
program dsra_truex2hat, eclass
	version 16.0	
	
	/*----- preserve START -----*/
	preserve

	/*----- parse syntax ----*/
	tempvar touse1 touse2
	_dsra_parse_syntax, touse1(`touse1') touse2(`touse2') : `0'

	/*----- compute -----*/
	Compute, depvar(`r(depvar)')		///
		x1_vars(`r(x1_vars)')		///
		x2_vars(`r(x2_vars)')		///
		x3i_vars(`r(x3i_vars)')		///
		x3e_vars(`r(x3e_vars)')		///
		touse1(`r(touse1)')		///
		touse2(`r(touse2)')		///
		lassoopts(`r(lassoopts)')	///
		true_gamma(`r(true_gamma)')	///
		`r(allgamma)'
	
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
		x2_vars(string)		///
		x3i_vars(string)	///
		x3e_vars(string)	///
		touse1(string)		///
		touse2(string)		///
		true_gamma(passthru)	///
		[lassoopts(string)	///
		allgamma]
	
	/*------- 1: In sample 2, lasso x2 on x3i and x3e  ---*/
	foreach var of local x2_vars {
		
		// lasso x2 on x3
		di as txt "Estimating lasso for {it:`var'}"
		qui lasso linear `var' `x3i_vars' `x3e_vars' 	///
			if `touse2', `lassoopts'

		// get selected vars for x2
		local x2_sel `x2_sel' `e(allvars_sel)'

		// get prediction for x2
		tempvar tmp

		GetTrueX2hat, var(`tmp') touse1(`touse1') ///
			`true_gamma' allvars_sel(`e(allvars_sel)') `allgamma'

		local x2hat `x2hat' `tmp'
	}

	/*-------- 2: get x3e not selected --------*/
	local x2_sel : list uniq x2_sel
	local x3e_nosel : list x3e_vars - x2_sel

	/*-------- 3: In sample 1, do lasso y, x1, x3i on x3e_nosel ---*/
	foreach var in `depvar' `x1_vars' `x3i_vars' {
		di as txt "Estimating lasso for {it:`var'}"
		qui lasso linear `var' `x3e_nosel' if `touse1', `lassoopts'
		local controls_sel `controls_sel' `e(allvars_sel)'
	}
	local controls_sel : list uniq controls_sel

	/*- 4: In sample 1, OLS of y on x1, x2hat, x3i and controls_sel --*/
	di as txt "Estimating regress for {it:`depvar'} on selected vars"
	qui regress `depvar' `x1_vars' `x2hat' `x3i_vars' `controls_sel' ///
		if `touse1', vce(robust)

	/*-- 5. post result -----*/
	local x3e_sel `controls_sel' `x2_sel'
	PostResult, 				///
		depvar(`depvar')		///
		x1_vars(`x1_vars')		///
		x2_vars(`x2_vars')		///
		x3i_vars(`x3i_vars')		///
		x3e_vars(`x3e_vars')		///
		touse1(`touse1') 		///
		touse2(`touse2') 		///
		x3e_sel(`x3e_sel')		///
		`allgamma'
end
					//----------------------------//
					//  post coeff
					//----------------------------//
program PostResult, eclass
	syntax , depvar(string)		///
		x1_vars(string)		///
		x2_vars(string)		///
		x3i_vars(string)	///
		x3e_vars(string)	///
		touse1(string)		///
		touse2(string)		///
		[x3e_sel(string)	///
		allgamma]

	local vars `x1_vars' `x2_vars' `x3i_vars'
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

	eret local title "Infeasible dsra `allgamma'"
	eret local cmd dsra_truex2hat`allgamma'
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

	eret post `b' `V'
end

					//----------------------------//
					// Get TrueX2hat
					//----------------------------//
program GetTrueX2hat
	syntax, var(string)		///
		touse1(string)		///
		true_gamma(string)	///
		allvars_sel(string)	///
		[allgamma]
	
	tempname b
	mat `b' = `true_gamma'

	local truevars `:colname `b''

	if (`"`allgamma'"' == "") {
		local allvars_sel : list allvars_sel & truevars
		local allvars_sel `allvars_sel' _cons
	}
	else {
		local allvars_sel `truevars'
	}

	qui gen double `var' =  0 if `touse1'

	foreach x of local allvars_sel {

		if (`"`x'"' == "_cons") {
			local xvar = 1
		}
		else {
			local xvar `x'
		}

		qui replace `var' = `var' + `xvar'*`b'[1, "`x'"] if `touse1'
	}
end
