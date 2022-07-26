*!version 1.0.0  06apr2020
/*-------------------------------------------------------------------
	naive regression adjustment 
-------------------------------------------------------------------*/ 
program naivera
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
		touse2(`r(touse2)')		
	
	/*----- preserve END -----*/
	restore

	/*--- Display -----*/
	Display

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
	
	/*------- 1: In sample 2, lasso x2 on x3i and x3e  ---*/
	foreach var of local x2_vars {
		
		// lasso x2 on x3
		di as txt "Estimating lasso for {it:`var'}"
		qui lasso linear `var' `x3i_vars' `x3e_vars' 	///
			if `touse2', `lassoopts'

		// get prediction for x2
		tempvar tmp
		predict double `tmp' if `touse1', postsel xb
		local x2hat `x2hat' `tmp'
	}

	/*- 2: In sample 1, OLS of y on x1, x2hat, x3i and controls_sel --*/
	di as txt "Estimating regress for {it:`depvar'} on vars"
	qui regress `depvar' `x1_vars' `x2hat' `x3i_vars' if `touse1', 	///
		vce(robust)

	/*-- 3. post result -----*/
	PostResult, 				///
		depvar(`depvar')		///
		x1_vars(`x1_vars')		///
		x2_vars(`x2_vars')		///
		x3i_vars(`x3i_vars')		///
		x3e_vars(`x3e_vars')		///
		touse1(`touse1') 		///
		touse2(`touse2') 		
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
		touse2(string)		

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

	eret scalar k_controls = `:list sizeof x3e_vars'

	eret local title "Post-lasso regression adjustment"
	eret local cmd naivera 
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
					// Display
					//----------------------------//
program Display
	Header
	_coef_table
end
					//----------------------------//
					// head
					//----------------------------//
program Header
//--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----
						// title
	local title _n as txt `"`e(title)'"'

	local col1 = 40
	local col2 = 68

						// Number of obs
	local n_obs1 as txt _col(`col1') `"Number of obs in sample 1"' 	///
		_col(`col2') "=" _col(69) as res %10.0fc e(N_1)

	local n_obs2 as txt _col(`col1') `"Number of obs in sample 2"' 	///
		_col(`col2') "=" _col(69) as res %10.0fc e(N_2)

	local k_controls as txt _col(`col1') `"Number of controls"' 	///
		_col(`col2') "=" _col(69) as res %10.0fc e(k_controls)


	di `title' `n_obs1'
	di `n_obs2'
	di `k_controls'
	di 
end
