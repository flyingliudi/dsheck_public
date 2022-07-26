*!version 1.0.0  27jul2021
program poheckman, eclass
	version 16.0

        _vce_parserun poheckman : `0'
        if "`s(exit)'" != "" {
                ereturn local cmdline `"poheckman `0'"'
                exit
        }

	if (replay()) {
		_display_heckman `0'
	}
	else {
		Estimate `0'
	}
end
					//----------------------------//
					//  Estimate
					//----------------------------//
program Estimate 

	/*---- parse syntax ----*/

	tempvar touse
	_parse_heckman `touse': `0'

	/*---- compute ----*/
	Compute, y1(`r(y1)')				///
		xvars(`r(xvars)')			///
		y2(`r(y2)')				///
		zvars(`r(zvars)')			///
		mlauto(`r(mlauto)')			///
		selvars(`r(selvars)')			///
		touse(`touse')				///
		selopts(`r(selopts)')			///
		`r(qui)'				///
		est_probit(`r(est_probit)')		///
		esample_probit(`r(esample_probit)')	///
		esample_main(`r(esample_main)')
end
					//----------------------------//
					//  compute
					//----------------------------//
program Compute, eclass
	syntax, y1(string)		///
		xvars(string)		///
		y2(string)		///
		zvars(string)		///
		touse(string)		///
		[mlauto(string) 	///
		selopts(string)		///
		selvars(string) 	///
		esample_probit(string)	///
		esample_main(string)	///
		est_probit(string)	///
		qui]

	di

	/* ----------------------------------------------------------- */
	// step 1: lasso probit to select vars

	di as txt "step 1: lasso probit to select vars"
	qui lasso probit `y2' `xvars' `zvars' if `touse', `selopts'
	local allvars `e(allvars)'
	local allvars_sel `e(allvars_sel)'
	local allvars_nosel : list allvars - allvars_sel

	ExtractCommon, zvars_sel(`allvars_sel') zvars_nosel(`allvars_nosel')
	local allvars_sel `r(zvars_sel)'
	local allvars_nosel `r(zvars_nosel)'

	local p : list sizeof allvars
	local k1 : list sizeof xvars
	local k : list sizeof allvars_sel

	/* ----------------------------------------------------------- */
	// step 2: post-lasso-probit estimation 

	di as txt "step 2: probit of y2 on selected zvars"
	qui probit `y2' `allvars_sel' if `touse'
	local N_total = e(N)

	tempname b_probit
	mat `b_probit' = e(b)
					// mytouse
	tempvar mytouse
	qui gen byte `mytouse' = e(sample)

					// esample probit
	if (`"`esample_probit'"' != "") {
		qui gen byte `esample_probit' = e(sample)
	}

					// verbose
	`qui' probit

					// est_probit
	if (`"`est_probit'"' != "") {
		qui est store `est_probit'
	}

	/* ----------------------------------------------------------- */
	// step 3: compute labmda

	di as txt "step 3: compute lambda"

					// zr
	tempvar zr
	qui _predict double `zr' if `touse', xb

					// lambda
	tempvar lambda
	qui gen double `lambda' = normalden(`zr')/normal(`zr')

	/* ----------------------------------------------------------- */
	// step 4: partialling out lambda
	di as txt "step 4: partialling-out estimation"

	tempvar touse_main
	qui gen byte `touse_main' = `touse' & `y2' == 1
	
	local othervars : list allvars - xvars
	qui poregress `y1' `xvars', controls(`lambda' `othervars')
	local N_sel = e(N)

	tempname b vs
	mat `b' = e(b)
	mat `vs' =e(V)

	/*---- post result ----*/
	PostResult , 				///
		b(`b')				///
		vs(`vs')			///
		xvars(`xvars')			///
		zvars(`zvars')			///
		y1(`y1')			///
		y2(`y2')			///
		p(`p')				///
		k1(`k1')			///
		myk(`k')			///
		b_probit(`b_probit')		///
		n_total(`N_total')		///
		n_sel(`N_sel')			///
		mytouse(`mytouse')		///
		allvars_sel(`allvars_sel')

	/*---- display ----*/
	_display_heckman `0'
end

program PostResult, eclass
	syntax, b(string)		///
		vs(string)		///
		xvars(string)		///
		mytouse(string)		///
		y2(string)		///
		y1(string)		///
		myk(string)		///
		k1(string)		///
		p(string)		///
		zvars(string)		///
		b_probit(string)	///
		n_sel(string)		///
		n_total(string)		///
		allvars_sel(string)


	local N = `n_total'
	local N_sel = `n_sel'

	local N_nonsel = `N' - `N_sel'

	eret post `b' `vs', esample(`mytouse') buildfvinfo

	eret scalar N = `N'
	eret scalar N_sel = `N_sel'
	eret scalar N_nonsel = `N_nonsel'
	eret scalar k = `myk'
	eret scalar k1 = `k1'
	eret scalar p = `p'
	eret local varsofinterest `xvars'
	eret local controls_all `zvars'
	eret local controls_sel `allvars_sel'
	eret matrix b_probit = `b_probit'
	eret local title "Partialling-out Heckman"
	eret local cmd poheckman
end
						//---------------------------//
						//	get fv vars
						//---------------------------//
program GetFvVars, rclass
	syntax [anything(name=vars)]

	if (`"`vars'"' == "") {
		exit
		// NotReached
	}

	local vars =  ustrregexra(`"`vars'"', "bn\.", ".")	
	local vars =  ustrregexra(`"`vars'"', "b\.", ".")	
	local vars : list uniq vars

	fvexpand `vars'
	local vars `r(varlist)'

	ret local varlist `vars'
end
					//----------------------------//
					// extract common
					//----------------------------//
program ExtractCommon, rclass
	syntax [, zvars_sel(string)	///
		zvars_nosel(string) ]
	
	local myzvars_sel `zvars_sel'
	local myzvars_nosel `zvars_nosel'
	
	foreach var of local zvars_sel {
		_ms_parse_parts `var'
		if (`"`r(type)'"' == "factor") {
			local base_sel `base_sel' `r(name)'
		}
	}

	foreach var of local zvars_nosel {
		_ms_parse_parts `var'

		if (`"`r(type)'"' == "factor") {
			local aa `r(name)'
			local in : list aa in base_sel

			if (`in') {
				local myzvars_sel `myzvars_sel' `var'
				local myzvars_nosel : list myzvars_nosel - var
			}
		}
	}

	fvexpand `myzvars_sel'
	local myzvars_sel `r(varlist)'

	fvexpand `myzvars_nosel'
	local myzvars_nosel `r(varlist)'
	
	ret local zvars_sel `myzvars_sel'
	ret local zvars_nosel `myzvars_nosel'
end
