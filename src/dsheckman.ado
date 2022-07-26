*!version 1.0.0  27jul2021
program dsheckman 
	version 16.0

        _vce_parserun dsheckman : `0'
        if "`s(exit)'" != "" {
                ereturn local cmdline `"dsheckman `0'"'
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
	if (`"`mlauto'"' != "") {
		/*---- step 1: lasso probit to select vars ----*/
		di as txt "step 1: lasso probit to select vars"
		qui lasso probit `y2' `xvars' `zvars' if `touse', `selopts'
		local zvars_sel `e(allvars_sel)'
		local zvars_sel : list zvars_sel | xvars

		GetFvVars `zvars_sel'
		local zvars_sel `r(varlist)'

		local allvars `xvars' `zvars'
		GetFvVars `allvars'
		local allvars `r(varlist)'

		local zvars_nosel : list allvars - zvars_sel
		GetFvVars `zvars_nosel'
		local zvars_nosel `r(varlist)'

		ExtractCommon, zvars_sel(`zvars_sel') zvars_nosel(`zvars_nosel')
		local zvars_sel `r(zvars_sel)'
		local zvars_nosel `r(zvars_nosel)'
	}
	else {
		di as text "step 1: set {it:varsofinterest} in "	///
			"selection equation"
		local zvars_sel `selvars'
		local zvars_nosel : list zvars - zvars_sel
	}

	_dsheckman_getvardim, xvars(`xvars') 	///
		zvars_sel(`zvars_sel') zvars_nosel(`zvars_nosel')
	local p = r(p)
	local k1 = r(k1)
	local k = r(k)

	/*---- step 2: dsprobit of y2 on selected zvars ----*/
	di as txt "step 2: dsprobit of y2 on selected zvars"
	qui dsprobit `y2' `zvars_sel' if `touse', controls(`zvars_nosel')

	local N_total = e(N)
	local k_controls_probit = e(k_controls)
	local k_controls_sel_probit = e(k_controls_sel)

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
	`qui' dsprobit

					// est_probit
	if (`"`est_probit'"' != "") {
		qui est store `est_probit'
	}

	// get zr
	tempvar zr
	qui _predict double `zr' if `touse', xb

	// get variance of r
	tempname vp
	mat `vp' = e(V)

	/*---- step 3: compute lambda ----*/
	di as txt "step 3: compute lambda"
	tempvar lambda
	qui gen double `lambda' = normalden(`zr')/normal(`zr')

	/*---- step 4: dsregress y1 on xvars, lambda with controls ----*/
	di as txt "step 4: dsregress y1 on xvars, lambda with controls"
	tempvar touse_main
	qui gen byte `touse_main' = `touse' & `y2' == 1

	qui dsregress `y1' `xvars' `lambda' if `touse_main', 	///
		controls(`zvars_nosel') `selopts'
	local controls_sel `e(controls_sel)'
	local N_sel = e(N)
	local k_controls_main = `k_controls_probit'
	local k_controls_sel_main = e(k_controls_sel)

	local bm = _b[`lambda']

	if (`"`esample_main'"' != "") {
		qui gen byte `esample_main' = e(sample)
	}

	qui regress `y1' `xvars' `lambda' `controls_sel' if `touse_main'

	// get u1
	tempvar u1
	qui predict double `u1' if `touse_main', res

	// get b 
	tempname b 
	mat `b' = e(b)

	/*---- step 5: adjust variance ----*/
	tempname vs
	adjust_heckman_V, 			///
		touse_main(`touse_main')	///
		zr(`zr')			///
		u1(`u1')			///
		bm(`bm')			///
		xvars(`xvars')			///
		zvars(`zvars_sel')		///
		lambda(`lambda')		///
		controls_sel(`controls_sel')	///
		vp(`vp')			///
		vs(`vs')

	/*---- post result ----*/
	PostResult , b(`b')					///
		vs(`vs')					///
		xvars(`xvars')					///
		zvars(`zvars')					///
		y1(`y1')					///
		y2(`y2')					///
		p(`p')						///
		k1(`k1')					///
		k(`k')						///
		zvars_sel(`zvars_sel')				///
		b_probit(`b_probit')				///
		n_total(`N_total')				///
		n_sel(`N_sel')					///
		mytouse(`mytouse')				///
		k_controls_sel_main(`k_controls_sel_main')	///
		k_controls_main(`k_controls_main')		///
		k_controls_sel_probit(`k_controls_sel_probit')	///
		k_controls_probit(`k_controls_probit')	

	/*---- display ----*/
	_display_heckman `0'
end

program PostResult, eclass
	syntax, b(string)			///
		vs(string)			///
		xvars(string)			///
		mytouse(string)			///
		y2(string)			///
		y1(string)			///
		k(string)			///
		k1(string)			///
		p(string)			///
		zvars_sel(string)		///
		zvars(string)			///
		b_probit(string)		///
		n_sel(string)			///
		n_total(string)			///
		k_controls_sel_main(string)	///
		k_controls_main(string)		///
		k_controls_sel_probit(string)	///
		k_controls_probit(string)	

	mata: permute_b_V("`b'", "`vs'", "`xvars' lambda")

	local bs `xvars' lambda 
	mat colname `b' = `bs'
	local k : list sizeof bs
	forvalues i = 1/`k' {
		local coleq `coleq' `y1'
	}
	mat coleq `b' = `coleq'
	mat colname `vs' = `bs'
	mat rowname `vs' = `bs'
	mat coleq `vs' = `coleq'
	mat roweq `vs' = `coleq'

	local N = `n_total'
	local N_sel = `n_sel'

	local N_nonsel = `N' - `N_sel'

	eret post `b' `vs', esample(`mytouse') buildfvinfo

	eret scalar N = `N'
	eret scalar N_sel = `N_sel'
	eret scalar N_nonsel = `N_nonsel'
	eret scalar k = `k'
	eret scalar k1 = `k1'
	eret scalar p = `p'
	eret local varsofinterest `xvars'
	eret local controls_all `zvars'
	eret matrix b_probit = `b_probit'
	eret local controls_sel `zvars_sel'
	eret scalar k_controls_sel_main = `k_controls_sel_main'
	eret scalar k_controls_main = `k_controls_main'
	eret scalar k_controls_sel_probit = `k_controls_sel_probit'
	eret scalar k_controls_probit = `k_controls_probit'
	eret local title "Double-selection-lasso Heckman"
	eret local cmd dsheckman
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


mata:
mata set matastrict on

void permute_b_V(		///
	string scalar	_bs,	///
	string scalar	_vs,	///
	string scalar	_vars)
{
	real matrix b, V, p

	b = st_matrix(_bs)
	V = st_matrix(_vs)

	_vars = tokens(_vars)
	p = 1..length(_vars)
	b = b[1, p]
	V = V[p, p]

	st_matrix(_bs, b)
	st_matrix(_vs, V)
}

end
