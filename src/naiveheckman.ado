*!version 1.0.0  17jun2020
program naiveheckman
	version 16.0

	if (replay()) {
		Display `0'
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
	Compute, y1(`r(y1)')		///
		xvars(`r(xvars)')	///
		y2(`r(y2)')		///
		zvars(`r(zvars)')	///
		touse(`touse')
end
					//----------------------------//
					//  compute
					//----------------------------//
program Compute, eclass
	syntax, y1(string)	///
		xvars(string)	///
		y2(string)	///
		zvars(string)	///
		touse(string)
	
	/*---- step 1: lasso probit y2 on zvars ----*/
	di
	di as txt "step 1: lasso probit {it:`y2'}"
	qui lasso probit `y2' `xvars' `zvars' if `touse', sel(plugin)
	local zvars_sel `e(allvars_sel)'
	local zvars_nosel : list zvars - zvars_sel

	_dsheckman_getvardim, xvars(`xvars') 	///
		zvars_sel(`zvars_sel') zvars_nosel(`zvars_nosel')
	local p = r(p)
	local k1 = r(k1)
	local k = r(k)

	/*---- step 2: post-lasso-probit y2----*/
	di as txt "step 2: post-lasso probit {it:`y2'}"
	qui probit `y2' `zvars_sel' if `touse'

	tempname Vp
	mat `Vp' = e(V)

	tempvar zr 
	predict double `zr', xb

	/*---- step 3: compute lambda ----*/
	di as txt "step 3: compute lambda"
	tempvar lambda
	qui gen double `lambda' = normalden(`zr')/normal(`zr')

	tempvar touse_main
	qui gen byte `touse_main' = ( (`y2' == 1) & `touse' )

	/*---- step 4: ols of y1 on xvars and lambda ----*/
	di as txt "step 4: ols of {it:`y1'} on {it:`xvars'} and lambda"
	qui regress `y1' `xvars' `lambda' if `touse_main'

	tempvar u1
	qui predict double `u1' if `touse_main', res
	local bm = _b[`lambda']

	/*---- adjust variance ----*/
	tempname Vs
	adjust_heckman_V, 			///
		touse_main(`touse_main')	///
		zr(`zr')			///
		u1(`u1')			///
		bm(`bm')			///
		xvars(`xvars')			///
		zvars(`zvars_sel')		///
		lambda(`lambda')		///
		vp(`Vp')			///
		vs(`Vs')

	/*---- post result ----*/
	tempname b 
	mat `b' = e(b)
	local bs `xvars' lambda _cons
	mat colname `b' = `bs'
	mat colname `Vs' = `bs'
	mat rowname `Vs' = `bs'

	qui count if `touse' 
	local N = r(N)

	qui count if `touse' & `y2' == 1
	local N_sel = r(N)

	qui count if `touse' & `y2' != 1
	local N_nonsel = r(N)

	eret post `b' `Vs', esample(`touse') buildfvinfo

	eret scalar N = `N'
	eret scalar N_sel = `N_sel'
	eret scalar N_nonsel = `N_nonsel'
	eret local title "High-dimensional Heckman naive"
	eret local cmd naiveheckman
	eret scalar p = `p'
	eret scalar k1 = `k1'
	eret scalar k = `k'

	/*---- display ----*/
	_display_heckman `0'
end
