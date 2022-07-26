*!version 1.0.0  17jun2020
program myheckman 
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
	
	/*---- step 1: probit y2----*/
	di
	di as txt "step 1: probit {it:`y2'}"
	qui probit `y2' `zvars' if `touse'

	tempname Vp
	mat `Vp' = e(V)

	tempvar zr 
	predict double `zr', xb

	/*---- step 2: ols of y1 on xvars and lambda ----*/
	di as txt "step 2: ols of {it:`y1'} on {it:`xvars'} and lambda"

	tempvar lambda
	qui gen double `lambda' = normalden(`zr')/normal(`zr')

	tempvar touse_main
	qui gen byte `touse_main' = ( (`y2' == 1) & `touse' )

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
		zvars(`zvars')			///
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

	eret post `b' `Vs', esample(`touse')

	eret scalar N = `N'
	eret scalar N_sel = `N_sel'
	eret scalar N_nonsel = `N_nonsel'
	eret local title "Heckman two-step"
	eret local cmd myheckman

	_display_heckman `0'
end
