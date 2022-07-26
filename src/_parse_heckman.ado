*!version 1.0.0  27jul2021
program _parse_heckman
					// parse
	Parse `0'
					// check
	Check
end
						//------------------------//
						//  parse
						//------------------------//
program Parse, rclass
	_on_colon_parse `0'
	local mytouse `s(before)'
	local 0 `s(after)'

	syntax anything(name=eq) 		///
		[if] [in],			///
		SELection(string)		///
		[selvars(passthru)		///
		AUTOmatic			///
		selopts(passthru)		///
		noverbose			///
		esample_probit(passthru)	///
		esample_main(passthru)		///
		est_probit(passthru)		///
		nocheck]
	
	marksample touse
	qui gen byte `mytouse' = `touse'

					//  parse main equation
	ParseEq `eq', touse(`touse')
	return add

					//  parse selection equation
	ParseSel `selection'
	return add

					// parse selection equation varsinterest
	ParseSelvars, `selvars' `automatic'
	return add

	ParseSelopts, `selopts'
	return add

	ParseVerbose, `verbose'
	return add

	ParseEsampleProbit, `esample_probit'
	return add

	ParseEsampleMain, `esample_main'
	return add

	ParseEstProbit, `est_probit'
	return add

	ParseCheck, `check'
	return add
end
					//----------------------------//
					// parse check
					//----------------------------//
program ParseCheck, rclass
	syntax [, nocheck]

	if (`"`check'"' == "nocheck") {
		local check = 0
	}
	else {
		local check = 1
	}

	ret scalar check = `check'
end
					//----------------------------//
					// Parse verbose	
					//----------------------------//
program ParseVerbose, rclass
	syntax [, noverbose]

	if (`"`verbose'"' == "noverbose") {
		local qui qui
	}

	ret local qui `qui'
end
					//----------------------------//
					// parse selection options	
					//----------------------------//
program ParseSelopts, rclass
	syntax [, selopts(string)]

	if (`"`selopts'"' == "") {
		local opts plugin
	}
	else if (`"`selopts'"' != "bic" &	///
		`"`selopts'"' != "cv" &		///
		`"`selopts'"' != "plugin" &	///
		`"`selopts'"' != "adaptive" &	///
		`"`selopts'"' != "") {
		di as err "option {bf:selopts()} allows only one of "	///
			"{bf:bic}, {bf:cv}, {bf:plugin}, or {bf:adaptive}"
		exit 198
	}
	else {
		local opts `selopts'
	}

	return local selopts selection(`opts')
end
					//----------------------------//
					//  parse main equation
					//----------------------------//
program ParseEq, rclass
	syntax varlist(numeric fv), touse(string)

	gettoken y1 xvars : varlist

	_rmcoll `xvars' if `touse', expand
	local xvars `r(varlist)'

	ret local y1 `y1'	
	ret local xvars `xvars'
end
					//----------------------------//
					// parse selection equation
					//----------------------------//
program ParseSel, rclass
	syntax anything(name=eq equalok)

	gettoken y2 next : eq, parse("=")	
	gettoken equal zvars: next, parse("=")

	if (`"`equal'"' != "=") {
		di as error "misspecified {bf:selection()}"
		di "{p 4 4 2}"
		di "the syntax for {bf:selection()} is "	///
			"{bf:selection({it:y2} = {it:zvars})}"
		di "{p_end}"
	}

	fvexpand `zvars'
	local zvars `r(varlist)'

	ret local y2 `y2'
	ret local zvars `zvars'
end
					//----------------------------//
					//  parse sel vars of interest
					//----------------------------//
program ParseSelvars, rclass
	syntax [, selvars(varlist numeric fv) 	///
		automatic]
	
	if (`"`selvars'"' == "" | `"`automatic'"' != "") {
		return local mlauto mlauto
		exit
		// NotReached
	}

	fvexpand `selvars'
	local selvars `r(varlist)'
	return local selvars `selvars'
end
						//------------------------//
						// CheckZvars
						//------------------------//
program Check
	if (r(check) == 0) {
		exit
		// NotReached
	}

	local xvars `r(xvars)'
	local zvars `r(zvars)'
	local inzvars : list xvars in zvars	

	if (!`inzvars') {
		di "{p 0 4 2}"
		di as err "option {bf:sel()} must contain "	///
			"variables in the main equation, which are {bf:`xvars'}"
		di "{p_end}"
		exit 198
	}


	local selvars `r(selvars)'
	if (`"`selvars'"' != "") {
		local inselvars : list xvars in selvars
		if (!`inselvars') {
			di "{p 0 4 2}"
			di as err "option {bf:selvars()} must contain "	///
				"variables in the main equation, which " ///
				"are {bf:`xvars'}"
			di "{p_end}"
			exit 198
		}
	}
	
end
					//----------------------------//
					// esample probit	
					//----------------------------//
program ParseEsampleProbit, rclass
	syntax [, esample_probit(namelist max=1)]

	if (`"`esample_probit'"' == "") {
		exit 
		// NotReached
	}

	return local esample_probit `esample_probit'
end

					//----------------------------//
					// esample main	
					//----------------------------//
program ParseEsampleMain, rclass
	syntax [, esample_main(namelist max=1)]

	if (`"`esample_main'"' == "") {
		exit 
		// NotReached
	}

	return local esample_main `esample_main'
end
					//----------------------------//
					// parse est probit	
					//----------------------------//
program ParseEstProbit, rclass
	syntax [, est_probit(string) ]

	if (`"`est_probit'"' == "") {
		exit
		// NotReached
	}

	ret local est_probit `est_probit'
end
