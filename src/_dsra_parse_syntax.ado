*!version 1.0.0  06apr2020
/*
syntax:

dsra depvar x1_varlist (x2_varlist = x3i_varlist) using, ///
	controls(x3e_varlist) [options]

note: must be called with preserve and restore
*/
program _dsra_parse_syntax, rclass
	version 16.0

	_on_colon_parse `0'
	local before `s(before)'
	local after `s(after)'

	local 0 `before'
	syntax , touse1(passthru) touse2(passthru)

	local 0 `after'
	syntax anything(equalok name=eq) 	///
		using 				///
		[if] [in],			///
		controls(passthru)		///
		[ lassoopts(passthru)		///
		true_gamma(passthru)		///
		allgamma			///
		noadjust]

	/*-- append data ---*/
	AppendData `using' `if' `in', `touse1' `touse2'
	return add

	/*-- parse equation --*/
	ParseEq `eq', `controls' `touse1' `touse2'
	return add

	/*--- parse lassoopts ---*/
	ParseLassoOpts, `lassoopts'
	return add

	/*--- parse true_gamma ---*/
	ParseTrueGamma, `true_gamma' `allgamma'
	return add

	/*--- parse noadjust -----*/
	ParseNoadjust, `adjust'
	return add
end

					//----------------------------//
					//  Append data
					//----------------------------//
program AppendData, rclass
	syntax using [if] [in],	///
		touse1(string)	///
		touse2(string)
	
	/*----- append data ----*/
	tempvar gr
	append `using', gen(`gr')
	
	/*-- marksample ---*/
	marksample touse

	/*----- touse1 -------*/
	gen byte `touse1' = 0
	qui replace `touse1' = `touse' if `gr' == 0

	/*----- touse2 -------*/
	gen byte `touse2' = 0
	qui replace `touse2' = `touse' if `gr' == 1

	ret local touse1 `touse1'
	ret local touse2 `touse2'
end

					//----------------------------//
					//  parse equation
					//----------------------------//
program ParseEq, rclass
	syntax anything(name=eq equalok)	///
		, controls(string)		///
		touse1(string)			///
		touse2(string)

	_iv_parse `eq'
	local depvar `s(lhs)'
	local x1_vars `s(exog)'
	local x2_vars `s(endog)'
	local x3i_vars `s(inst)'
	
	_fv_check_depvar `depvar'
	ret local depvar `depvar'

	ParseVars `x1_vars' if `touse1'
	ret local x1_vars `s(vars)'

	ParseVars `x2_vars' if `touse2'
	local tmp `s(vars)'
	if (`:list sizeof tmp' > 1) {
		di as err "x2 must be a scalar"
		exit 198
	}
	ret local x2_vars `s(vars)'


	ParseVars `x3i_vars' if `touse1'
	ret local x3i_vars `s(vars)'

	ParseVars `controls' if `touse2', normcoll
	ret local x3e_vars `s(vars)'
end
					//----------------------------//
					//  Parse vars
					//----------------------------//
program ParseVars, sclass
	syntax varlist(numeric fv) if	///
		[, normcoll]

	if (`"`rmcoll'"' == "normcoll") {
		fvexpand `varlist' `if'
	}
	else {
		_rmcoll `varlist' `if', expand
	}
	local vars `r(varlist)'

	sret local vars `vars'
end
					//----------------------------//
					//  parse lassoopts
					//----------------------------//
program ParseLassoOpts, rclass
	syntax [, lassoopts(string) ] 

	if (`"`lassopts'"' == "") {
		local lassoopts sel(plugin)
	}

	ret local lassoopts `lassoopts'
end

					//----------------------------//
					//  parse lassoopts
					//----------------------------//
program	ParseTrueGamma, rclass
	syntax [, true_gamma(string) 	///
		allgamma ]

	if (`"`true_gamma'"' == "") {
		exit
		// NotReached
	}

	confirm matrix `true_gamma'
	ret local true_gamma `true_gamma'
	ret local allgamma `allgamma'
end
					//----------------------------//
					//  parse noadjust
					//----------------------------//
program ParseNoadjust, rclass
	syntax [, noadjust]

	if (`"`adjust'"' == "") {
		local adjust adjust
	}

	ret local adjust `adjust'
end
