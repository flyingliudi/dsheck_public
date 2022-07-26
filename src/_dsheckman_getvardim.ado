						//------------------------//
						// get variable dimension
						//------------------------//
program _dsheckman_getvardim, rclass
	syntax, xvars(string)		///
		zvars_sel(string)	///
		zvars_nosel(string)
	
	fvexpand `xvars'
	local tmp `r(varlist)'
	local k1 : list sizeof tmp

	fvexpand `zvars_sel'
	local tmp `r(varlist)'
	local k : list sizeof tmp

	fvexpand `zvars_nosel'
	local tmp `r(varlist)'
	local p : list sizeof tmp
	local p = `p' + `k'
	

	return scalar p = `p'
	return scalar k1 = `k1'
	return scalar k = `k'
end

