*! version  1.0.0
program _display_heckman
	
	local title as txt "`e(title)'"

	
	local col = 39
						//  nobs
	local nobs _col(`col') as txt "Number of obs" _col(67) "="	///
		_col(69) as res %10.0fc e(N)

	local col2 = `col' + 7
	local n_sel _col(`col2') as txt "Selected" _col(67) "="	///
		_col(69) as res %10.0fc e(N_sel)

	local n_nonsel _col(`col2') as txt "Nonselected" _col(67) "="	///
		_col(69) as res %10.0fc e(N_nonsel)
	
	if (e(p)!=.) {
		local p _col(`col') as txt "Number of variables" ///
			_col(67) "=" _col(69) as res %10.0fc e(p)
	}

	if (e(k)!=.) {
		local k _col(`col') as txt 			///
			"Number of selected controls" 		///
			_col(67) "=" _col(69) as res %10.0fc e(k)
	}
	
	if (e(k1)!=.) {
		local k1 _col(`col') as txt 			///
			"Number of main variables" 		///
			_col(67) "=" _col(69) as res %10.0fc e(k1)
	}

	di
	di `title' `nobs'
	di `n_sel'
	di `n_nonsel'
	if (e(k1)!=.) {
		di `p'
		di `k'
		di `k1'
	}
	di
	_coef_table
	Footnote
end

						//------------------------//
						//  Footnote
						//------------------------//
program Footnote

	if (`"`e(cmd)'"' != "dsheckman" & `"`e(cmd)'"' != "naiveheckman") {
		exit
		// NotReached
	}

	local k1 = e(k1)
	local k = e(k)
	local p = e(p)
	di as txt "{p 0 6 2}Note: in the main equation, there are "	///
		"{res:`k1'} variables; in the selection equation, "	///
		"{res:`k'} among {res:`p'} variables are used to "	///
		"predict inverse mills ratio.{p_end}"
end
