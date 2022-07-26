*! version 1.0.0  06apr2020
					//----------------------------//
					// Display
					//----------------------------//
program _dsra_display
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

	local k_controls_sel as txt _col(`col1') 	///
		`"Number of selected controls"' 	///
		_col(`col2') "=" _col(69) as res %10.0fc e(k_controls_sel)
	
	local n_reps as txt _col(`col1') 	///
		`"Number of replications"' 	///
		_col(`col2') "=" _col(69) as res %10.0fc e(N_reps)


	di `title' `n_obs1'
	di `n_obs2'
	di `k_controls'
	di `k_controls_sel'
	if (e(N_reps) != .) {
		di `n_reps'
	}
	di 
end
