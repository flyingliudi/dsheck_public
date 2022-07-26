*! version 1.0.0  21apr2021
program poprobit_gmm_probit
	syntax varlist if , 		///
		at(name) 		///
		depvar(string)		///
		tt_vars(string)		///
		xbvar(string)		///
		[derivatives(varlist) ]
	
	tempvar xb
	qui matrix score double `xb' = `at' `if' 

	tempvar myxb
	qui gen double `myxb' = `xb' + `xbvar' `if'

	tempvar pr
	qui gen double `pr' = normal(`myxb')

	qui replace `varlist' = `depvar' - `pr' `if'

	if (`"`derivatives'"' == "") {
		exit
		// NotReached
	}

	qui replace `derivatives' = -normalden(`myxb') `if'
end
