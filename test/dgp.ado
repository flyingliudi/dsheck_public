*! version  1.0.0
/*
	DGP for two sample
*/
program dgp
	syntax , nobs1(string)	///
		nobs2(string)	///
		kx1(string)	///
		kx2(string)	///
		kx3i(string)	///
		kx3e(string)	///
		s1(string)	///
		s2(string)	///
		s0(string)

	/*---- set obs ----*/
	clear
	local nobs = `nobs1' + `nobs2'
	set obs `nobs'

	/*---- set x3e -----*/
	mata: mk_toeplitz(`kx3e', `nobs', "x3e")
	forvalues i=1/`kx3e' {
		Normalize x3e`i'
	}

	/*---- set x3i ----*/
	forvalues i=1/`kx3i' {
		local k = 10*(`i'-1) + 1
		gen double x3i`i' = x3e`k'
		forvalues m=1/5 {
			local k = `k' + `m'
			local b = runiform(0.1, 3)
			replace x3i`i' = `b'*x3e`k' + x3i`i'
		}
		
		Normalize x3i`i'

		replace x3i`i' = x3i`i' + rnormal()
	}

	/*---- set x1 ----*/
	forvalues i=1/`kx1' {
		local k = 15*(`i'-1) + 1
		gen double x1`i' = x3e`k'
		forvalues m=1/5 {
			local k = `k' + `m'
			local b = runiform(0.1, 3)
			replace x1`i' = `b'*x3e`k' + x1`i'
		}
		Normalize x1`i'

		replace x1`i' = x1`i' + rnormal()
	}

	/*----- set x2 -----*/
	forvalues i=1/`kx2' {
		local k = 11*(`i'-1) + 1
		gen double x2`i' = x3e`k'
		forvalues m=1/5 {
			local k = `k' + `m'
			local b = runiform(0.1, 1)
			replace x2`i' = `b'*x3e`k' + x2`i'
		}

		local b = runiform(0.1, 3)
		replace x2`i' = x2`i' + x3i1*`b' 
		Normalize x2`i'
		replace x2`i' = x2`i' + rnormal(0, 1)
	}
		
	/*------- set y -------*/
	gen double y =  0
	forvalues i=1/`kx1' {
		replace y  = y + x1`i'
	}

	forvalues i=1/`kx2' {
		replace y  = y + x2`i'
	}

	forvalues i=1/`kx3i' {
		replace y  = y + x3i`i'
	}

	replace y = y + runiform(1,1.5)*(rchi2(1) + rnormal(0, 1))

	/*------ split the original data into two samples */
	gen double os = rnormal()
	sort os
	gen gr = 2
	replace gr = 1 in 1/`nobs1'

	/*------- first sample -------*/
	preserve
	drop x2*
	keep if gr == 1
	save `s1', replace
	restore

	/*------- second sample -----*/
	preserve
	keep if gr == 2
	save `s2', replace
	restore

	/*------- original sample ----*/
	save `s0', replace

	use `s1', clear
end

					//----------------------------//
					//  normalize
					//----------------------------//
program Normalize
	syntax varname

	sum `varlist'
	replace `varlist' = `varlist' - r(mean)
	replace `varlist' = `varlist'/r(sd)
end

mata:

mata set matastrict on 
void mk_toeplitz(		///
	real scalar	p,	///
	real scalar	n,	///
	string scalar	xs)	
{
	real vector	d, r, idx 
	real matrix	V, L, X
	real scalar	df, toeplitz
	
	df = 15
	toeplitz = 1.3

	d = (1..p)
	r = d:^(-toeplitz)	
	V = Toeplitz(r', r)
	L = cholesky(V)
	X    = (1/(sqrt(2*df)))*(rchi2(n, p, df) :- df)
	X    = X*L'			

	idx = st_addvar("double", xs:+strofreal((1..p)))
	st_store(., idx, X)
}

end
