program dgp_heckman, rclass
	syntax [, px(integer 2) 	///
		pz(integer 100)		///
		nobs(integer 1000)	///
		gamma(real 0.2)		///
		beta(real 1) ]
	
	/*----- set obs -------*/
	clear
	qui set obs `nobs'

	/*---- set x and z-----*/
	foreach var in x z {
		mata: mk_toeplitz(`p`var'', `nobs', "`var'")
		forvalues i=1/`p`var'' {
			qui Normalize `var'`i'
		}
	}

	/*----- set v2 and u1 ----*/
	qui gen double v2 = rnormal()
	qui gen double u1 = `gamma'*v2 + rnormal()

	/*----- set y1 and y2 ----*/
	replace x1 = x1 + z3 +z10 + z11 +z12 + z15
	replace x2 = x2 + 0.5*(z3 +z10 + z11 +z12 + 2*z15)
	qui gen y1 = `beta' + `beta'*x1 + `beta'*x2 + u1

	local true_zvars z1 z2 z3 z5 z10 z11 z12 z15 z20

	foreach var of local true_zvars {
		replace `var' = `var' + (x1 - x2) + rchi2(1)
		qui Normalize `var'
	}

	local bmin = 0.3*sqrt(4*log(`pz')/`nobs')
	qui gen y2 = (-1.5 + z1 - z2 + `bmin'*z3 		///
		+ 1*z5 - `bmin'*z10  - `bmin'*z11 + `bmin'*z12  ///
		-`bmin'*z15 + z20 + v2 > 0)

	qui replace y1 = . if y2 == 0

	ret local true_zvars `true_zvars'
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
