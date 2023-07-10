program dgp2_heckman, rclass
	syntax [, px(integer 2) 	///
		pz(integer 100)		///
		nobs(integer 1000)	///
		gamma(real 0.2)		///
		bbig(real 0.4)		///
		bsmall(real 0.2)	///
		beta(real 1) ]
	
					//  set obs
	clear
	qui set obs `nobs'

					//  set x and z
	local pall = `px' + `pz'
	mata: mk_toeplitz(`pall', `nobs', "var")

	forvalues i = 1/`px' {
		rename var`i'0 x`i'
	}

	rename var# z#, renumber
	
					// v2 and u1
	qui gen double v2 = rnormal()
	qui gen double u1 = `gamma'*v2 + rnormal(0, 1)

					// y1
	qui gen y1 = `beta' + `beta'*x1 + `beta'*x2 + u1

					// y2
	local true_zvars x1 x2 z1 z2 z3 z5 z10 z11 z12 z15 z20
	qui gen y2 = (`bbig' + `bbig'*x1 + `bbig'*z1 + `bbig'* z2 	///
		+ `bsmall'*z3 + `bbig'*z5 + `bsmall'*z10   		///
		+ `bsmall'*z11 + `bsmall'*z12  ///
		-`bsmall'*z15 + `bbig'*z20 + `bsmall'*x2 + v2 > 0)

	qui replace y1 = . if y2 == 0

	ret local true_zvars `true_zvars'

	order *, sequential
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
	
	df = 20
	toeplitz = 1.1

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
