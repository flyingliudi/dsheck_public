*! version 1.0.0  14apr2020
/*----------------------------------------------------------------------
	bootstrap -dsra-
----------------------------------------------------------------------*/
program dsra_bootstrap
	version 16.0

	syntax  anything(name=eq) 	///
		using/			///
		[if] [in]		///
		[, reps(integer 100)	///
		* ]

	//*------ preserve START ----------*/
	preserve

	/*------- 1: dmse using original data -----*/
	tempname est_orig
	di as txt "dsra using original data"
	qui dsra `eq' `if' `in' using `using', `options'
	est store `est_orig'

	/*------- 2: get original data -------*/
	tempfile orig1 orig2	
	qui save `orig1', replace
	
	use `using', clear
	qui save `orig2', replace

	/*------- 3: bootstrap dsra ---------*/
	tempfile b2
	tempname B 
	di as txt _n "start bootstrap {cmd:dsra} ..."
	forvalues i = 1/`reps' {

		use `orig2', clear
		bsample
		qui save `b2', replace

		use `orig1', clear
		bsample

		qui dsra `eq' `if' `in' using `b2', `options'
		mat `B' = (nullmat(`B') \ e(b))

		di_dot, i(`i') 
	}
	di as txt "bootstrap {cmd:dsra} finished"

	/*------ preserve END ----------*/
	restore

	/*------ post result ----------*/
	PostResult, bmat(`B') est_orig(`est_orig') reps(`reps')

	/*------ Display ------------*/
	_dsra_display
end

					//----------------------------//
					//  display dot
					//----------------------------//
program di_dot
	syntax, i(string)
	
	if (mod(`i', 50) == 0) {
		di as txt "   `i'"
	}
	else {
		di as txt _c "."
	}
end

					//----------------------------//
					//  Post result
					//----------------------------//
program PostResult, eclass
	syntax, est_orig(string)	///
		bmat(string)		///
		reps(string)
	
	qui est restore `est_orig'
	
	tempname Vs
	mata: get_V(`"`bmat'"', `"`Vs'"')

	eret repost V = `Vs'
	eret scalar N_reps = `reps'
	eret local vce bootstrap
	eret local vcetype Bootstrap 
	eret local cmd dsra_bootstrap
end

mata :
mata set matastrict on
					//----------------------------//
					//  get V
					//----------------------------//
void get_V(			///
	string scalar	bs,	///
	string scalar	Vs)
{
	real matrix	B, V
	real scalar	k
	real rowvector	b

	B = st_matrix(bs)
	k = rows(B)
	b = mean(B)

	V = B:-b
	V = cross(V, V)/(k-1)

	st_matrix(Vs, V)
	st_matrixrowstripe(Vs, st_matrixcolstripe("e(b)"))
	st_matrixcolstripe(Vs, st_matrixcolstripe("e(b)"))
}
end

