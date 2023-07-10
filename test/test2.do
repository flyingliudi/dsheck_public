cscript

mata:
void mk_toeplitz(		///
	real scalar	p,	///
	real scalar	n)
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

	C = correlation(X)
	C
}
end

mata :
mk_toeplitz(p =10, n=1000)
end
