*!version 1.0.0  17jun2020
program adjust_heckman_V
	syntax, touse_main(string)	///
		zr(string)		///
		u1(string)		///
		bm(string)		///
		xvars(string)		///
		zvars(string)		///
		lambda(string)		///
		vp(string)		///
		vs(string)		///
		[controls_sel(string)]

	mata: _heckman_adjust_V( 	///
		`"`touse_main'"',	///
		`"`zr'"', 		///
		`"`u1'"',		///
		`bm',			///
		`"`xvars'"',		///
		`"`zvars'"',		///
		`"`lambda'"',		///
		`"`controls_sel'"',	///
		`"`vp'"',		///
		`"`vs'"')

end

mata:
mata set matastrict on

void _heckman_adjust_V(			///
	string scalar	_touse_main,	///
	string scalar	_zr,		///
	string scalar	_u1,		///
	real scalar	_bm,		///
	string scalar	_xvars,		///
	string scalar	_zvars,		///
	string scalar	_lambda,	///
	string scalar	_controls_sel,	///
	string scalar	_Vp,		///
	string scalar	_Vs)
{
	real colvector	zr, m, delta, u1, R
	real scalar	sig2, N, rho
	real matrix	W, invWpW, WpRW, Q, V

	zr = st_data(., _zr, _touse_main)
	m = normalden(zr):/normal(zr)
	delta = m:*(m+zr)
	u1 = st_data(., _u1, _touse_main)

	N = length(u1)
	sig2 = (cross(u1, u1) + _bm^2*sum(delta))/N

	rho = _bm/sqrt(sig2)

	W = st_data(., _xvars + " " + _lambda+ " " + _controls_sel, _touse_main)
	W = W, J(N, 1, 1)
	invWpW = invsym(cross(W, W))
	R = 1:-rho^2*delta
	WpRW = cross(W, R, W)

	Z = st_data(., _zvars, _touse_main)
	Z = Z, J(N, 1, 1)
	WpDZ = cross(W, delta, Z)
	Vp = st_matrix(_Vp)
	Q = rho^2*WpDZ*Vp*(WpDZ')

	V = sig2*invWpW*(WpRW + Q)*invWpW
	st_matrix(_Vs, V)
}

end
