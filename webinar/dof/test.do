cscript

set seed 12345671

adopath ++ ../../src

use psid2013_analysis, clear
global logs ../logs

/*-----------------------------------------------------------------------------
	define the control variables in the selection equation	
-----------------------------------------------------------------------------*/

sjlog using $logs/log1, replace
local vars_sel exper exper2 educ_level	childcare_expen_2012 		///
	i.if_kidsle15 num_kids wage_husband exp_appl i.wtr_enrolled	///
	i.wtr_grad_hs i.wtr_attend_college i.wtr_cert_educ 		///
	i.wtr_educ_usa	i.father_educ_usa i.mother_educ_usa		/// 
	i.rural_urban i.own_vehicle i.current_state
sjlog close, replace
	

sjlog using $logs/log2, replace
dsheckman lnwage educ_level exper, selection(inlf = `vars_sel') 	
sjlog close, replace

sjlog using $logs/log3, replace
dsheckman lnwage educ_level exper, selection(inlf = `vars_sel')	///
	selvars(num_kids educ_level exper)
sjlog close, replace

