cscript

adopath ++ ../src


local nobs = 1000
local px = 2
local pz = 50
local gamma = 0.5
local beta = 1
local bsmall = 1
local bbig = 1.5
local nreps = 35
local vars x1 x2 z1 
local true_b  `bbig' `bsmall' `bbig' 

					// set controls
qui dgp2_heckman, nobs(`nobs') px(`px') pz(`pz') gamma(`gamma') ///
	beta(`beta') bsmall(`bsmall') bbig(`bbig')


dsheckman y1 x*, sel(y2 = x* z*)	

