cscript

adopath ++ ../src

local seed 12345671
set seed `seed'

local nobs = 1000
local px = 2
local pz = 500
local gamma = 0.7
local beta = 1

dgp_heckman, nobs(`nobs') px(`px') pz(`pz') gamma(`gamma') beta(`beta')

dsheckman y1 x*, sel(y2 = x* z*) selvars(x* z1-z20)
