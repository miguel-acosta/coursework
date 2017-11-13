cap program drop textable
program define textable

local fname `1'

local outstr = ""

ds
local vl = r(varlist)
local N = _N

file open hi using "`fname'.tex", write text replace


foreach ii of numlist 1/`N' {
    local outstr = ""
    foreach vv of varlist `vl' {
        local vstr = string(`vv') in `ii'
        local outstr = "`outstr' `vstr'&"
    }
    local outstr = substr("`outstr'", 1,strlen("`outstr'")-1)
    local outstr = "`outstr'\\"
    file write hi "`outstr'"  _n
}
file close hi

end
