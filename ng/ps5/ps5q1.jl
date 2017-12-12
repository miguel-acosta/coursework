function GLSdetrend(y; trend = true)
    TT     = length(y)
    ytilde = copy(y)
    xtilde = ones(TT)    

    x = trend ? [ones(TT) 1:TT] : ones(TT)
    α = trend ? 1-13.5/TT       : 1-7/TT

    for tt in 2:TT
        xtilde[tt,:] = x[tt,:]  - α * x[tt-1,:]
        ytilde[tt] = y[tt]  - α * y[tt-1]
    end
    δ = (xtilde.' * xtilde) \ xtilde.' * ytilde
    yd = y - x.' * δ
end
