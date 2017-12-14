include("../jllib/TSfun.jl")
include("summaryPlots.jl")
using Distributions, PyPlot, KernelDensity

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function ols(y,X; β0 = 0, cons = true, trend = false, trendbreak = false)
    TT = length(y)
    if cons & trend & trendbreak > false
        tbreak = [ones(trendbreak); zeros(TT-trendbreak)]
        X  = [ones(TT) collect(1:TT) copy(X) tbreak.*collect(1:TT) tbreak]
    elseif cons & trend
        X  = [ones(TT) collect(1:TT) copy(X)]
    elseif cons
        X  = [ones(TT) copy(X)]
    end
    β  = (X.' * X)\ X.' * y
    Σ  = inv((X.' * X)) *  (sum((y - X * β).^2) / length(y))
    se = length(size(X)) == 1 ? sqrt(Σ) : sqrt.(diag(Σ))
    t = (β-β0)./se
    return(β, se, t, Σ)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function GLSdetrend(y; trend = true)
    TT     = length(y)
    ytilde = copy(y)
    xtilde = trend ? ones(TT,2)   : ones(TT,1)    

    x = trend ? [ones(TT) 1:TT] : ones(TT)
    α = trend ? 1-13.5/TT       : 1-7/TT

    for tt in 2:TT
        xtilde[tt,:] = x[tt,:]  - α * x[tt-1,:]
        ytilde[tt] = y[tt]  - α * y[tt-1]
    end
    δ = (xtilde.' * xtilde) \ xtilde.' * ytilde
    yd = y - x * δ

    return(yd)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function buildY(TT)
    e = rand(Normal(0,1),TT)
    y = zeros(TT)
    for tt in 2:TT
        y[tt] = y[tt-1] + e[tt]
    end
    return(y)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function buildYbreak(TT)
    e = rand(Normal(0,1),TT)
    y = zeros(TT)
    for tt in 2:TT
        y[tt] = 0.92 * y[tt-1] + (tt<=75 ? 0 : -0.03) * tt + e[tt]
    end
    return(y)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function sim(TT, tbreak, realTbreak)
    if realTbreak
        y    = buildYbreak(TT)
    else
        y    = buildY(TT)
    end
    yd   = GLSdetrend(y)
    Ly   = lagmatrix(y,1)
    Lyd  = lagmatrix(yd,1)
    
    αols  = ols(Ly[:,1],  Ly[:,2],  cons = true,  trend = true, trendbreak = tbreak,
                β0 = 1)[1][3]
    αgls  = ols(Lyd[:,1], Lyd[:,2], cons = false, trend = false, trendbreak = false,
                β0 = 1)[1][1]
    
    t1ols =  ols(Ly[:,1],  Ly[:,2],  cons = true,  trend = true, trendbreak = tbreak,
                 β0 = 1)[3][3]
    t1gls =  ols(Lyd[:,1], Lyd[:,2], cons = false, trend = false, trendbreak = false,
                 β0 = 1)[3][1]

    t9ols =  ols(Ly[:,1],  Ly[:,2],  cons = true,  trend = true, trendbreak = tbreak,
                 β0 = 0.9)[3][3]
    t9gls =  ols(Lyd[:,1], Lyd[:,2], cons = false, trend = false, trendbreak = false,
                 β0 = 0.9)[3][1]

    t98ols = ols(Ly[:,1],  Ly[:,2],  cons = true,  trend = true, trendbreak = tbreak,
                 β0 = 0.98)[3][3]
    t98gls = ols(Lyd[:,1], Lyd[:,2], cons = false, trend = false, trendbreak = false,
                 β0 = 0.98)[3][1]

    return(αols,αgls,t1ols,t1gls,t9ols,t9gls,t98ols,t98gls)    
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function sims(nsim; tbreak = false, realTbreak = false)
    α    = zeros(nsim,2)
    t100 = zeros(nsim,2)
    t098 = zeros(nsim,2)
    t090 = zeros(nsim,2)
    for ss in 1:nsim 
        αols,αgls,t1ols,t1gls,t9ols,t9gls,t98ols,t98gls = sim(200,tbreak,realTbreak)
        α[ss,:]    = [αols   αgls  ]
        t100[ss,:] = [t1ols  t1gls ]
        t098[ss,:] = [t98ols t98gls]
        t090[ss,:] = [t9ols  t9gls ]
    end    
    return(α, t100, t098, t090)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
NSIM = 20000
## Part a
srand(6413)
α, t100, t098, t090 = sims(NSIM)
summaryPlots(α, t100, t098, t090)

## Part b
critval(t100, t098, t090, "a"; GLS = true)

## Part c
srand(6413)
α_c, t100_c, t098_c, t090_c = sims(NSIM, tbreak = 75)
critval(t100_c, t098_c, t090_c, "c"; GLS = false)

## Part d
srand(6413)
α_d, t100_d, t098_d, t090_d = sims(NSIM, tbreak = 75, realTbreak = true)
critval(t100_d, t098_d, t090_d, "d"; GLS = false)
