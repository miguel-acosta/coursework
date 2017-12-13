include("../jllib/TSfun.jl")
using Distributions

function olsDF(y,P; trend = true, chop = 0, cons = true)
    ## Allowing us to trim some observations for IC tests
    if chop > 0
        y = y[(chop+1):end,:]
    end

    ## Construct differences and lags 
    TT   = length(y)
    Δy   = y[2:TT] - y[1:TT-1]
    LΔy  = lagmatrix(Δy, P)
    TTreg = size(LΔy)[1]
    yreg = LΔy[:,1]
    
    ## Construct X matrix depending on whether there is a constant/trend
    if trend
        xreg = [ones(TTreg) 1:TTreg LΔy[:,2:end]]
    elseif !trend & cons
        xreg = [ones(TTreg)         Δy[:,2:end]]
    elseif !trend & !cons
        reg = Δy[:,2:end]
    else
        print("WHAT ARE YOU DOING?! A TREND AND NO CONSTANT?!")
    end

    ## Run the regressions 
    β  = (xreg.' * xreg)\ xreg.' * yreg
    Σ  = inv((xreg.' * xreg)) *  (sum((yreg - xreg * β).^2) / (TTreg^2))
    se = sqrt.(diag(Σ))
    return(β, se, Σ)
end
    

function BIC(y ; trend = true, pmax = 4)
    TT= length(y)
    K = trend ? 2 : 1
    Σ(m)   = olsDF(y,m;chop=(pmax-m))[3]
    BIC(m) = log(det(Σ(m))) + (log(TT)/TT) * (m * K^2 + K)
    BICall = [BIC(m) for m in 1:pmax]
    BICmin = findmin(BICall)[2]
    return(BICmin)
end


function dickeyfuller(y; trend = true, cons = true)
    TT = length(y)
    p  = BIC(y)
    β, se, Σ = olsDF(y,p)
    df = trend ? β[3] : (cons ? β[2] : β[1])
    return(df)
end

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

    dfgls = dickeyfuller(yd; trend = false, cons = false)
    return(dfgls)
end

function buildY(TT)
    e = rand(Normal(0,1),TT)
    y = zeros(TT)
    for tt in 2:TT
        y[tt] = y[tt-1] + e[tt]
    end
    return(y)
end

ysamp = buildY(200)
