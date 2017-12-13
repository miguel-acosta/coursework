include("../jllib/TSfun.jl")
using Distributions, PyPlot, KernelDensity
srand(6413)
function ols(y,X; β0 = 0, cons = true, trend = false)
    TT = length(y)
    if cons & trend
        X  = [ones(1:TT) 1:TT copy(X)]
    elseif cons
        X  = [ones(1:TT) copy(X)]
    end
    β  = (X.' * X)\ X.' * y
    Σ  = inv((X.' * X)) *  (sum((y - X * β).^2) / length(y))
    se = length(size(X)) == 1 ? sqrt(Σ) : sqrt.(diag(Σ))
    t = (β-β0)./se
    return(β, se, t, Σ)
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

    return(yd)
end

function buildY(TT)
    e = rand(Normal(0,1),TT)
    y = zeros(TT)
    for tt in 2:TT
        y[tt] = y[tt-1] + e[tt]
    end
    return(y)
end

function sim(TT)
    y    = buildY(TT)
    yd   = GLSdetrend(y)
    Ly   = lagmatrix(y,1)
    Lyd  = lagmatrix(yd,1)
    Xols = [ones(TT-1) 1:(TT-1) Ly[:,2]]
    
    αols  = ols(Ly[:,1],  Ly[:,2],  cons = true,  trend = true,  β0 = 1)[1][3]
    αgls  = ols(Lyd[:,1], Lyd[:,2], cons = false, trend = false, β0 = 1)[1][1]
    
    t1ols =  ols(Ly[:,1],  Ly[:,2],  cons = true,  trend = true,  β0 = 1)[3][3]
    t1gls =  ols(Lyd[:,1], Lyd[:,2], cons = false, trend = false, β0 = 1)[3][1]

    t9ols =  ols(Ly[:,1],  Ly[:,2],  cons = true,  trend = true,  β0 = 0.9)[3][3]
    t9gls =  ols(Lyd[:,1], Lyd[:,2], cons = false, trend = false, β0 = 0.9)[3][1]

    t98ols = ols(Ly[:,1],  Ly[:,2],  cons = true,  trend = true,  β0 = 0.98)[3][3]
    t98gls = ols(Lyd[:,1], Lyd[:,2], cons = false, trend = false, β0 = 0.98)[3][1]

    return(αols,αgls,t1ols,t1gls,t9ols,t9gls,t98ols,t98gls)    
end

function sims(nsim)
    α    = zeros(nsim,2)
    t100 = zeros(nsim,2)
    t098 = zeros(nsim,2)
    t090 = zeros(nsim,2)
    for ss in 1:nsim 
        αols,αgls,t1ols,t1gls,t9ols,t9gls,t98ols,t98gls = sim(200)
        α[ss,:]    = [αols   αgls]
        t100[ss,:] = [t1ols  t1gls]
        t098[ss,:] = [t98ols t98gls]
        t090[ss,:] = [t9ols t9gls]
    end    
    return(α, t100, t098, t090)
end

function summary(α,t100,t098,t090)
    ## Plot the distribution of α
    KDEαOLS = kde(α[:,1])
    KDEαGLS = kde(α[:,2])
    plot(KDEαOLS.x, KDEαOLS.density/sum(KDEαOLS.density), label = "OLS")
    plot(KDEαGLS.x, KDEαGLS.density/sum(KDEαGLS.density), label = "GLS")
    xlabel("α")
    legend()
    title(string("Distribution of α"))
    savefig("output/alph.pdf")
    close()

    ## Plot the distribution of t-statistics for the various
    ## tests: OLS
    tsims = [t100, t098, t090]
    tnames = ["1", "0.98", "0.90"]
    for tt in 1:length(tsims)
        KDEtOLS = kde(tsims[tt][:,1])
        plot(KDEtOLS.x, KDEtOLS.density/sum(KDEtOLS.density),
             label = string("OLS, α = ", tnames[tt]))
    end
    xlabel("t")
    legend()
    title(string("Distribution of t-stats: OLS"))
    savefig("output/tOLS.pdf")
    close()

    ## Plot the distribution of t-statistics for the various
    ## tests: OLS
    for tt in 1:length(tsims)
        KDEtGLS = kde(tsims[tt][:,2])
        plot(KDEtGLS.x, KDEtGLS.density/sum(KDEtGLS.density),
             label = string("GLS, α = ", tnames[tt]))        
    end
    xlabel("t")
    legend()
    title(string("Distribution of t-stats: GLS"))
    savefig("output/tGLS.pdf")
    close()

    ## Plot the distribution of t-statistics for the various
    ## tests: OLS and GLS
    for tt in 1:length(tsims)
        KDEtOLS = kde(tsims[tt][:,1])
        KDEtGLS = kde(tsims[tt][:,2])
        plot(KDEtOLS.x, KDEtOLS.density/sum(KDEtOLS.density),
             label = string("OLS, α = ", tnames[tt]))
        plot(KDEtGLS.x, KDEtGLS.density/sum(KDEtGLS.density),
             label = string("GLS, α = ", tnames[tt]))
        xlabel("t")
        legend()
        title(string("Distribution of t-stats"))
        savefig(string("output/t", tnames[tt], ".pdf"))
        close()            
    end

end

α, t100, t098, t090 = sims(20000)
summary(α, t100, t098, t090)
