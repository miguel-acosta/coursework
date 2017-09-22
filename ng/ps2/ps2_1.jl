using PyPlot
## A few functions will be helpful first.
@everywhere function OLS(y,x,; se = "OLS", cons = true, qHAC = 10)
    TT    = length(y)
    if cons
        x = [repmat([1],TT,1) x]
    end
    β = (x.'*x)\x.' *y

    if se == "OLS"
        return β, OLSse(y,x,β)
    elseif se == "White"
        return β, WhiteSE(y,x,β)
    elseif se == "HAC"
        return β, HACse(y,x,β,qHAC)
    else
        return("Invalid std. error name")
    end
end

@everywhere function OLSse(y,x,β)
    TT    = length(y)
    error = y - x * β
    s2    = sum(error.^2)/(TT-2)
    V     = inv(x.'*x) * s2
    return(sqrt.(diag(V)))
end

@everywhere function WhiteSE(y,x,β)
    xtxi = inv(x.'*x)
    u    = y - x * β
    V = xtxi * (x.' * diagm(u.^2) * x) * xtxi
    return(sqrt.(diag(V)))
end

## This is equation 10.5.21 of Hamilton
@everywhere function HACse(y,x,β,q)
    k    = size(x,2)
    TT    = length(y)    
    if q > TT
        error("your q is too big... q > T")
    end
    xtxi = inv(x.'*x)
    u    = y - x * β

    ## Build the terms 1 by 1
    M1   = zeros(k,k)
    for tt = 1:TT
        M1 = M1 + u[tt]^2 * x[tt,:] * (x[tt,:].')
    end
    M2   = zeros(k,k)
    for vv = 1:q
        M2inner = zeros(k,k)
        for tt = (vv+1):TT
            M2inner = M2inner + x[tt,:] * u[tt] * u[tt-vv] * x[tt-vv,:].'
                              + x[tt-vv,:] * u[tt-vv] * u[tt] * x[tt,:].'
        end
        M2 = M2 + (1-vv/(q+1)) * M2inner
    end
    V = xtxi * (M1 + M2) * xtxi
    
    return(sqrt.(diag(V)))
end

srand(6413)

beta0 = 0
beta1 = 0
kappa = 2
alpha = 0.35
delta = 0.6
TT    = 500
nsim  = 5000

v = randn(TT,nsim)
se_OLS   = SharedArray{Float64}(nsim)
se_White = SharedArray{Float64}(nsim)
se_HAC   = SharedArray{Float64}(nsim)

beta_hat   = SharedArray{Float64}(nsim)

@parallel (+)  for nn = 1:nsim
    y = zeros(TT,1)
    h = zeros(TT,1)
    u = zeros(TT,1)
    h[1] = kappa/(1-alpha-delta)
    u[1] = sqrt(h[1]) * v[1,nn]

    for tt in 2:TT
        h[tt] = kappa + alpha * u[tt-1]^2 + delta * h[tt-1]
        u[tt] = sqrt(h[tt]) * v[tt,nn]
        y[tt] = beta0 + beta1 * y[tt-1] + u[tt]
    end
    ols   = OLS(y[2:end], y[1:end-1], se = "OLS")
    white = OLS(y[2:end], y[1:end-1], se = "White")
    hac   = OLS(y[2:end], y[1:end-1], se = "HAC", qHAC = 5)

    se_OLS[nn]   =  ols[2][2]
    se_White[nn] =  white[2][2]
    se_HAC[nn]   =  hac[2][2]

    beta_hat[nn] =  ols[1][2]
end

reject_ols   = abs.(beta_hat./se_OLS)   .> 1.96
reject_White = abs.(beta_hat./se_White) .> 1.96
reject_hac   = abs.(beta_hat./se_HAC)   .> 1.96

size_ols   = sum(reject_ols)/nsim
size_white = sum(reject_White)/nsim
size_hac   = sum(reject_hac)/nsim

include("textable.jl")
textable(["OLS", "White", "HAC, \$q=5\$"], [size_ols,size_white,size_hac], precision = "%2.3f", fname = "q1")


function hist(series, name, TITLE)
    ioff()
    f = figure(figsize = (4.5, 4.5))
    plt[:hist](series,100)
    title(TITLE)
    savefig(string(name,".pdf"))
    close(f)
end
hist(beta_hat, "beta", "Estimated Coefficients Across Simulations")
hist(se_OLS, "seOLS", "Estimated OLS Std. Errors Across Simulations")
hist(se_White, "seWhite", "Estimated White Std. Errors Across Simulations")
hist(se_HAC, "seHAC", "Estimated HAC Std. Errors Across Simulations")
            
