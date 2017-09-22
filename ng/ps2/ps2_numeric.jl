using PyPlot, Optim, SymPy
include("textable.jl")

function simulateARMA(α,θ,σ,TT)
    e = zeros(TT)
    y = zeros(TT)
    randns = randn(TT-1,1)
    e[2:TT] = randns * σ
    for tt in 2:TT
        y[tt] = α * y[tt-1] + e[tt] + θ * e[tt-1]
    end
    return y,e
end

function plotts(y, name, snames)
    f = figure(figsize = (4.5, 4.5))
    plot(1:length(y), repmat([0],length(y),1), color = "black", linewidth = 1)    
    for ii = 1:length(snames)
        plot(1:size(y)[1], y[:,ii], linewidth = 2,label=snames[ii])
    end
    xlabel("t")
    ylabel(snames[1])
    legend()
    savefig(string(name,".pdf"))
    close(f)
end

function errorARMA(αθ,data)
    TT    = length(data)
    α     = αθ[1]
    θ     = αθ[2]
    M     = zeros(TT)
    error = zeros(TT)    
    for tt in 2:TT
        M[tt] = α*data[tt-1] + θ*(error[tt-1])
        error[tt] = data[tt] - M[tt]
    end
    return(error)
end

function g(αθ,data)
    error = errorARMA(αθ,data)
    α     = αθ[1]
    θ     = αθ[2]
    σ2    = var(error)
    TT    = length(data)
    yt    = data[3:end]
    ytm1  = data[2:end-1]
    ytm2  = data[1:end-2]
    γ0    = σ2 * (2*α*θ+θ^2+1)/(1-α^2)
    γ1    = σ2 * (α^2*θ + θ^2*α+θ+α)/(1-α^2)
    γ2    = α * γ1
    barg = [sum(yt.^2)/(TT-2)    - γ0;
            sum(yt.*ytm1)/(TT-2) - γ1;
            sum(yt.*ytm2)/(TT-2) - γ2]
    
    return(barg)
end

function W_oGMM(αθ,data)
    error = errorARMA(αθ,data)
    TT = length(data)
    gbar = g(αθ,data)
    Ωhat = zeros(length(gbar),length(gbar))    
    α     = αθ[1]
    θ     = αθ[2]
    σ2    = var(error)
    γ0    = σ2 * (2*α*θ+θ^2+1)/(1-α^2)
    γ1    = σ2 * (α^2*θ + θ^2*α+θ+α)/(1-α^2)
    γ2    = α * γ1
    for tt = 3:TT
        g_i   = [data[tt]^2          - γ0;
                 data[tt]*data[tt-1] - γ1;
                 data[tt]*data[tt-2] - γ2]
        Ωhat  = copy(Ωhat) + g_i * (g_i.') - gbar * (gbar.')
    end
    return (Ωhat/(TT-2))
end

function g_mom(αθ,data)
    srand(55)
    error  = errorARMA(αθ,data)
    α      = αθ[1]
    θ      = αθ[2]
    σ2     = var(error)
    TT     = length(data)
    yLarge = simulateARMA(α, θ, sqrt(σ2),Integer(1e5))[1]
    k2     = σ2 * (2*α*θ+θ^2+1)/(1-α^2)
    k3     = sum(yLarge.^3)/1e5
    k4     = sum(yLarge.^4)/1e5    
    barg   = [sum(data.^2)/TT    - k2;
              sum(data.^3)/TT    - k3;
              sum(data.^4)/TT    - k4]
    return(barg)
end

function W_oGMM_mom(αθ,data)
    srand(55)    
    error = errorARMA(αθ,data)
    TT = length(data)
    gbar = g_mom(αθ,data)
    Ωhat = zeros(length(gbar),length(gbar))    
    α     = αθ[1]
    θ     = αθ[2]
    σ2    = var(error)
    yLarge = simulateARMA(α, θ, sqrt(σ2),Integer(1e5))[1]    
    k2     = σ2 * (2*α*θ+θ^2+1)/(1-α^2)
    k3     = sum((yLarge-mean(yLarge)).^3)/1e5
    k4     = sum((yLarge-mean(yLarge)).^4)/1e5
    meand  = mean(data)
    for tt = 1:TT
        g_i   = [(data[tt]-meand)^2 - k2;
                 (data[tt]-meand)^3 - k3;
                 (data[tt]-meand)^4 - k4]
        Ωhat  = copy(Ωhat) + g_i * (g_i.') - gbar * (gbar.')
    end
    return (Ωhat/(TT))
end


function oGMM(objective, vcov, data, start)
    ## First, get GMM vcov matrix using identity weighting
    J(theta) = objective(theta, data).' * objective(theta,data)
    αθgmm = Optim.optimize(J,start, BFGS()).minimizer
    
    J_opt(theta) = objective(theta, data).' * inv(vcov(αθgmm, data)) * objective(theta,data)
    αθogmm = Optim.optimize(J_opt,start, LBFGS(), Optim.Options(g_tol = 1e-12))
    sigma  = var(errorARMA(αθogmm.minimizer, data))
    append!(αθogmm.minimizer, sigma)
    return αθogmm.minimizer
end

function momentFinder()
    yl,e,el,b,a = symbols("y_{t-1},e_{t},e_{t-1},\\beta,\\alpha", real=true)
    y = a*yl + b * el + e
    write("k3.tex", SymPy.latex(expand(y^3)))
    write("k4.tex", SymPy.latex(expand(y^4)))
end
    

    

##----------------------------------------------------------------------------##
##---------------------------- Simulation ------------------------------------##
##----------------------------------------------------------------------------##
srand(6413)
y,e = simulateARMA(0.8, 0.5, sqrt(0.5),500)

plotts(y,"arma", "y")

##----------------------------------------------------------------------------##
##---------------------------- GMM: autocorr ---------------------------------##
##----------------------------------------------------------------------------##

αθogmm = oGMM(g, W_oGMM, y, [0.5, 0.5])

textable(["\$\\alpha\$", "\$\\theta\$", "\$\\sigma^2\$"],
         αθogmm,
         precision = "%4.4f", fname = "ogmm_ac")

##----------------------------------------------------------------------------##
##---------------------------- GMM: centered moments -------------------------##
##----------------------------------------------------------------------------##

αθogmm_mom = oGMM(g_mom, W_oGMM_mom, y, [0.55, 0.5])


textable(["\$\\alpha\$", "\$\\theta\$", "\$\\sigma^2\$"],
         αθogmm_mom,
         precision = "%4.4f", fname = "ogmm_mom")
