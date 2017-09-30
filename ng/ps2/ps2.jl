using PyPlot, Optim, SymPy
include("textable.jl")
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function simulateARMA(α,θ,σ,TT; randns = "create")
    e = zeros(TT)
    y = zeros(TT)
    if randns == "create"
        randns = randn(TT-1,1) * σ
    end
    e[2:TT] = randns
    for tt in 2:TT
        y[tt] = α * y[tt-1] + e[tt] + θ * e[tt-1]
    end
    return y,e
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
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

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
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

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function g(αθ,data)
    error = errorARMA(αθ,data)
    α     = αθ[1]
    θ     = αθ[2]
    σ2    = var(error)
    TT    = length(data)
    yt    = data[3:end]-mean(data)
    ytm1  = data[2:end-1]-mean(data)
    ytm2  = data[1:end-2]-mean(data)
    γ0    = σ2 * (2*α*θ+θ^2+1)/(1-α^2)
    γ1    = σ2 * (α^2*θ + θ^2*α+θ+α)/(1-α^2)
    γ2    = α * γ1
    barg = [sum(yt.^2)/(TT-2)    - γ0;
            sum(yt.*ytm1)/(TT-2) - γ1;
            sum(yt.*ytm2)/(TT-2) - γ2]
    
    return(barg)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
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
    data  = copy(data)-mean(data)
    for tt = 3:TT
        g_i   = [data[tt]^2          - γ0;
                 data[tt]*data[tt-1] - γ1;
                 data[tt]*data[tt-2] - γ2]
        Ωhat  = copy(Ωhat) + g_i * (g_i.') - gbar * (gbar.')
    end
    return (Ωhat/(TT-2))
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function kappa4(α, β, σ, TT)
    c  = 3 * σ^4 *(β^4 + 4 * α * β^3 + 6 * α^2 * β^2 + 4 * α^3 * β + 1)
    K4 = zeros(TT,1)
    for tt = 2:TT
        K4[tt] = α^4 * K4[tt-1] + c
    end
    return K4
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function g_mom(αθ,data)
    error  = errorARMA(αθ,data)
    α      = αθ[1]
    θ      = αθ[2]
    σ2     = var(error)
    TT     = length(data)
    k2     = σ2 * (2*α*θ+θ^2+1)/(1-α^2)
    k3     = 0
    k4     = kappa4(α, θ, sqrt(σ2), TT)
    data   = copy(data)-mean(data)
    barg   = [sum(data.^2    - k2);
              sum(data.^3    - k3);
              sum(data.^4    - k4)]/TT
    return(barg)
end


##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function W_oGMM_mom(αθ,data)
    error = errorARMA(αθ,data)
    TT = length(data)
    gbar = g_mom(αθ,data)
    Ωhat = zeros(length(gbar),length(gbar))    
    α     = αθ[1]
    θ     = αθ[2]
    σ2    = var(error)
    k2     = σ2 * (2*α*θ+θ^2+1)/(1-α^2)
    k3     = 0
    k4     = kappa4(α, θ, sqrt(σ2), TT)
    meand  = mean(data)
    for tt = 1:TT
        g_i   = [(data[tt]-meand)^2 - k2;
                 (data[tt]-meand)^3 - k3;
                 (data[tt]-meand)^4 - k4[tt]]
        Ωhat  = copy(Ωhat) + g_i * (g_i.') - gbar * (gbar.')
    end
    return (Ωhat/(TT))
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function g_mom_ac(αθ,data)
    error  = errorARMA(αθ,data)
    α      = αθ[1]
    θ      = αθ[2]
    σ2     = var(error)
    TT     = length(data)
    k2     = σ2 * (2*α*θ+θ^2+1)/(1-α^2)
    k3     = 0
    k4     = kappa4(α, θ, sqrt(σ2), TT)[3:end]
    γ1    = σ2 * (α^2*θ + θ^2*α+θ+α)/(1-α^2)
    γ2    = α * γ1
    meand = mean(data[3:end])
    yt    = data[3:end]-meand
    ytm1  = data[2:end-1]-meand
    ytm2  = data[1:end-2]-meand
    data  = data[3:end]-meand
    barg   = [sum(data.^2    - k2);
              sum(data.^3    - k3);
              sum(data.^4    - k4);
              sum(yt.*ytm1   - γ1);
              sum(yt.*ytm2   - γ2)]/(TT-3)
    return(barg)
end


##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function W_oGMM_mom_ac(αθ,data)
    error = errorARMA(αθ,data)
    TT = length(data)
    gbar = g_mom_ac(αθ,data)
    Ωhat = zeros(length(gbar),length(gbar))    
    α     = αθ[1]
    θ     = αθ[2]
    σ2    = var(error)
    k2     = σ2 * (2*α*θ+θ^2+1)/(1-α^2)
    k3     = 0
    k4     = kappa4(α, θ, sqrt(σ2), TT)
    meand  = mean(data)
    γ1    = σ2 * (α^2*θ + θ^2*α+θ+α)/(1-α^2)
    γ2    = α * γ1    
    for tt = 3:TT
        g_i   = [(data[tt]-meand)^2 - k2;
                 (data[tt]-meand)^3 - k3;
                 (data[tt]-meand)^4 - k4[tt]
                 (data[tt]-meand)*(data[tt-1]-meand) - γ1;
                 (data[tt]-meand)*(data[tt-2]-meand) - γ2]
        Ωhat  = copy(Ωhat) + g_i * (g_i.') - gbar * (gbar.')
    end
    return (Ωhat/(TT-3))
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function g_mom_ac_no3(αθ,data)
    error  = errorARMA(αθ,data)
    α      = αθ[1]
    θ      = αθ[2]
    σ2     = var(error)
    TT     = length(data)
    k2     = σ2 * (2*α*θ+θ^2+1)/(1-α^2)
    k4     = kappa4(α, θ, sqrt(σ2), TT)[3:end]
    γ1    = σ2 * (α^2*θ + θ^2*α+θ+α)/(1-α^2)
    γ2    = α * γ1
    meand = mean(data[3:end])
    yt    = data[3:end]-meand
    ytm1  = data[2:end-1]-meand
    ytm2  = data[1:end-2]-meand
    data  = data[3:end]-meand
    barg   = [sum(data.^2    - k2);
              sum(data.^4    - k4);
              sum(yt.*ytm1   - γ1);
              sum(yt.*ytm2   - γ2)]/(TT-3)
    return(barg)
end


##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function W_oGMM_mom_ac_no3(αθ,data)
    error = errorARMA(αθ,data)
    TT = length(data)
    gbar = g_mom_ac_no3(αθ,data)
    Ωhat = zeros(length(gbar),length(gbar))    
    α     = αθ[1]
    θ     = αθ[2]
    σ2    = var(error)
    k2    = σ2 * (2*α*θ+θ^2+1)/(1-α^2)
    k4    = kappa4(α, θ, sqrt(σ2), TT)
    meand = mean(data)
    γ1    = σ2 * (α^2*θ + θ^2*α+θ+α)/(1-α^2)
    γ2    = α * γ1    
    for tt = 3:TT
        g_i   = [(data[tt]-meand)^2 - k2;
                 (data[tt]-meand)^4 - k4[tt]
                 (data[tt]-meand)*(data[tt-1]-meand) - γ1;
                 (data[tt]-meand)*(data[tt-2]-meand) - γ2]
        Ωhat  = copy(Ωhat) + g_i * (g_i.') - gbar * (gbar.')
    end
    return (Ωhat/(TT-3))
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
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

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function momentFinder()
    yl,e,el,b,a = symbols("y_{t-1},e_{t},e_{t-1},\\beta,\\alpha", real=true)
    y = a*yl + b * el + e
    write("k3.tex", SymPy.latex(expand(y^3)))
    write("k4.tex", SymPy.latex(expand(y^4)))
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function HACse(y,x,β,q)
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
    
    return(sqrt.(diag(V)), V)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function gbarMD(αθσ2, ϕ1ϕ2σ2hat)
    α  = αθσ2[1]
    θ  = αθσ2[2]
    σ2 = αθσ2[3]

    ϕ1  = α + θ
    ϕ2  = - α * θ
    σ2v = σ2 * (1+θ^4)

    return([ϕ1 ϕ2 σ2v] - ϕ1ϕ2σ2hat)
end

    

    
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function gbarSMM(αθσ2, shocks, ψsamp)
    α  = αθσ2[1]
    θ  = αθσ2[2]
    σ2 = αθσ2[3]
    S = size(shocks)[2]
    TT = size(shocks)[1] + 1    
    ψ = zeros(3,1)
    for s = 1:S
        sim = simulateARMA(α,θ,[],TT; randns = shocks[:,s] * sqrt(σ2))[1]
        m   = mean(sim[3:end])
        ψ[1] = copy(ψ[1]) + sum((sim[3:end]-m).*(sim[3:end-0]-m))/(TT-2)
        ψ[2] = copy(ψ[2]) + sum((sim[3:end]-m).*(sim[2:end-1]-m))/(TT-2)
        ψ[3] = copy(ψ[3]) + sum((sim[3:end]-m).*(sim[1:end-2]-m))/(TT-2)
    end
    return(ψsamp - ψ/S)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function AR2toARMA11(X, β)
    proj = (X.' * X)
    αAUX =  (proj * β)[2]/(proj * β)[1]
    θsym = Sym("theta")
    θAUX = Float64(minimum(solve(proj[1,1]/proj[1,2] -
                                 (2*αAUX*θsym + θsym^2 + 1)/
                                 (αAUX^2*θsym + θsym^2*αAUX + θsym+αAUX), θsym)))
    sig2 = (proj[1,1]/size(X)[1] * (1-αAUX^2))/(2 * αAUX * θAUX + θAUX^2 + 1)
    return([αAUX, θAUX, sig2])
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function gbarII(αθσ2, shocks, βsamp)
    α  = αθσ2[1]
    θ  = αθσ2[2]
    σ2 = αθσ2[3]
    S = size(shocks)[2]
    TT = size(shocks)[1] + 1    
    βols = zeros(length(βsamp),1)
    for s = 1:S
        sim = simulateARMA(α,θ,[],TT; randns = shocks[:,s] * sqrt(σ2))[1]
        X     = [sim[2:end-1] sim[1:end-2]]
        βols  = copy(βols) + (X.' * X)\X.' * sim[3:end]
    end
    return(βols/S - βsamp)
end

##----------------------------------------------------------------------------##
##---------------------------- Simulation ------------------------------------##
##----------------------------------------------------------------------------##
srand(6413)
y,e = simulateARMA(0.8, 0.5, sqrt(0.5),500)

plotts(y,"arma", ["y"])

##----------------------------------------------------------------------------##
##---------------------------- GMM: autocorr ---------------------------------##
##----------------------------------------------------------------------------##
αθogmm = oGMM(g, W_oGMM, y, [0.5, 0.5])

textable(["\$\\alpha\$", "\$\\beta\$", "\$\\sigma^2\$"],
         αθogmm,
         precision = "%4.4f", fname = "ogmm_ac")

##----------------------------------------------------------------------------##
##---------------------------- GMM: centered moments -------------------------##
##----------------------------------------------------------------------------##
αθogmm_mom = oGMM(g_mom, W_oGMM_mom, y, [0.8, 0.5])

textable(["\$\\alpha\$", "\$\\beta\$", "\$\\sigma^2\$"],
         αθogmm_mom,
         precision = "%4.4f", fname = "ogmm_mom")

##----------------------------------------------------------------------------##
##---------------------------- GMM: centered moments + autocorr --------------##
##----------------------------------------------------------------------------##
αθogmm_mom_ac = oGMM(g_mom_ac, W_oGMM_mom_ac, y, [0.8, 0.5])

textable(["\$\\alpha\$", "\$\\beta\$", "\$\\sigma^2\$"],
         αθogmm_mom_ac,
         precision = "%4.4f", fname = "ogmm_mom_ac")

##----------------------------------------------------------------------------##
##-----------------GMM: centered moments + autocorr, no 3rd moment -----------##
##----------------------------------------------------------------------------##
αθogmm_mom_ac_no3 = oGMM(g_mom_ac_no3, W_oGMM_mom_ac_no3, y, [0.8, 0.5])

textable(["\$\\alpha\$", "\$\\beta\$", "\$\\sigma^2\$"],
         αθogmm_mom_ac_no3,
         precision = "%4.4f", fname = "ogmm_mom_ac_no3")


##----------------------------------------------------------------------------##
##---------------------------- Auxiliary Model: OLS --------------------------##
##----------------------------------------------------------------------------##
X     = [y[2:end-1] y[1:end-2]]
βols  = (X.' * X)\X.' * y[3:end]
σ2v   = sum((y[3:end] - X * βols).^2)/(length(y)-2) #variance of the errors
SEols = HACse(y[3:end], X, βols, 5)[1]
ϕ1    = sum(βols)
ϕ1se  = sqrt.([1 1] * HACse(y[3:end], X, βols, 5)[2] * [1;1])
write("aux_ols.tex", string("The estimate of \$\\phi(1)\$ is ",
                            round(ϕ1,2), " ", round.(ϕ1se,3), "."))

textable(["\$\\phi\_1\$", "\$\\phi\_2\$"],
         [βols SEols],
         precision = "%4.4f", fname = "aux_ols_full")


textable(["\$\\alpha\$", "\$\\beta\$", "\$\\sigma^2\$"],
         AR2toARMA11(X, βols),
         precision = "%4.4f", fname = "aux_ols_trans")
    

##----------------------------------------------------------------------------##
##---------------------------- Auxiliary Model: MD ---------------------------##
##----------------------------------------------------------------------------##
# I don't know what this means............. 
J_auxMD(αθσ2) = sum(gbarMD(αθσ2,[βols[1] βols[2] σ2v]).^2)
αθσ_auxMD = Optim.optimize(J_auxMD,[0.8, 0.5, 0.5], BFGS())
textable(["\$\\alpha\$", "\$\\beta\$", "\$\\sigma^2\$"],
         αθσ_auxMD.minimizer,
         precision = "%4.4f", fname = "aux_md")

##----------------------------------------------------------------------------##
##---------------------------- SMM -------------------------------------------##
##----------------------------------------------------------------------------##
S  = 501
TT = 500
shocks = randn(TT-1, S)
m = mean(y[3:end])
ACsamp = [sum((y[3:end]-m).*(y[3:end-0]-m))/(TT-2)
          sum((y[3:end]-m).*(y[2:end-1]-m))/(TT-2)
          sum((y[3:end]-m).*(y[1:end-2]-m))/(TT-2)]
J(αθσ2) = sum(gbarSMM(αθσ2, shocks, ACsamp).^2)
γsmm = Optim.optimize(J,[0.8, 0.5, 0.5], BFGS())

textable(["\$\\alpha\$", "\$\\beta\$", "\$\\sigma^2\$"],
         γsmm.minimizer,
         precision = "%4.4f", fname = "smm")


##----------------------------------------------------------------------------##
##------------------------INDIRECT INFERENCE ---------------------------------##
##----------------------------------------------------------------------------##
JII(αθσ2) = sum(gbarII(αθσ2, shocks, βols).^2)
γii = Optim.optimize(JII,[0.8, 0.5, 0.5], BFGS())
textable(["\$\\alpha\$", "\$\\beta\$", "\$\\sigma^2\$"],
         γii.minimizer,
         precision = "%4.4f", fname = "ii")
