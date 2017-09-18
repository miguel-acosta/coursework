using StatsBase, Distributions, PyPlot, Optim
include("textable.jl")
srand(6413)
#------------------------------------------------------------------------------#
#----------------------------------- Q1.a -------------------------------------#
#------------------------------------------------------------------------------#


function genLambdaInvCDF(u, λ_1, λ_2, λ_3, λ_4)
    return λ_1 + (u.^λ_3 - (1-u).^λ_4)./λ_2
end



# Generate some uniform random numbers
N = 500
x = genLambdaInvCDF(rand(N,1), 0, -1, -0.0075, -0.03)

sk = skewness(x) + 3 ## pretty sure julia kurtosis is excess. 
ku = kurtosis(x)

# Bera-Jarque Test
jb = (N/6)*(sk^2 + (1/4)*(ku-3)^2)


# Internet says JB is chi-squared(2). So, what level do we reject at?
#α = 100*(1-cdf(Chisq(2), jb))
α = quantile(Chisq(2), 0.95)
#jbout = string("The Bera-Jarque statistic was $jb; this means "
#             "we can reject the null at the $α\% level, since Bera-Jarque "
#             "is distributed \$\Chi^2\_2\$"))
#write("jb.tex", jbout)
jbout = textable(["mean", "variance", "skewness", "kurtosis", "Bera-Jarque", "5\\% Critical Value"],
         [mean(x), var(x), sk, ku, jb, α], ["%4.4f", "%4.4f", "%4.2f", "%4.2f", "%4.2f", "%4.2f"])
write("jb.tex", jbout)         


# qq-plot time...
quants  = collect(0.01:0.01:0.99)
normQuants = quantile(Normal(mean(x),sqrt(var(x))), quants)
lambQuants = quantile(vec(x),quants)

ioff()
f = figure()
ymin = minimum([minimum(normQuants), minimum(lambQuants)])
ymax = maximum([maximum(normQuants), maximum(lambQuants)])
plot([ymin,ymax],[ymin,ymax], color = "black",label = "45 Deg. Line")
scatter(normQuants, lambQuants, color="red", label = "Quantiles")
xlabel("Normal")
ylabel("Generalized Lambda")
title("QQ-Plot of Generalized Lambda and Normal Distributions")
legend()
savefig("qq.pdf")
close(f)

#------------------------------------------------------------------------------#
#----------------------------------- Q1.d -------------------------------------#
#------------------------------------------------------------------------------#
include("newtonRaphson.jl")
g(μ,σ2) = 1/N * [sum(x-μ), sum((x-μ).^2 - σ2),sum((x-μ).^3),
                     sum((x-μ).^4 - 3 * σ2^2)]

J(theta) = g(theta[1],theta[2]).' * g(theta[1],theta[2])

θhat,exitStatus = newtonRaphson(J, [0.02,0.001],1e-7,1e5)
θhat_Optim      = Optim.optimize(J, [0.02,0.001], LBFGS())
nlsout = textable(["\$\\mu\$", "\$\\sigma^2\$"], [θhat  θhat_Optim.minimizer],
                  "%4.4f")
write("nls.tex", nlsout)         

#------------------------------------------------------------------------------#
#----------------------------------- Q2 ---------------------------------------#
#------------------------------------------------------------------------------#
α = 0.8
θ = 0.5
σ = sqrt(0.05)
e = zeros(N)
y = zeros(N)
e[2:N] = randn(N-1,1) * σ
for tt in 2:N
    y[tt] = α * y[tt-1] + e[tt] + θ * e[tt-1]
end

f = figure()
plot(1:N, y, color = "black", label = "y", linewidth = 2 )
#plot(1:N, y-e, linestyle = "dashed", color = "blue", label = "Excluding innovations")
plot(1:N, e, color = "blue", label = "Innovations",linewidth=1)
xlabel("t")
ylabel("y")
legend()
savefig("y.pdf")
close(f)


#------------------------------------------------------------------------------#
#----------------------------------- Q2.i -------------------------------------#
#------------------------------------------------------------------------------#
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

function m(αθ, data)
    error = errorARMA(αθ, data)
    return(error.' * error)
end
    
function g!(grad, αθ, data)
    error = errorARMA(αθ, data)
    grad[2] = -error[2:end].' * error[1:end-1]
    grad[1] = -error[2:end].' * data[1:end-1]
end

m_obj(αθ) = m(αθ, y)
g_obj!(grad, αθ) = g!(grad, αθ, y)
θhat_2 = Optim.optimize(m_obj, g_obj!, [0.5,0.5], LBFGS())
sigma = var(errorARMA(θhat_2.minimizer, y))
append!(θhat_2.minimizer, sigma)
nlsout = textable(["\$\\alpha\$", "\$\\theta\$", "\$\\sigma^2\$"], θhat_2.minimizer,
                  "%4.4f")
write("nls2.tex", nlsout)

#------------------------------------------------------------------------------#
#----------------------------------- Q2.iii -----------------------------------#
#------------------------------------------------------------------------------#
include("getFred.jl")
dates, cpi = getFred("CPIAUCSL")
Π = [log(cpi[tt]) - log(cpi[tt-12]) for tt in 13:length(cpi) if Dates.month.(dates)[tt] == 12]
dates_Π = [date for date in dates[13:end] if Dates.month.(date) == 12]
TT = length(Π)

f = figure()
plot(dates_Π, Π, color = "black", linewidth = 2 )
xlabel("year")
ylabel("Π")
title("Annual Inflation")
savefig("inflation.pdf")
close(f)

## CLS
m_Q2iii(xx) = m(xx,Π)
g_Q2iii!(grad,xx) = g!(grad, xx, Π)
αθ_Q2iii_SA = Optim.optimize(m_Q2iii, g_Q2iii!, [0.5, 0.5], SimulatedAnnealing(),Optim.Options(iterations = 100000))
sigma_SA = var(errorARMA(αθ_Q2iii_SA.minimizer, Π))
append!(αθ_Q2iii_SA.minimizer, sigma_SA)

αθ_Q2iii = Optim.optimize(m_Q2iii, g_Q2iii!, [0.5,0.5], LBFGS())
sigma = var(errorARMA(αθ_Q2iii.minimizer, Π))
append!(αθ_Q2iii.minimizer, sigma)

clsout = textable(["\$\\alpha\$", "\$\\theta\$", "\$\\sigma^2\$"],
                  [αθ_Q2iii_SA.minimizer αθ_Q2iii.minimizer],
                  "%4.4f")
write("cls.tex", clsout)


## (Optimal) GMM
## Need some moment conditions. Impose that e's are serially uncorrelated,
## and that e is uncorrelated lagged y. 
function g(αθ,data)
    error = errorARMA(αθ,data)
    TT = length(data)
    return([(data[1:TT-1].' * error[2:TT])^2, (error[2:TT].' * error[1:TT-1])^2]/(TT-1))
end

function W_oGMM(αθ,data)
    error = errorARMA(αθ,data)
    TT = length(data)
    Ωhat = zeros(2,2)
    gbar = g(αθ,data)
    for tt = 2:TT
        g_i   = [data[tt-1] * error[tt], error[tt] * error[tt-1]]
        Ωhat  = copy(Ωhat) + g_i * (g_i.') - gbar * (gbar.')
    end
    return (Ωhat/(TT-1))
end
        

# First, get GMM vcov matrix using identity weighting
J(theta) = g([theta[1],theta[2]], Π).' * g([theta[1],theta[2]],Π)
αθ_Q2iiiGMM = newtonRaphson(J, [0.5,0.5],1e-12,1e5) #
Optim.optimize(J,[0.5,0.5], BFGS())


J_opt(theta) = g([theta[1],theta[2]], Π).' * inv(W_oGMM(αθ_Q2iiiGMM[1], Π)) * g([theta[1],theta[2]],Π)
αθ_Q2iiiGMMopt = newtonRaphson(J_opt, [0.5,0.5],1e-12,1e5)
sigma = var(errorARMA(αθ_Q2iiiGMMopt[1], Π))
append!(αθ_Q2iiiGMMopt[1], sigma)

αθ_Q2iiiGMMopt_Optim = Optim.optimize(J_opt,[0.5,0.5], LBFGS(), Optim.Options(g_tol = 1e-12))
sigma_Optim = var(errorARMA(αθ_Q2iiiGMMopt[1], Π))
append!(αθ_Q2iiiGMMopt_Optim.minimizer, sigma_Optim)

gmmout = textable(["\$\\alpha\$", "\$\\theta\$", "\$\\sigma^2\$"],
                  [αθ_Q2iiiGMMopt[1] αθ_Q2iiiGMMopt_Optim.minimizer],
                  "%4.4f")
write("ogmm.tex", gmmout)
