using StatsBase, Distributions, PyPlot, Optim
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

sk = skewness(x)
ku = kurtosis(x)

# Bera-Jarque Test
jb = (N/6)*(sk^2 + (1/4)*(ku-3)^2)

# Internet says JB is chi-squared(2). So, what level do we reject at?
α = 1-cdf(Chisq(2), jb)

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


#------------------------------------------------------------------------------#
#----------------------------------- Q1.d -------------------------------------#
#------------------------------------------------------------------------------#
include("newtonRaphson.jl")
g(μ,σ2) = 1/N * [sum(x-μ), sum((x-μ).^2 - σ2),sum((x-μ).^3),
                     sum((x-μ).^4 - 3 * σ2^2)]

J(theta) = g(theta[1],theta[2]).' * g(theta[1],theta[2])

θhat,exitStatus = newtonRaphson(J, [0.02,0.001],1e-7,1e5)
# Optim.optimize(J, [0.02,0.001])

#------------------------------------------------------------------------------#
#----------------------------------- Q2 ---------------------------------------#
#------------------------------------------------------------------------------#
α = 0.8
θ = 0.5
σ = sqrt(0.05)
e = zeros(N)
y = zeros(N)
e[2:N] = randn(N-1,1) * σ
y[2:N] = [α * y[tt-1] + e[tt] + θ * e[tt-1] for tt in 2:N]

#------------------------------------------------------------------------------#
#----------------------------------- Q2.i -------------------------------------#
#------------------------------------------------------------------------------#
function m(αθ)
    α = αθ[1]
    θ = αθ[2]
    M = zeros(N)
    e = zeros(N)    
    for tt in 2:N
        M[tt] = α*y[tt-1] + θ*(e[tt-1])
        e[tt] = y[tt] - M[tt]
    end
    return(sum(e.^2))
end
    
    
Optim.optimize(m,[0.9,0.9], BFGS())


#------------------------------------------------------------------------------#
#----------------------------------- Q2.iii -----------------------------------#
#------------------------------------------------------------------------------#
include("getFred.jl")
dates, cpi = getFred("CPIAUCSL")
Π = [log(cpi[tt]) - log(cpi[tt-12]) for tt in 13:length(cpi)]
