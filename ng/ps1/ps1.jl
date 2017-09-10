using StatsBase
using Distributions
using PyPlot
srand(6413)
#------------------------------------------------------------------------------#
#----------------------------------- Q1 ---------------------------------------#
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
normQuants = quantile(Normal(mean(x),var(x)), quants) # should this be normal 0,1?? 
lambQuants = quantile(vec(x),quants)

ioff()
f = figure()
plot([0,maximum(normQuants)],[0,maximum(normQuants)], color = "black",label = "45 Deg. Line")
scatter(normQuants, lambQuants, color="red", label = "Quantiles")
xlabel("Normal")
ylabel("Generalized Lambda")
title("QQ-Plot of Generalized Lambda and Normal Distributions")
legend()

savefig("qq.pdf")
