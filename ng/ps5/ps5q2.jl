using Distributions, PyPlot, KernelDensity

srand(6413)
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function genSampI(ϕ, TT) 
    u = rand(Normal(0,0.5),TT)
    e = rand(Normal(0,1),TT)
    β = 1
    x = copy(u)
    for tt in 2:TT
        x[tt] = ϕ * x[tt-1] + u[tt]
    end
    y = β * x + e
    return(y,x)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function genSampII(ϕ, TT) 
    u = rand(Normal(0,0.5),TT)
    v = rand(Normal(0,0.5),TT)
    e = copy(v)
    for tt in 2:TT
        e[tt] = e[tt-1] + v[tt]
    end
    β = 1
    x = copy(u)
    for tt in 2:TT
        x[tt] = ϕ * x[tt-1] + u[tt]
    end
    y = β * x + e
    return(y,x)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function sampleCoefficients(ϕ, TT, part)
    if part == 1
        y,x = genSampI(ϕ, TT)
    else
        y,x = genSampII(ϕ, TT)
    end
    X = [ones(TT) 1:TT x]
    βhat = (X.' * X) \ X.' * y
    return(βhat[2], βhat[3])
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function finiteSampleDistribution(ϕ, part; NN=10000, TT=200)
    srand(6413)
    γdist = zeros(NN)
    βdist = zeros(NN)
    for nn in 1:NN
        γsamp, βsamp = sampleCoefficients(ϕ, TT)
        γdist[nn] = γsamp
        βdist[nn] = βsamp
    end
    KDEγ = kde(γdist)
    KDEβ = kde(βdist)    
    return(mean(γdist), mean(βdist), KDEγ.x, KDEγ.density, KDEβ.x, KDEβ.density)
    
end

##----------------------------------------------------------------------------##
##----------------------------- PART I ---------------------------------------##
##----------------------------------------------------------------------------##
for ϕi in [0.8, 0.9, 0.95, 0.99, 1]
    meanγ, meanβ, xγ,yγ, xβ, yβ = finiteSampleDistribution(ϕi,1)
    ## Plot β
    plot(xβ,yβ/sum(yβ),label = string("ϕ = ", ϕi))

end
xlabel("β")
legend()
title(string("Distribution of β"))
savefig("output/bet_kde_1.pdf")
close()

for ϕi in [0.8, 0.9, 0.95, 0.99, 1]
    meanγ, meanβ, xγ,yγ, xβ, yβ = finiteSampleDistribution(ϕi,1)
    ## Plot γ
    plot(xγ,yγ/sum(yγ),label = string("ϕ = ", ϕi))
end
xlabel("γ")
legend()
title(string("Distribution of γ"))
savefig("output/gam_kde_1.pdf")
close()

##----------------------------------------------------------------------------##
##----------------------------- PART II --------------------------------------##
##----------------------------------------------------------------------------##
for ϕi in [-0.8, 0.8]
    meanγ, meanβ, xγ,yγ, xβ, yβ = finiteSampleDistribution(ϕi,2)
    ## Plot β
    plot(xβ,yβ/sum(yβ),label = string("ϕ = ", ϕi))
end
xlabel("β")
legend()
title(string("Distribution of β"))
savefig("output/bet_kde_2.pdf")
close()

for ϕi in [-0.8, 0.8]
    meanγ, meanβ, xγ,yγ, xβ, yβ = finiteSampleDistribution(ϕi,2)
    ## Plot β
    plot(xγ,yγ/sum(yγ),label = string("ϕ = ", ϕi))
end
xlabel("γ")
legend()
title(string("Distribution of γ"))
savefig("output/gam_kde_2.pdf")
close()
