using Distributions,PyPlot
srand(6413)


##----------------------------------------------------------------------------##
## Question 3, part 1
##----------------------------------------------------------------------------##
invDoubleExp(u) = log(2*u) * (u < 1/2) - log(2*(1-u)) * (u >= 1/2)

function drawNormDubExp(N)
    b = 1
    dubExpDraws = zeros(N)
    niter = 1
    while b <= N && niter <  (N*10)
        x      = invDoubleExp(rand(Uniform(0,1)))
        u      = rand(Uniform(0,1))
        thresh = exp(-.5 * x^2 + abs(x) - .5)
        if u <= thresh
            dubExpDraws[b] = x
            b += 1
        end
        niter += 1
    end
    return(dubExpDraws,niter)
end

nde200,niter200   = drawNormDubExp(200)
nde10000,niter10000 = drawNormDubExp(10000)

write("output/dubExpAccept.tex", string(niter200))
PyPlot.plt[:hist](nde200,50,normed = true,color="dodgerblue")
title("200 Normal draws: Accept/Reject Double Exponential")
savefig("output/q3i.pdf")
close()

write("output/dubExpAcceptBig.tex", string(niter10000))
PyPlot.plt[:hist](nde10000,100,normed = true,color="maroon")
title("10,000 Normal draws: Accept/Reject Double Exponential")
savefig("output/q3iBig.pdf")
close()


##----------------------------------------------------------------------------##
## Question 3, part 2
##----------------------------------------------------------------------------##
function drawTruncNormBF(N,μmin;μ=0,σ=1)
    b = 1
    truncNormDraws = zeros(N)
    niter = 1
    while b <= N && niter <  (N*10)
        x      = rand(Normal(μ,σ))
        if x >= μmin
            truncNormDraws[b] = x
            b += 1
        end
        niter += 1
    end
    return(truncNormDraws,niter)
end


function drawTruncNormICDF(N,μmin;μ=0,σ=1)
    u = rand(Uniform(0,1),N)
    CDFn(x)  = cdf(Normal(μ,σ), x)
    ICDFn(u) = quantile(Normal(μ,σ), u)
    x = σ * ICDFn.(u * (1-CDFn((μmin-μ)/σ)) + CDFn((μmin-μ)/σ)) + μ
    return(x)
end

function drawTruncNormAR(N,μmin)
    b = 1
    truncNormDraws = zeros(N)
    Ginv(u) = μmin - log(1-u)
    niter = 1
    while b <= N && niter <  (N*10)
        x      = Ginv.(rand(Uniform(0,1)))
        u      = rand(Uniform(0,1))
        thresh = exp(-(x^2)/2 + x-μmin - 1.5)
        if u <= thresh
            truncNormDraws[b] = x
            b += 1
        end
        niter += 1
    end
    return(truncNormDraws,niter)
end



μbar   = -1
Ndraws = 1000
## Brute Force
drawsBF, niterBF = drawTruncNormBF(Ndraws,μbar)
PyPlot.plt[:hist](drawsBF,75,normed = true,color="darkcyan")
write("output/q3iiBF.tex", string(niterBF))
title("1000 Truncated Normal Draws: Brute Force")
savefig("output/q3iiBF.pdf")
close()

## The old fashined way: Inverse CDF
drawsCDF         = drawTruncNormICDF(Ndraws,μbar)
PyPlot.plt[:hist](drawsCDF,75,normed = true,color="red")
title("1000 Truncated Normal Draws: Inverse CDF")
savefig("output/q3iiICDF.pdf")
close()

### Accept reject
## First, get a feeling for what's going on 
CDFn(x)  = cdf.(Normal(0,1), x)
PDFn(x)  = pdf.(Normal(0,1), x)

f(x) = PDFn(x)/(1-CDFn(-1))
g(x) = exp.(-(x+1)).*(x .>= -1)

θgrid = collect(μbar:0.01:2)
plot(θgrid, f.(θgrid),label="f")
plot(θgrid, g.(θgrid),label="g")
plot(θgrid, f.(θgrid)./g.(θgrid), label = "f/g")
M,ind = findmax(f(θgrid)./g(θgrid))
plot(θgrid, f(θgrid)./(M*g(θgrid)),label = "f/Mg")
plot(θgrid, M*g(θgrid),label = "Mg")
xlabel("θ")
legend()
savefig("output/q3iiifg.pdf")
close()

## Do the drawing: 
drawsAR,niterAR = drawTruncNormAR(Ndraws,μbar)
write("output/q3iiAR.tex", string(niterAR))
PyPlot.plt[:hist](drawsAR,75,normed = true,color="red")
title("1000 Truncated Normal Draws: Accept-Reject")
savefig("output/q3iiIAR.pdf")
close()
##----------------------------------------------------------------------------##
## Question 3, part 3
##----------------------------------------------------------------------------##
σ1       = 0.5
σ2       = 4
γ        = 1
μ1       = 10
μ2       = 0
ρ        = γ/(σ1*σ2)
θ1bounds = [9.5,100]
θ2bounds = [-10,3]
Ndraws   = 1000
burnin   = 250


N = burnin + Ndraws
draws  = zeros(N, 2)
draws[1,:] = [μ1,μ2]
b = 2
niter = 1
while b <= N && niter <  (N*10)
    candidate1 = rand(Normal(μ1 + (σ1/σ2)* ρ * (draws[b-1,2]-μ2),
                             sqrt((1-ρ^2) * σ1^2)))
    candidate2 = rand(Normal(μ2 + (σ2/σ1) * ρ * (candidate1-μ1),
                             sqrt((1-ρ^2) * σ2^2)))
    if candidate1 >= θ1bounds[1] && candidate1 <= θ1bounds[2] &&
       candidate2 >= θ2bounds[1] && candidate2 <= θ2bounds[2]
        draws[b,1] = candidate1
        draws[b,2] = candidate2
        b+=1
    end
    niter+=1
end
draws = draws[burnin:end,:]
plot(draws[:,1],draws[:,2],".",color = "lightcoral")
title("1000 Bivariate Normal Draws: MCMC")
xlabel(L"$\theta_1$")
ylabel(L"$\theta_2$")
savefig("output/q3iii.pdf")
close()

write("output/q3iii.tex",
      string("I set \$\\sigma\_1 = ", σ1, "\$, ",
             "\$\\sigma\_2 = ", σ2, "\$, ",
             "\$\\gamma = ", γ, "\$, ",
             "\$\\mu\_1 = ", μ1, "\$, ",
             "\$\\mu\_2 = ", μ2, "\$, ",
             " and required \$\\theta\_1\\in\[",
             θ1bounds[1], ",", θ1bounds[2],"\]\$", 
             " and \$\\theta\_2\\in\[",
             θ2bounds[1], ",", θ2bounds[2],"\]\$.",
             " It took ", niter, " draws to get the ",
             burnin, " burn-in draws and ", Ndraws,
             " kept-draws."))
             
