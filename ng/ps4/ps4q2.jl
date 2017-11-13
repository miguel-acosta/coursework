using Distributions, PyPlot, Combinatorics
srand(6413)
##----------------------------------------------------------------------------##
##                                 FUNCTIONS                                  ##
##----------------------------------------------------------------------------##
function plotdist(x,y;label="",mode=false,xlab="θ")
    ymax, xmax = findmax(y)
    offsety = 0.01*(maximum(y)-minimum(y))
    offsetx = 0.01*(maximum(x)-minimum(x))    
    plot(x,y,label = label)
    plot(x[xmax], ymax, "ko")
    text(x[xmax]+offsetx, ymax+offsety,
         string(mode ? " mode = " : " max = ", x[xmax]))
    xlabel(xlab)
    ylabel(mode ? "Posterior" : "Log Likelihood")
    tight_layout()
end


##----------------------------------------------------------------------------##
## Question 2, part 1
##----------------------------------------------------------------------------##
Nd = 50
β = 0.9
σ2 = 1
τ = 1/σ2

x = rand(Uniform(10,20),Nd)
e = rand(Normal(0,σ2),Nd)

y = x*β + e

l(b) = -sum((y-x*b).^2)
plotdist(0:0.001:2, l.(0:0.001:2))
savefig("output/q2i.pdf")
close()

##----------------------------------------------------------------------------##
## Question 2, part 2
##----------------------------------------------------------------------------##
n = 100
yb2 = rand(Binomial(n,0.2),Nd)
yb5 = rand(Binomial(n,0.5),Nd)
l1(θ) = sum(yb2)*log(θ) + sum(n-yb2)*log(1-θ)
l2(θ) = sum(yb5)*log(θ) + sum(n-yb5)*log(1-θ)
plotdist(0.01:0.001:.99, l1.(0.01:0.001:.99),label="θ=0.2")
plotdist(0.01:0.001:.99, l2.(0.01:0.001:.99),label="θ=0.5")
legend()
savefig("output/q2ii.pdf")
close()


θ = 0.5
s = 7
nck(NN,KK) = length(combinations(1:NN,KK))
l3(NN) = log(nck(NN,s)) + s*log(θ) + (NN-s) * log(1-θ)
ngrid = collect(7:1:30)
l3eval = l3.(ngrid)
plotdist(ngrid, l3eval,xlab = "n")
savefig("output/q2iiN.pdf")
close()


##----------------------------------------------------------------------------##
## Question 2, part 3
##----------------------------------------------------------------------------##
n = 5
s = 2
postConj(θ,a,b) = (a-1+s) * log(θ) + (b-1+n-s) * log(1-θ)
postJeff(θ)     = s*log(θ) + (n-s) * log(1-θ) + log(sqrt(1/θ + 1/(1-θ)))

θgrid = 0.01:0.001:0.99
plotdist(θgrid,postConj.(θgrid,1,1),mode=true,label="Conjugate, 1")
plotdist(θgrid,postConj.(θgrid,3,3),mode=true,label="Conjugate, 3")
plotdist(θgrid,postJeff.(θgrid),mode=true,label="Jeffrey's")
legend()
savefig("output/q2iii.pdf")
close()
