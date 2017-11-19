using Distributions, PyPlot
include("../jllib/textable.jl")
srand(6413)

dat = readdlm("input/probit.dat")
y = dat[:,1]
X = [ones(y) dat[:,2]]

βprior = [0,0]
g      = 10 
yprior = zeros(y)
σ      = 1

nburn = 1000
nsamp = 9000

ndat       = size(X)[1]
k          = size(X)[2]
N          = nburn + nsamp
βsamp      = zeros(N,k)
ysamp      = zeros(N,ndat)
βsamp[1,:] = βprior
ysamp[1,:] = yprior
xxi        = 
B          = inv(g*inv(X.'*X) + X.'*X)
Ny1        = sum(y .> 0.5) ## afraid numerical precision with == 
Ny0        = ndat = Ny1

for ii in 2:N
    βpost = B * ((X.'*X) * βprior/g + X.'* ysamp[ii-1,:])
    βsamp[ii,:] = rand(MvNormal(βpost,B * σ^2))

    y1 = rand.(Truncated.(Normal.(X * βsamp[ii,:],1), 0, Inf),1)
    y0 = rand.(Truncated.(Normal.(X * βsamp[ii,:],1), -Inf,0),1)

    yy1 = [yy[1]  for yy in y1] # Drawing from truncated normal 
    yy0 = [yy[1]  for yy in y0] # this way return array of lists
    
    ysamp[ii,y.>0.5] = yy1[y.>0.5]
    ysamp[ii,y.<0.5] = yy0[y.<0.5]

end

βsamp = βsamp[nburn:end,:]
ysamp = ysamp[nburn:end,:]

βmean = mean(βsamp,1)
βstd  = std(βsamp,1)
ymean = mean(ysamp,1)

textable(["\$\\beta\_0\$", "\$\\beta\_1\$"],
         [βprior βmean.' βstd.'];
         fname = "output/betaprobit")

write("output/probprobit.tex", string(round(cdf.(Normal(0,1), βmean * [1,2])[1],2)))

plot(X[y.>0.5,2],ymean[y.>0.5],  "o",label = "y=1")
plot(X[y.<0.5,2],ymean[y.<0.5],  "o",color = "red",label = "y=0")
xlabel("x")
ylabel("y")
legend()
savefig("output/q4.pdf")
close()
