using Distributions, PyPlot
include("../jllib/textable.jl")
srand(6413)

n     = 400
npar  = 4
β     = [-6.0; 3.0*ones(npar-1)]
X     = [ones(n) rand(Normal(),n,npar-1)]
U     = rand(Normal(),n).*(X[:,2].^2)
Ystar = X * β + U
Y     = Ystar.*(Ystar.>0)
DATA  = [Y X]
b_ols = (X.'*X)\  X.' *Y

function obj_lad(b,DATA)
    y = DATA[:,1]
    x = DATA[:,2:end]
    n = length(y)
    fit = x * b
    fit = fit.*(fit.>0)

    Ln = y - fit
    return(-sum(abs.(Ln)))
end
##----------------------------------------------------------------------------##
function qbe(objective, b, DATA; burn = Int64(2e4), maxit = Int64(burn + 2e4),
             lb = b-10, ub = b+10, lv = 0.002, hv = 3, t = 1, vm = sqrt.(.15*abs.(b)),
             c = 0.5*ones(length(b)), updateFreq = 200)

    # Setup and priors
    npar = length(b)
    
    # Counters and storage
    naccept        = zeros(npar,1)
    nreject        = zeros(npar,1)    
    xb             = copy(b)
    ratio          = zeros(npar,1)
    store          = zeros(maxit, npar)
    
    f = objective(b,DATA)

    for jj in 1:maxit
        for ii = 1:npar
            ## Draw candidate
            p = Uniform(lb[ii],ub[ii])
            q(m) = Normal(m,sqrt(vm[ii]))
            xb[ii] = b[ii] - rand(q(0))
            
            ## Throw out if it's outside of the bounds
            if (xb[ii] < lb[ii]) || (xb[ii] > ub[ii])
                xb[ii] = copy(b[ii])
            end
            
            ## Calculate acceptance probability/posterior change
            fb = objective(xb,DATA)
            py  = pdf(p,xb[ii])
            px  = pdf(p,b[ii] )
            qxy = pdf(q(xb[ii]), b[ii])
            qyx = pdf(q(b[ii]), xb[ii])            
            ρ   = exp(fb)*py*qxy/(exp(f)*px*qyx)
            
            if ρ >= 1
                f = copy(fb)
                b = copy(xb)
                naccept[ii] += 1
            else
                u = rand(Uniform())
                if u <= ρ
                    f = copy(fb)
                    b = copy(xb)
                    naccept[ii] += 1
                else
                    nreject[ii] += 1
                end
            end
        end ## end loop over parameters
        store[jj,:] = copy(b)

        ## Update variance
        if (jj % updateFreq) == 0
            ratio = copy(naccept)/copy(updateFreq)
            print(ratio)
            print("\n")
            print(vm)
            print("\n")            
            
            for ii in 1:npar
                vm[ii] = ratio[ii] > .6 ? vm[ii]*(1+c[ii]*((ratio[ii]-.6)/.4)) : vm[ii]
                vm[ii] = ratio[ii] < .4 ? vm[ii]/(1+c[ii]*((.4-ratio[ii])/.4)) : vm[ii]
                vm[ii] = vm[ii] > (ub[ii]-lb[ii]) ? ub[ii] - lb[ii]            : vm[ii]
                vm[ii] = vm[ii] < lv*abs(b[ii])   ? lv*abs(b[ii])              : vm[ii]
                vm[ii] = vm[ii] > hv*abs(b[ii])   ? abs(b[ii])/hv              : vm[ii]
            end
            naccept = zeros(npar,1)
            nreject = zeros(npar,1)
        end
    end

    # Summarize it all 
    store0 = store[burn+1:maxit,:]

    mean0 = mean(store0,1)
    median0 = median(store0,1)    
    std0  = std(store0,1)    
    cr5   = cr10 = zeros(2,npar)
    for ii in 1:npar
        cr5[:,ii]   = [quantile(store0[:,ii],0.025), quantile(store0[:,ii],0.975)]
        cr10[:,ii]  = [quantile(store0[:,ii],0.05),  quantile(store0[:,ii],0.95 )]        
    end

    for ii in 1:npar
        print(string(mean0[ii], " (", std0[ii], ") ~~~ [", cr5[1,ii], ",", cr5[2,ii],"]\n"))
    end
    
    return(mean0, median0, std0, cr5, cr10,store, store0, ratio)
end
    
    
             
## Run the thing
βqbe, βqbeMed, SEqbe, CI95, CI90, allb, b_kept, RATIO = qbe(obj_lad, b_ols, DATA)



## Make graphs and tables
textable(["Mean", "Median", "SD", "95 CI Low", "95 CI High", "Acceptance"],
         hcat([βqbe.',βqbeMed.', SEqbe.', CI95[1,:],CI95[2,:], RATIO]...).',
         fname = "output/VH")

PyPlot.plt[:hist](b_kept[:,1],100,normed = true,color="dodgerblue")
title(L"Draws of $\beta_1$")
savefig("output/VHb1.pdf")
close()

PyPlot.plt[:hist](b_kept[:,2],100,normed = true,color="red")
title(L"Draws of $\beta_2$")
savefig("output/VHb2.pdf")
close()

PyPlot.plt[:hist](b_kept[:,3],100,normed = true,color="palegreen")
title(L"Draws of $\beta_3$")
savefig("output/VHb3.pdf")
close()

PyPlot.plt[:hist](b_kept[:,3],100,normed = true,color="lightslategray")
title(L"Draws of $\beta_4$")
savefig("output/VHb4.pdf")
close()

plot(allb[:,1], label = L"$\beta_1$",color = "dodgerblue")
plot(allb[:,2], label = L"$\beta_2$",color = "red")
plot(allb[:,3], label = L"$\beta_3$",color = "palegreen")
plot(allb[:,3], label = L"$\beta_4$",color = "lightslategray")
axvline(x=2e4, label = "Burn-in")
title("Monte Carlo Draws")
legend()
savefig("output/VHMC.pdf")
close()
