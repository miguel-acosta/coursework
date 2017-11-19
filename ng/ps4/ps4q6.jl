using Distributions

n     = 400
npar  = 3
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

    fit = x * b
    fit = fit.*(fit.>0)

    return(-sum(abs.(y-fit)))
end
##----------------------------------------------------------------------------##
function qbe(objective, b, DATA; burn = 2e4, maxit = burn + 2e4,
             lb = b-5, ub = b+5, lv = 0.01, t = 1, vm = sqrt(.15*abs(b)),
             c = 0.5*ones(length(b)), updateFreq = 200)

    # Setup
    npar = length(b)
    
    # Counters and storage
    nup = nreject = new = ndown = nout = 0
    naccept = zeros(npar,1)
    store = zeros(niter, npar)

    for jj in 1:niter
        for ii = 1:npar

            xb[ii] = b[ii] - rand(Normal()) * sqrt(vm[i])
            if (xb[ii] < lb[ii]) | (xb[ii] > ub[ii])
                xb[ii] = b[ii]
                nout += 1
            end

            fb = objective(xb,DATA)

            if fb >= f
                f = fb
                b = xb
                naccept[ii] += 1
            else
                [SOME STUFF]
            end
        end ## end loop over parameters
        store[jj,:] = b

        ## Update variance
        if (jj % 200) == 0
            ratio = naccept/ns
            for ii in 1:npar
                vm[ii] = ratio[ii] > .6 ? vm[ii]*(1+c[ii]*((ratio[ii]-.6)/.4)) : vm[ii]
                vm[ii] = ratio[ii] < .6 ? vm[ii]/(1+c[ii]*((.4-ratio[ii])/.4)) : vm[ii]
                vm[ii] = vm[ii] > (ub[ii]-lb[ii]) ? ub[ii] - lb[ii]            : vm[ii]
                vm[ii] = vm[ii] < lv*abs(b[ii])   ? lv[ii]*abs(b[ii])          : vm[ii]
                vm[ii] = vm[ii] > hv*abs(b[ii])   ? abs(b[ii])/hv              : vm[ii]
            end
            nup = nreject = new = ndown = nout = 0
            naccept = zeros(npar,1)
        end
    end

    store0 = store[burn+1:maxit,:]

    mean0 = mean(store0b,1)
    std0  = std (store0b,1)    
    cr5   = cr10 = zeros(2,npar)
    for ii in 1:npar
        cr5 [:,ii]  = [percentile(store0b[:,ii],0.025), percentile(store0b[:,ii],0.975)]
        cr10[:,ii]  = [percentile(store0b[:,ii],0.05),  percentile(store0b[:,ii],0.95 )]        
    end
    

    
    return(mean0, std0, cr5, cr10, ratio)
end
    
    
             
