include("../jllib/getFred.jl")
using IterableTables, DataFrames, PyPlot

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function VARols(P, variables)
    V        = length(variables)
    TT       = length(variables[1])
    A        = zeros(V,V,P)
    μ        = zeros(V,1)

    for depvar = 1:V
        ## Get y variable
        y = variables[depvar][(P+1):TT]

        ## Create x matrix
        X = zeros(TT-P, P*V+1)
        X[:,1] = 1
        col = 2
        for pp = 1:P
            for vv = 1:V
                x = variables[vv]
                X[:,col] = x[(P-pp+1):(TT-pp)]
                col += 1
            end
        end
        β = (X.' * X) \ X.' * y
        μ[depvar] = β[1]

        for pp = 1:P
            inds = (V*(pp-1) + 2):(V*(pp) + 1)
            print(inds)
            A[depvar,:,pp] = β[inds]
        end
    end
        return(A, μ)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function moments(Y)
    LY = Y[1:(end-1),:]
    Y  = Y[2:end,:]
    TT = size(Y)[1]

    EY = repmat(mean(Y,1), TT,1)
    print(size(Y))
    print(size(EY))
    Γ0 = (Y-EY).' * (Y-EY) /TT
    Γ1 = (Y-EY).' * (LY-EY) /TT
    return EY[1,:], Γ0, Γ1
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function VMA(A,μ,j,shock)
    N = length(μ)
    P = size(A)[3]
    Y = zeros(N, j)
    Y[:,1] = shock # + μ
    for jj = 2:j
        y = zeros(N,1) #μ
        for kk = max(1,jj-P):(jj-1)
#            print(string(kk, "--", jj-kk,"\n"))
            y += A[:,:,(jj-kk)] * Y[:,kk]
        end
        Y[:,jj] = y
    end
#    Y += repmat(μ,1,j)
    return(Y)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function plotts(y, name, snames)
    N = size(y)[2]
    T = size(y)[1]    
    print(N)
    f = figure(figsize = (12,4))
    for ii = 1:length(snames)
        subplot(100 + N*10 + ii)
        plot(1:T, repmat([0],T,1), color = "black", linewidth = 1)
        plot(1:T, y[:,ii], linewidth = 2,label=snames[ii])
        xlim(0,T)
        xlabel("t")
        ylabel(snames[ii])        
    end
    savefig(string(name,".pdf"))
    close(f)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## ESTIMATE THE VAR
gdp = getFred("GDP",freq = "q")
M3  = getFred("MABMM301USM189S",freq = "q", agg = mean)
FF  = getFred("FEDFUNDS",freq = "q", agg = mean)

all = merge(gdp, merge(M3,FF),colnames = ["GDP", "M3", "FF"])[1:end-1]

Y      = DataFrame(all)[2:end]
T      = size(Y)[1]
Y[2:T,1] = 400*(log.(Y[2:T,1]) - log.(Y[1:(T-1),1]))
Y[2:T,2] = 400*(log.(Y[2:T,2]) - log.(Y[1:(T-1),2]))
Y = Y[2:T,:]

Aols,μols = VARols(2,Y)

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## Compute moments
μ, Γ0, Γ1 = moments(Array(Y))

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## MA Representation 
Φ = VMA(Aols,μols,100,[0,0,1])
plotts(Φ.', "shockFFR", ["GDP", "M3", "FFR"])
Φ = VMA(Aols,μols,100,[0.01,0,0])
plotts(Φ.', "shockM3", ["GDP", "M3", "FFR"])
Φ = VMA(Aols,μols,100,[0.01,0,0])
plotts(Φ.', "shockGDP", ["GDP", "M3", "FFR"])

