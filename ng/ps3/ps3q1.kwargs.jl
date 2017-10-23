include("../jllib/getFred.jl")
import Quandl

function VARols(P; kwargs...)
    variables = Dict(kwargs)
    print(keys(variables))
    print('\n')
    V         = length(variables)
    TT = length(variables["v1"])
    A  = zeros(V,V,P)
    μ  = zeros(V,1)
    
    for depvar = 1:V
        ## Get y variable 
        y = variables[string("v",depvar)][(P+1):TT]

        ## Create x matrix 
        X = zeros(TT-P, P*V+1)
        X[:,1] = 1
        col = 2
        for pp = 1:P
            for vv = 1:V
                x = variables[string("v",vv)]
                X[:,col] = x[(P-p+1):(T-p+1)]
                col += 1
            end
        end
        β = (X.' * X) \ X.' * y
        μ[depvar] = β[1]

        for pp = 1:P
            A[depvar,:,p] = β[(P*(pp-1) +2):(P*pp +1)]
        end
    end
        return(A, μ)
end


gdp = getFred("GDP",freq = "q")
M3  = getFred("MABMM301USM189S",freq = "q", agg = mean)
FF  = getFred("FEDFUNDS",freq = "q", agg = mean)

all = merge(gdp, merge(M3,FF),colnames = ["GDP", "M3", "FF"])[1:end-1]



x = VARols(2; v1 = all["GDP"], v2 = all["M3"], v3 = all["FF"])
