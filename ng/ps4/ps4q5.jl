using Distributions
srand(649)

## Generate data 
N = 100
x = rand(Uniform(0,1),N)
e = rand(Normal(0,0.5),N)
y = x * 1 + e

## Calculations 
X       = [ones(N) x]
βhat    = (X.'*X) \ X.' * y
s2      = (y-X*βhat).' * (y-X*βhat)/(N-2)
k       = length(βhat)

## Function for posterior
function post(βprior, Vprior, vprior, s2prior)
    
    ## Posteriors
    Vpost  = inv(inv(Vprior) + X.'* X)
    βpost  = Vpost*(X.'*X * βhat + inv(Vprior)*βprior)
    vpost  = vprior + N
    vspost = vprior*s2prior + (N-k) * s2 + (βhat - βprior).' *inv(Vprior + inv(X.'*X)) *(βhat - βprior)

    return(βpost, Vpost * vspost/(vpost - 2))
end


## Find posterior: NG
βpost, varβpost = post([0,1], eye(2), 1,1)

## Find posterior: Flat
βpostFlat, varβpostFlat = post([0,1], 6413 * eye(2), 0,1)
