## Gradient of f at x0
function gradient(f,x0,eps = 1e-6)
    n    = length(x0)
    grad = Array{Float64}(n)
    for ii = 1:n
        plus = zeros(n)
        plus[ii] = eps
        grad[ii] = (f(x0 + plus) - f(x0))/eps
    end
    return(grad)
end

## Hessian 
function hessian(f,x0,eps = 1e-6)
    n    = length(x0)
    H = Array{Float64}(n,n)
    for ii = 1:n
        for jj = 1:n
            plus = zeros(n)
            plus[jj] = eps            
            H[ii,jj] = (gradient(f,x0+plus,eps)[ii] - gradient(f,x0,eps)[ii])/eps
        end
    end
    return(H)
end

## Newton Raphson 
function newtonRaphson(f, x0, tol=1e-4, maxiter=1e3)
    guess = x0 
    diff  = sum(abs(f(guess)))
    iter  = 0
    while diff > tol && iter < maxiter
        guess = guess - hessian(f,guess)\gradient(f,guess)
        diff = sum(abs(f(guess)))
        iter = iter + 1
    end
    if diff <= tol
        stopReason = "Tolerance satisfied."
    elseif iter >= maxiter
        stopReason = "Maximum number of iterations achieved"
    else
        stopReason = "This is not doing what you want it to, Miguel."
    end
    return(guess,stopReason)
end    
