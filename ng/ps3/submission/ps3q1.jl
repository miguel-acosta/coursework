include("../jllib/getFred.jl")
include("../jllib/textable.jl")
include("VARfuncs.jl")
using IterableTables, DataFrames, PyPlot, HypothesisTests
    
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## SET UP DATA
gdp = getFred("GDPC1",freq = "q")
M3  = getFred("MABMM301USM189S",freq = "q", agg = mean)
FF  = getFred("FEDFUNDS",freq = "q", agg = mean)

all = merge(gdp, merge(M3,FF),colnames = ["GDP", "M3", "FF"])[1:end-1]

Y      = DataFrame(all)[2:end]
T      = size(Y)[1]
Y[2:T,1] = 400*(log.(Y[2:T,1]) - log.(Y[1:(T-1),1]))
Y[2:T,2] = 400*(log.(Y[2:T,2]) - log.(Y[1:(T-1),2]))
Y = Y[2:T,:]
plotts(Y, "data", ["GDP", "M3", "FF"])

n  = size(Y)[2]
P  = 2
T  = size(Y)[1]

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## ESTIMATE THE VAR
Aols,μols,Σols,resids = VARols(P,Y)
textable([], μols, "%.2g", "MUols")
textable([], Aols[:,:,1], "%.2g", "A1ols")
textable([], Aols[:,:,2], "%.2g", "A2ols")
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## Stability
M,C = companion(Aols,μols)
textable([], abs.(eig(C)[1]).', "%.2g", "eigenA")

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## Compute moments
μ, Γ0, Γ1 = moments(Array(Y))
textable([], μ, "%.2g", "muY")
textable([], Γ0, "%.2g", "g0Y")
textable([], Γ1, "%.2g", "g1Y")
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## MA Representation 
plotts(IRF(Aols,μols,20,[0,0,1]).',
       "shockFFR", ["Δ GDP", "Δ M3", "FFR"])
plotts(IRF(Aols,μols,20,[0,1,0]).',
       "shockM3", ["Δ GDP", "Δ M3", "FFR"])
plotts(IRF(Aols,μols,20,[1,0,0]).',
       "shockGDP", ["Δ GDP", "Δ M3", "FFR"])

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## AIC and BIC
K = n
pmax = 8

Σ(m)   = VARols(m,Y;chop=(pmax-m))[3]
AIC(m) = log(det(Σ(m))) + (2/T) * (m * K^2 + K)
BIC(m) = log(det(Σ(m))) + (log(T)/T) * (m * K^2 + K)

AICall = [AIC(m) for m in 1:pmax]
BICall = [BIC(m) for m in 1:pmax]
BICmin = findmin(BICall)[2]
AICmin = findmin(AICall)[2]

textable([], AICall.', "%.2g", "AICall")
textable([], BICall.', "%.2g", "BICall")
write(string("AICmin.tex"), string(AICmin))
write(string("BICmin.tex"), string(BICmin))
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## ESTIMATE THE VARS TWO OTHER WAYS, THEN COMPARE COEFFICIENTS
LY = lagmatrix(Y,2)
YY = LY[:,1:K].'
XX = [LY[:,K+1:end] ones(T-P,1)].'

Π = YY * XX.' * inv(XX * XX.')
π = kron((XX * XX.')\ XX ,eye(K)) * vec(YY)

[[vec(reshape(Aols,K,K*P));μols] vec(Π) π] ## Looks good!
textable([], Π, "%.2g", "PIhat")
textable([], π, "%.2g", "pihat")
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## GET DATA BACK BY SIMULATING WITH COMPANION MATRIX
diffSim = maximum(abs.(Array(Y[3:end,:]) - simulateVAR(Y,K,P,C,M, resids)))
write(string("diffSim.tex"), string(diffSim))
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## FORECAST DATA
fcast = zeros(100,3)
for aa = 1:100
    fcast[aa,:] = forecast(aa, Aols, μols, resids)
end
plotts([Array(Y) ;fcast], "forecast", ["GDP", "M3", "FFR"])

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## PREDICTION ERROR COVARIANCE
Yvar  = zeros(12,3)
M3var = zeros(12,3)
FFvar = zeros(12,3)
for h = 1:12
    fev        = FEV(Aols,μols,h)[2]
    Yvar[h,:]  = fev[1,:]
    M3var[h,:] = fev[2,:]
    FFvar[h,:] = fev[3,:]    
end

plotFEV("GDP", Yvar) 
plotFEV("M3" , M3var)
plotFEV("FF" , FFvar)

fevINF   = FEV(Aols,μols,10000)[2]
textable(["Variance of \$\\Delta\$ GDP","Variance of \$\\Delta\$ M3",
          "Variance of FF"],
          fevINF, "%.2g", "FEVLR")

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## BOOTRSTRAPPING IRFS A LA RUNKLE
b025, t975 = runkle(Y, Aols, μols, resids, K, P, 499, 20)

irf = IRF(Aols,μols,20,[1,0,0]).'
plotIRF(irf[:,1], "Δ GDP", "IRF_gdp_gdp",CIlow = b025[:,1,1], CIhigh=t975[:,1,1],Clevel=95)
plotIRF(irf[:,2], "Δ M3" , "IRF_gdp_m3" ,CIlow = b025[:,2,1], CIhigh=t975[:,2,1],Clevel=95)
plotIRF(irf[:,3], "FF"   , "IRF_gdp_ff" ,CIlow = b025[:,3,1], CIhigh=t975[:,3,1],Clevel=95)

irf = IRF(Aols,μols,20,[0,1,0]).'
plotIRF(irf[:,1], "Δ GDP", "IRF_m3_gdp",CIlow = b025[:,1,2], CIhigh=t975[:,1,2],Clevel=95)
plotIRF(irf[:,2], "Δ M3" , "IRF_m3_m3" ,CIlow = b025[:,2,2], CIhigh=t975[:,2,2],Clevel=95)
plotIRF(irf[:,3], "FF"   , "IRF_m3_ff" ,CIlow = b025[:,3,2], CIhigh=t975[:,3,2],Clevel=95)


irf = IRF(Aols,μols,20,[0,0,1]).'
plotIRF(irf[:,1], "Δ GDP", "IRF_ff_gdp",CIlow = b025[:,1,3], CIhigh=t975[:,1,3],Clevel=95)
plotIRF(irf[:,2], "Δ M3" , "IRF_ff_m3" ,CIlow = b025[:,2,3], CIhigh=t975[:,2,3],Clevel=95)
plotIRF(irf[:,3], "FF"   , "IRF_ff_ff" ,CIlow = b025[:,3,3], CIhigh=t975[:,3,3],Clevel=95)

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## AR(2) for GDP
Lgdp = lagmatrix(Y[1],2)
ygdp = Lgdp[:,1]
Xgdp = [ones(size(Lgdp)[1],1) Lgdp[:,2:end]]
βgdp = (Xgdp.' * Xgdp)\Xgdp.' * ygdp
egdp = ygdp - Xgdp * βgdp
μgdp = βgdp[1]
Agdp = reshape(βgdp[2:end],1,1,2)

fcast_gdp_1_ar2 = [forecast(1,Agdp,μgdp,egdp[1:tt])[1] for tt in 1:length(egdp)]
ar2_fcast_error_1 = fcast_gdp_1_ar2[1:end-1] - Y[1][3:end-1]

fcast_gdp_1_var   = [forecast(1,Aols,μols,resids[1:tt,:])[1] for tt in 1:size(resids)[1]]
var_fcast_error_1 = fcast_gdp_1_var[1:end-1] - Y[1][3:end-1]
plotts([var_fcast_error_1 ar2_fcast_error_1], "fcast_error_gdp",  ["VAR(2)", "AR(2)"];sub=false, LEGEND=true,
       TITLE = "One-period-ahead Forecast Errors for GDP")

TTf = length(var_fcast_error_1)
pval = pvalue(UnequalVarianceTTest(var_fcast_error_1, ar2_fcast_error_1))
write(string("pvalYfcast.tex"), string(pval))
