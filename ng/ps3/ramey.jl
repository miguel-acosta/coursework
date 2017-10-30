using ExcelReaders, TimeSeries, IterableTables, DataFrames
include("../jllib/getFred.jl")
include("VARfuncs.jl")
""" 
From Ramey: 
I use monthly data and include log industrial production, the unemployment rate, the log of the
CPI, and the log of a commodity price index in the first block. The second block consists of the
federal funds rate. The third block consists of the logs of nonborrowed reserves, total reserves
and M1.

"""

##  Pull in FRED Data
IP   = log.(getFred("INDPRO",freq = "m"))
UR   =      getFred("UNRATE",freq = "m")
CPI  = log.(getFred("CPIAUCSL",freq = "m"))
FF   =      getFred("FEDFUNDS",freq = "m")
NBR  = log.(getFred("BOGNONBR",freq = "m")[Date(1959,1,1):Date(2007,12,31)])
TR   = log.(getFred("RESBALNS",freq = "m"))
M1   = log.(getFred("M1SL",freq = "m"))

## Commodity prices from Ramey's database
ramey   = readxl(DataFrame, "Ramey_HOM_monetary/Monetarydat.xlsx", "Monthly!A1:I682", header= true)
LPCOM   = ramey[:LPCOM]
nonNA   = maximum(!isna(LPCOM) .* collect(1:length(LPCOM)))
LPCOM   = Array(LPCOM[1:nonNA])
rameyD  = Date.(Date(ramey[:DATES][1],1,1):Dates.Month(1):Date(2100,1,1))[1:nonNA]
COMMODITY = TimeArray(rameyD, LPCOM)

## Merge
all = merge(IP, merge(UR, merge(CPI, merge(COMMODITY, merge(FF, merge(NBR, merge(TR, M1)))))),
            colnames = ["IP", "UR", "CPI", "COMMODITY", "FF", "NBR", "TR", "M1"])

## Trim samples for CEE replication 
dataSamp1 = DataFrame(all[Date(1965,6,1):Date(1995,6,1)])[2:end]
dataSamp2 = DataFrame(all[Date(1983,1,1):Date(2007,12,1)])[2:end]
dataSamp3 = DataFrame(all[Date(1983,1,1):Date(2007,12,1)])[[:IP,:UR,:CPI,:COMMODITY,:FF]]

## Estimate VARs. 12 lags. It was not easy to figure out that 12 lags were used... 
Aols1,μols1,Σols1,resids1 = VARols(12, dataSamp1)
Aols2,μols2,Σols2,resids2 = VARols(12, dataSamp2)
Aols3,μols3,Σols3,resids1 = VARols(12, dataSamp3)

## Get the choleskies and triangulars
function LDL(X)
    Pchol = chol(X).'
    S  = diagm(diag(Pchol))
    D  = S^2
    L  = X * inv(S)
    return(L,D)
end
P1 = chol(Σols1).'
P2 = chol(Σols2).'
P3 = chol(Σols3).'
L1, D1 = LDL(Σols1)
L2, D2 = LDL(Σols2)
L3, D3 = LDL(Σols3)

## IRFs
Tirf = 48
irf1 = IRF(Aols1, μols1, Tirf, P1  * [0,0,0,0,1,0,0,0])
irf2 = IRF(Aols2, μols2, Tirf, P2  * [0,0,0,0,1,0,0,0])
irf3 = IRF(Aols3, μols3, Tirf, P3  * [0,0,0,0,1])

function IRFcee(irf,ID, vnames)
    for ii = 1:length(vnames)
        f = figure(figsize = (4,4))
        plot(1:Tirf, repmat([0],Tirf,1), color = "black", linewidth = 1)        
        plot(1:Tirf, irf1[ii,:], linewidth = 2)
        title(vnames[ii])
        savefig(string("cee", ID, "_", vnames[ii], ".pdf"))
        close()
    end
end



IRFcee(irf1, 1, ["IP", "UR", "CPI", "COMMODITY", "FF", "NBR", "TR", "M1"])
IRFcee(irf2, 2, ["IP", "UR", "CPI", "COMMODITY", "FF", "NBR", "TR", "M1"])
IRFcee(irf3, 3, ["IP", "UR", "CPI", "COMMODITY", "FF"])
