using ExcelReaders, TimeSeries, IterableTables, DataFrames, PyPlot
include("../jllib/getFred.jl")
include("VARfuncs.jl")
## Ramey's  first database
ramey   = readxl(DataFrame, "Ramey_HOM_monetary/Monetarydat.xlsx", "Monthlydat6996!A1:K337", header= true)
ramey[:DATE] = Date(ramey[:DATES][1],1,1):Dates.Month(1):Date(1996,12,1)
rameyShort = ramey[ramey[:DATE] .>= Date(1969,3,1),:] ## , :RRORIG]]
rameyShortVAR = rameyShort[:,[:LIP, :UNEMP, :LCPI, :LPCOM , :FFR]]

## Ramey's  second database
ramey   = readxl(DataFrame, "Ramey_HOM_monetary/Monetarydat.xlsx", "Monthlydat6907!A1:K469", header= true)
ramey[:DATE] = Date(ramey[:DATES][1],1,1):Dates.Month(1):Date(2007,12,1)
rameyLong = ramey[ramey[:DATE] .>= Date(1969,3,1),:] #, :RRSHOCK]]
rameyLongVAR = rameyLong[:,[:LIP, :UNEMP, :LCPI, :LPCOM , :FFR]]


## Trim samples for CEE replication 

## Estimate VARs. 12 lags. It was not easy to figure out that 12 lags were used... 
Aols1,μols1,Σols1,resids1 = VARols(12, rameyShortVAR)
Aols2,μols2,Σols2,resids2 = VARols(12, rameyLongVAR)


#### Proxy SVAR
function PROXY(resids, b, ramey, vname)
    TT = size(resids)[1] 
    X = [ones(TT,1) resids[:,5]] ## Fed funds rate shock 
    Z = [ones(TT,1) ramey[end-TT+1:end, vname]]   ## RR shock
    Pz = Z/(Z.'*Z)*Z.'
    for ii in 1:4
        y = resids[:,ii]
        β =  (X.'*Pz * X) \ X.' * Pz * y
        b[ii] = β[2]
    end
    return(b)
end
bSHORT = PROXY(resids1, [0,0,0,0,1.0], rameyShort, :RRORIG)
bLONG  = PROXY(resids2, [0,0,0,0,1.0], rameyLong, :RRSHOCK)

Tirf = 48
irfSHORT = IRF(Aols1, μols1, Tirf, bSHORT)
irfLONG  = IRF(Aols2, μols2, Tirf, bLONG)

## scale log variables by 100:
irfSHORT[[1,3,4],:] = irfSHORT[[1,3,4],:] * 100
irfLONG[[1,3,4],:] = irfLONG[[1,3,4],:] * 100

function IRFproxy(irf1, irf2,ID, vnames)
    for ii = 1:length(vnames)
        f = figure(figsize = (4,4))
        plot(1:Tirf, repmat([0],Tirf,1), color = "black", linewidth = 1)        
        plot(1:Tirf, irf1[ii,:], linewidth = 2, label = "Short Sample")
        plot(1:Tirf, irf2[ii,:], linewidth = 2,linestyle = "--",label="Long Sample")        
        title(vnames[ii])
        legend()
        savefig(string("proxyVAR_", ID, "_", vnames[ii], ".pdf"))
        close()
    end
end
IRFproxy(irfSHORT, irfLONG, "COMBINED", ["IP", "UR", "CPI", "COMMODITY", "FF"])
#IRFcee(irfLONG , "LONG", ["IP", "UR", "CPI", "COMMODITY", "FF"])

