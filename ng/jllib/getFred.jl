using TimeSeries
##----------------------------------------------------------------------------##
function getFred(seriesname; freq = "q", agg = mean)
    fname = string(seriesname, ".txt")
    if !isfile(fname)
        fredurl = string("https://fred.stlouisfed.org/data/",
                         seriesname, ".txt")
        download(fredurl, string(seriesname, ".txt"))
    end
    lines = readlines(open(fname))

    top = maximum([contains(line, "VALUE") for line in lines].*(1:length(lines)))
    
    dates  = [Date(split(line)[1], "y-m-d") for line in lines[top+1:end]]
    values = [float(split(line)[2]) for line in lines[top+1:end]]

    TS     = TimeArray(dates, values)
    TSQ = []
    if freq == "q"
        for qq = 1:4
            Q = when(TS, quarterofyear, qq)
            Q = collapse(Q, year, first, agg)
            if qq == 1
                TSQ = Q
            else
                TSQ = vcat(TSQ, Q)
            end
        end
        TS = TSQ
    end
    return(TS)
end
