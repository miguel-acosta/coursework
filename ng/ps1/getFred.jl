function getFred(seriesname)
    fname = string(seriesname, ".txt")
    if !isfile(fname)
        fredurl = string("https://fred.stlouisfed.org/data/", seriesname, ".txt")
        download(fredurl, string(seriesname, ".txt"))
    end
    lines = readlines(open(fname))
    close(f)
    top = maximum([contains(line, "VALUE") for line in lines].*(1:length(lines)))
    
    dates  = [Date(split(line)[1], "y-m-d") for line in lines[top+1:end]]
    values = [float(split(line)[2]) for line in lines[top+1:end]]

    return(dates,values)
end
