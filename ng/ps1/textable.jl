using Formatting
function textable(labels, data, precision = "%.2g")
    outstr = ""

    ## Get information about size of table 
    vector = length(size(data)) == 1
    nrows = size(data)[1]
    if vector
        ncols = 1
    else
        ncols = size(data)[2]
    end

    ## Precision

    if typeof(precision) == String
        pre_mat = repmat([precision], nrows, ncols)
    else
        pre_mat = copy(precision)
    end

    ## Write the thing
    linesep = "\\\\ \n "
    colsep = "\&"    
    for row in 1:nrows
        outstr = string(outstr, labels[row])
        for col in 1:ncols
            formatted = sprintf1(pre_mat[row,col], data[row,col])
            outstr = string(outstr, colsep, "\$ ", formatted,   "\$")
        end
        outstr = string(outstr, linesep)
    end
    return(outstr)
end


