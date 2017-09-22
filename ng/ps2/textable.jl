## Creates the inside of a table that can be imported into tex.
## If no filename is given, returns the string to be printed.
## 
## Inputs
##   labels    : A string array of row labels
##   data      : Data to be printed to the table
##   precision : Either a string describing the format of every
##               element in the table, or a table of such strings,
##               the same size as 'data', specifying a format
##               for each variable
##   fname     : The name of the output file

## Note: I did nothing to check for bad inputs... So you get out
##       what you put in. Also, this doesn't yet support column
##       names. 

using Formatting
function textable(labels, data; precision = "%.2g", fname = "")
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


    ## Write to a file if you want. 
    if fname != ""
        write(string(fname, ".tex"), outstr)
    end

    return(outstr)
end


