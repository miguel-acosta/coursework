include("../jllib/textable.jl")
function summaryPlots(α,t100,t098,t090)
    ## Plot the distribution of α
    KDEαOLS = kde(α[:,1])
    KDEαGLS = kde(α[:,2])
    plot(KDEαOLS.x, KDEαOLS.density, label = "OLS")
    plot(KDEαGLS.x, KDEαGLS.density, label = "GLS")
    xlabel("α")
    legend()
    title(string("Distribution of α"))
    savefig("output/alph.pdf")
    close()

    ## Plot the distribution of t-statistics for the various
    ## tests: OLS
    tsims = [t100, t098, t090]
    tnames = ["1", "0.98", "0.90"]
    for tt in 1:length(tsims)
        KDEtOLS = kde(tsims[tt][:,1])
        plot(KDEtOLS.x, KDEtOLS.density,
             label = string("OLS, α = ", tnames[tt]))
    end
    xlabel("t")
    legend()
    title(string("Distribution of t-stats: OLS"))
    savefig("output/tOLS.pdf")
    close()

    ## Plot the distribution of t-statistics for the various
    ## tests: OLS
    for tt in 1:length(tsims)
        KDEtGLS = kde(tsims[tt][:,2])
        plot(KDEtGLS.x, KDEtGLS.density,
             label = string("GLS, α = ", tnames[tt]))        
    end
    xlabel("t")
    legend()
    title(string("Distribution of t-stats: GLS"))
    savefig("output/tGLS.pdf")
    close()

    ## Plot the distribution of t-statistics for the various
    ## tests: OLS and GLS
    for tt in 1:length(tsims)
        KDEtOLS = kde(tsims[tt][:,1])
        KDEtGLS = kde(tsims[tt][:,2])
        plot(KDEtOLS.x, KDEtOLS.density,
             label = string("OLS, α = ", tnames[tt]))
        plot(KDEtGLS.x, KDEtGLS.density,
             label = string("GLS, α = ", tnames[tt]))
        xlabel("t")
        legend()
        title(string("Distribution of t-stats"))
        savefig(string("output/t", tnames[tt], ".pdf"))
        close()            
    end

end

function critval(t100, t098, t090, suffix; GLS = true)
    tsims = [t100, t098, t090]
    tnames = ["1.00", "0.98", "0.90"]
    critOLS = zeros(3)
    critGLS = zeros(3)    
    for tt in 1:length(tsims)
        critOLS[tt] = quantile(abs.(tsims[tt][:,1]),0.05)
        critGLS[tt] = quantile(abs.(tsims[tt][:,2]),0.05)
    end

    dataOut = GLS ? [critOLS critGLS] : critOLS
    textable(tnames, dataOut; precision = "%.4g", fname = string("output/tdist_", suffix))
end
