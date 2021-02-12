using Test
#using Documenter
using Distributions, LogNormals 

include("fitstats.jl")

include("normal.jl")
include("lognormal.jl")
include("normal.jl")
include("logitnormal.jl")

include("distributionvector.jl")
include("sumlognormals.jl")



# TODO move to documentation
if (false) # only interactively
    using StatsPlots
    plot(d); 
    plot!(dn, linetype = :line)
    vline!([mean(d)])
end

# testing doctests
# make sure to not test for error. This does not work in test, because error compromises former output
# better keep docu and tests separated, but test execution of examples too.
#DocMeta.setdocmeta!(LogNormals, :DocTestSetup, :(using Distributions,LogNormals); recursive=true)
#doctest(LogNormals, manual = false)


