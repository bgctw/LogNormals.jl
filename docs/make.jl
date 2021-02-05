using Documenter
using LogNormals, Distributions

push!(LOAD_PATH,"../src/")
makedocs(sitename="LogNormals.jl",
         pages = [
            "Home" => "index.md",
            "Fit to statistic" => "fitstats.md",
            "Sum LogNormals" => "sumlognormals.md",
            #"misc" => "misc.md",
         ],
         format = Documenter.HTML(prettyurls = false)
)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/bgctw/LogNormals.jl.git",
    devbranch = "main"
)
