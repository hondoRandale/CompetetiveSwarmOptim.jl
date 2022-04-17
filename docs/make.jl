using CompetetiveSwarmOptim
using Documenter

DocMeta.setdocmeta!(CompetetiveSwarmOptim, :DocTestSetup, :(using CompetetiveSwarmOptim); recursive=true)

makedocs(;
    modules=[CompetetiveSwarmOptim],
    authors="Jules Rasetaharison <Jules.Rasetaharison@tutanota.com>",
    repo="https://github.com/hondoRandale/CompetetiveSwarmOptim.jl/blob/{commit}{path}#{line}",
    sitename="CompetetiveSwarmOptim.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://hondoRandale.github.io/CompetetiveSwarmOptim.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/hondoRandale/CompetetiveSwarmOptim.jl",
    devbranch="main",
)
