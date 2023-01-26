using YangShallowWater
using Documenter

DocMeta.setdocmeta!(YangShallowWater, :DocTestSetup, :(using YangShallowWater); recursive=true)

makedocs(;
    modules=[YangShallowWater],
    authors="Nathanael Wong <natgeo.wong@outlook.com>",
    repo="https://github.com/natgeo-wong/YangShallowWater.jl/blob/{commit}{path}#{line}",
    sitename="YangShallowWater.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://natgeo-wong.github.io/YangShallowWater.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/natgeo-wong/YangShallowWater.jl",
    devbranch="main",
)
