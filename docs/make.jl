using AstroUtils
using Documenter

DocMeta.setdocmeta!(AstroUtils, :DocTestSetup, :(using AstroUtils); recursive=true)

makedocs(;
    modules=[AstroUtils],
    authors="Grant Hecht",
    repo="https://github.com/GrantHecht/AstroUtils.jl/blob/{commit}{path}#{line}",
    sitename="AstroUtils.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://GrantHecht.github.io/AstroUtils.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/GrantHecht/AstroUtils.jl",
    devbranch="main",
)
