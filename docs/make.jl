using StylizedFacts
using Documenter

DocMeta.setdocmeta!(StylizedFacts, :DocTestSetup, :(using StylizedFacts); recursive=true)

makedocs(;
    modules=[StylizedFacts],
    authors="aaron-wheeler",
    repo="https://github.com/aaron-wheeler/StylizedFacts.jl/blob/{commit}{path}#{line}",
    sitename="StylizedFacts.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://aaron-wheeler.github.io/StylizedFacts.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/aaron-wheeler/StylizedFacts.jl",
    devbranch="main",
)
