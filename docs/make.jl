using BijectiveHilbert
using Documenter


makedocs(;
    modules=[BijectiveHilbert],
    warnonly = [:missing_docs],
    authors="Andrew Dolgert <adolgert@cmu.edu>",
    repo=Remotes.GitHub("adolgert", "BijectiveHilbert.jl"),
    sitename="BijectiveHilbert.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://computingkitchen.com/BijectiveHilbert.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Usage" => "usage.md",
        "Background" => "hilbert.md",
        "reference.md"
    ],
)

deploydocs(;
    devbranch = "main",
    repo="github.com/adolgert/BijectiveHilbert.jl.git",
    deploy_config=Documenter.GitHubActions()
)
