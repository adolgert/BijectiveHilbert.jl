using BijectiveHilbert
using Documenter

makedocs(;
    modules=[BijectiveHilbert],
    authors="Andrew Dolgert <adolgert@uw.edu>",
    repo="https://github.com/adolgert/BijectiveHilbert.jl/blob/{commit}{path}#L{line}",
    sitename="BijectiveHilbert.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://adolgert.github.io/BijectiveHilbert.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Usage" => "usage.md",
        "Algorithms" => [
            "simple2d.md",
            "globalgray.md",
            "compact.md"
        ]
    ],
)

deploydocs(;
    devbranch = "main",
    repo="github.com/adolgert/BijectiveHilbert.jl",
    deploy_config=Documenter.GitHubActions()
)
