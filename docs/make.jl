using BijectiveHilbert
using Documenter


makedocs(;
    # modules=[BijectiveHilbert],
    authors="Andrew Dolgert <adolgert@uw.edu>",
    repo=Remotes.GitHub("adolgert", "BijectiveHilbert.jl"),
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
            "hilbert.md",
            "simple2d.md",
            "globalgray.md",
            "compact.md",
            "facecontinuous.md"
        ]
    ],
)

deploydocs(;
    devbranch = "main",
    repo="github.com/adolgert/BijectiveHilbert.jl.git",
    deploy_config=Documenter.GitHubActions()
)
