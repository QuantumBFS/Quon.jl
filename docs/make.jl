using Quon
using Documenter

DocMeta.setdocmeta!(Quon, :DocTestSetup, :(using Quon); recursive=true)

makedocs(;
    modules=[Quon],
    authors="Roger-Luo <rogerluo.rl18@gmail.com> and contributors",
    repo="https://github.com/Roger-luo/Quon.jl/blob/{commit}{path}#{line}",
    sitename="Quon.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Roger-luo.github.io/Quon.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Roger-luo/Quon.jl",
)
