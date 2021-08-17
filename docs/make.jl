using Documenter
using Quon
using DocThemeIndigo

indigo = DocThemeIndigo.install(Quon)

makedocs(;
    modules = [Quon],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical="https://QuantumBFS.github.io/Quon.jl",
        assets=String[indigo],
    ),
    pages = [
        "Home" => "index.md",
    ],
    repo = "https://github.com/QuantumBFS/Quon.jl",
    sitename = "Quon.jl",
)

deploydocs(; repo = "github.com/QuantumBFS/Quon.jl")
