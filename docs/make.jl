using PolyaGammaSamplers
using Documenter

DocMeta.setdocmeta!(PolyaGammaSamplers, :DocTestSetup, :(using PolyaGammaSamplers); recursive=true)

makedocs(;
    modules=[PolyaGammaSamplers],
    authors="Iván Gutiérrez <ivangutierrez1988@gmail.com> and contributors",
    repo="https://github.com/igutierrezm/PolyaGammaSamplers.jl/blob/{commit}{path}#{line}",
    sitename="PolyaGammaSamplers.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://igutierrezm.github.io/PolyaGammaSamplers.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Getting Started" => "getting_started.md",
        "Samplers" => "samplers.md",
    ],
)

deploydocs(;
    repo="github.com/igutierrezm/PolyaGammaSamplers.jl",
)
