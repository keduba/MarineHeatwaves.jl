using MarineHeatwaves
using Documenter

DocMeta.setdocmeta!(MarineHeatwaves, :DocTestSetup, :(using MarineHeatwaves); recursive=true)

makedocs(
    modules = [MarineHeatwaves],
    authors = "keduba <jfk.uba@pm.me>, Manal Elawady <mh.elawady@uliege.be> and contributors",
    sitename = "MarineHeatwaves.jl",
    format = Documenter.HTML(;
        canonical="https://keduba.github.io/MarineHeatwaves.jl",
        edit_link="main",
        assets=String[],
    ),
    pages = [
        "Home" => "index.md",
        "Example" => "example.md",
        "Extremes" => [
            "Calculating extremes" => "mextremes.md",
            "Detecting events" => "mevents.md",
        ],
        "Metrics" => "metrics.md"
    ],
    remotes = nothing,
)

deploydocs(repo="github.com/keduba/MarineHeatwaves.jl",
    devbranch = "main",
)
