using Documenter
using Dante

makedocs(
    sitename = "Dante",
    format = Documenter.HTML(),
    modules = [Dante],
    pages = [
        "Home" => "index.md",
        "User Guide" => Any[
            "Getting Started" =>  "getting_started.md",
            "Input Parameters" => "input.md"
            ],
        "Method" => Any[
            "MHD" => "mhd.md",
            "Finite Volume" => "fv.md",
            "Standard Test" => "test.md"
            ],
        "API Reference" => "api.md"
    ]
)

deploydocs(;
    repo="github.com/henry2004y/Dante.jl",
)
