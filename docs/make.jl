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

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
