using Documenter
using PsiTK 


using Documenter
using PsiTK 

makedocs(
    sitename = "PsiTK.jl",
    modules  = [PsiTK],
    authors = "Tobias Schäfer",
    checkdocs = :exports,
    format = Documenter.HTML(;
        sidebar_sitename = false, 
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => "manual.md",
        "API Reference" => "api.md",
    ]
)

deploydocs(
    repo = "github.com/toschaefer/PsiTK.jl.git",
)
