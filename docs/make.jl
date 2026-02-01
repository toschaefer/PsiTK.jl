using Documenter
using PsiTK 


using Documenter
using PsiTK 

makedocs(
    sitename = "PsiTK.jl",
    modules  = [PsiTK],
    authors = "Tobias Schäfer",
    pages = [
        "Home" => "index.md",
        "Manual" => "manual.md"
    ]
)

deploydocs(
    repo = "github.com/toschaefer/PsiTK.jl.git",
)
