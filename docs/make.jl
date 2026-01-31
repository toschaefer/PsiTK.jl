using Documenter
using PsiTK 

makedocs(
    sitename = "PsiTK.jl",
    modules  = [PsiTK],
    pages = [
        "Home" => "index.md",
        "Manual" => "manual.md"
    ]
)

deploydocs(
    repo = "github.com/toschaefer/PsiTK.jl.git",
)
