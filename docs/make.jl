using Documenter
using PsiTK 

makedocs(;
    modules = [PsiTK],
    authors = "Tobias Schäfer",
    sitename = "PsiTK.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://toschaefer.github.io/PsiTK.jl",
        edit_link = "main",
        logo = "logo/PsiTK.png", 
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => "manual.md"
    ],
)

deploydocs(
    repo = "github.com/toschaefer/PsiTK.jl.git",
)
