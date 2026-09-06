@testitem "Aqua.jl Quality Assurance" begin
    using Aqua
    using PsiTK

    # test_all encompasses ambiguities, unbound args, undefined exports, etc.
    Aqua.test_all(PsiTK)
end
