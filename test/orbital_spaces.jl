@testitem "Orbital Counts" tags=[:orbital_spaces] begin
    using PsiTK

    # Two k-points with 3 and 2 orbitals; the basis is not needed for counting
    ψ = [rand(ComplexF64, 5, 3), rand(ComplexF64, 5, 2)]
    eigenvalues = [zeros(3), zeros(2)]
    space = OrbitalSpace(nothing, ψ, eigenvalues, eigenvalues, 0.0, true)
    @test n_orbitals(space, 1) == 3
    @test n_orbitals(space, 2) == 2
    @test_throws ArgumentError n_orbitals(space)
    @test n_orbitals(select_orbitals(space, [[1, 2], [1, 2]])) == 2
end
