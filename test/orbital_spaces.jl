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

@testitem "Save and Load Orbitals" setup=[TestSystems] tags=[:io] begin
    using PsiTK
    using DFTK
    using JLD2

    scfres = TestSystems.setup_h_chain_hf(; Ecut=15)
    basis = scfres.basis
    kernel = only(t.kernel for t in basis.model.term_types if t isa ExactExchange)
    occ_space, _ = split_occupied_virtual(OrbitalSpace(scfres))
    target = DensitySpecificVirtuals(n_orbitals=4, interaction_kernel=kernel)
    solver = LOBPCG(tol=1e-7, maxiter=500, callback=identity)
    dsv_space = generate_orbitals(target, occ_space, scfres.ham; solver)
    merged = merge_spaces(occ_space, dsv_space)
    @test occ_space.is_orthonormal && !merged.is_orthonormal

    mktempdir() do dir
        file = joinpath(dir, "orbitals.jld2")
        for space in (occ_space, merged)
            save_orbitals(file, space)
            @test !isfile(file * ".new")
            loaded = load_orbitals(file, basis)
            @test loaded isa typeof(space)
            @test loaded.basis === basis
            for field in (:ψ, :eigenvalues, :occupation, :εF, :is_orthonormal)
                @test getfield(loaded, field) == getfield(space, field)
            end
        end

        basis_other = PlaneWaveBasis(basis.model; Ecut=20, kgrid=basis.kgrid)
        err = @test_throws ArgumentError load_orbitals(file, basis_other)
        @test occursin("Ecut", err.value.msg)
        @test_throws ArgumentError save_orbitals(joinpath(dir, "orbitals.h5"), merged)
    end
end
