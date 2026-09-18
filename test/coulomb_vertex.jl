@testitem "Coulomb Vertex Computation" setup=[TestSystems] tags=[:coulomb_vertex] begin
    using LinearAlgebra

    scfres = TestSystems.setup_water_hf(n_bands_converge=8)
    space = OrbitalSpace(scfres)
    nkpt = length(scfres.basis.kpoints)
    nbands = scfres.n_bands_converge

    fitting = compute_coulomb_vertex(space; n_bands = nbands)
    (; Γ, G_vectors, kernel_fourier) = fitting
    @test isnothing(fitting.singular_vectors)

    # Dimensions (nkpt, nbands, nkpt, nbands, nG_reduced)
    @test size(Γ)[1:4] == (nkpt, nbands, nkpt, nbands)
    @test size(Γ, 5) == length(G_vectors) == length(kernel_fourier) > 0

    # Fingerprint for regression testing
    @test isapprox(norm(Γ), 2.319359783263448, rtol = 1e-6)

    # Bare overlap densities: Γ = √v ⊙ ρ on the same G vectors
    ρmnG, G_vectors_ρ = compute_overlap_densities(space; n_bands = nbands, Ecut_ratio = 2/3)
    @test G_vectors_ρ == G_vectors
    @test ρmnG .* reshape(sqrt.(kernel_fourier), 1, 1, 1, 1, :) ≈ Γ

    # Orthonormality: ρ_mn(G=0) ∝ δ_mn
    iG0 = findfirst(iszero, G_vectors)
    ρ0 = ρmnG[1, :, 1, :, iG0]
    @test ρ0 ≈ ρ0[1, 1] * I(nbands) atol = 1e-8

    # Hermiticity: ρ_nm(-G) = conj(ρ_mn(G))
    G_to_idx = Dict(G => i for (i, G) in enumerate(G_vectors))
    idx_minus_G = [G_to_idx[-G] for G in G_vectors]
    @test ρmnG[1, 2, 1, 1, idx_minus_G] ≈ conj.(ρmnG[1, 1, 1, 2, :])

    # Callback is called once per unique orbital pair (upper triangle for symmetric spaces)
    steps = Int[]
    compute_overlap_densities(space; n_bands = nbands, callback = info -> push!(steps, info.step))
    @test steps == 1:(nbands * (nbands + 1) ÷ 2)

    # Default Ecut_ratio = 1.0 reproduces the full plane-wave grid of the basis
    ρ_full, G_full = compute_overlap_densities(space; n_bands = nbands)
    @test length(G_full) == length(scfres.basis.kpoints[1].G_vectors)
    @test size(ρ_full, 5) > size(ρmnG, 5)

    # Ecut_ratio = 4 (DFTK's default supersampling of 2) holds the exact overlap densities;
    # beyond the FFT grid an error is raised instead of silently truncating
    ρ_exact, G_exact = compute_overlap_densities(space; n_bands = nbands, Ecut_ratio = 4)
    @test size(ρ_exact, 5) == length(G_exact) > size(ρ_full, 5)
    @test all(G -> -G in G_exact, G_exact)
    @test_throws ErrorException compute_overlap_densities(space; n_bands = nbands, Ecut_ratio = 8)

    # The Cc4s dump refuses non-orthonormal (non-canonical) spaces
    space_nonortho = OrbitalSpace(space.basis, space.ψ, space.eigenvalues, space.occupations, space.εF, false)
    @test_throws ErrorException dump_cc4s_files(space_nonortho, fitting; folder = mktempdir())

    # CoulombGramian compression: Γ_F = Γ_G U with the returned singular vectors
    fitting_cg = compress_coulomb_vertex(fitting, CoulombGramian(thresh = 1e-3))
    val_cg = norm(fitting_cg.Γ)
    @test isapprox(val_cg, 2.318834927236268, rtol = 1e-6)
    @test size(fitting_cg.Γ)[1:4] == (nkpt, nbands, nkpt, nbands)
    NG, NF = size(fitting_cg.singular_vectors)
    @test NG == size(Γ, 5) && NF == size(fitting_cg.Γ, 5) < NG
    Γmat = reshape(Γ, :, NG)
    @test reshape(Γmat * fitting_cg.singular_vectors, size(fitting_cg.Γ)) ≈ fitting_cg.Γ
    @test fitting_cg.G_vectors === G_vectors && fitting_cg.kernel_fourier === kernel_fourier

    # AdaptiveRandomizedSVD compression
    fitting_svd = compress_coulomb_vertex(fitting, AdaptiveRandomizedSVD(thresh = 1e-3))
    val_svd = norm(fitting_svd.Γ)
    # The randomized subspace can only lose spectral weight w.r.t. the exact Gramian result
    @test val_svd <= val_cg * (1 + 1e-10)
    @test isapprox(val_svd, val_cg, rtol = 1e-3)
    @test size(fitting_svd.Γ)[1:4] == (nkpt, nbands, nkpt, nbands)
    @test size(fitting_svd.Γ, 5) < size(Γ, 5)

    # Compressing twice accumulates the transformations
    fitting_cg2 = compress_coulomb_vertex(fitting_cg, CoulombGramian(thresh = 1e-2))
    @test size(fitting_cg2.singular_vectors, 1) == NG
    @test reshape(Γmat * fitting_cg2.singular_vectors, size(fitting_cg2.Γ)) ≈ fitting_cg2.Γ
end
