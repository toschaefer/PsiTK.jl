@testitem "Coulomb Vertex Computation" setup=[TestSystems] tags=[:coulomb_vertex] begin
    using LinearAlgebra

    scfres = TestSystems.setup_water_hf(n_bands_converge=8)
    space = OrbitalSpace(scfres)

    ΓmnG, G_vectors, kernel_fourier = compute_coulomb_vertex(space; n_bands = scfres.n_bands_converge)

    # Check dimensions
    nkpt = length(scfres.basis.kpoints)
    nbands = scfres.n_bands_converge

    # For a gamma point only calculation (nkpt=1), or generally:
    # Dimensions should be (nkpt, nbands, nkpt, nbands, nG_reduced)
    @test size(ΓmnG)[1:4] == (nkpt, nbands, nkpt, nbands)
    @test size(ΓmnG, 5) > 0 # Some G vectors must exist

    # Calculate a fingerprint scalar for regression testing
    val = norm(ΓmnG)

    # Check against reference value
    @test isapprox(val, 2.318769223791925, rtol = 1e-6)

    # Bare overlap densities: Γ = √v ⊙ ρ on the same G vectors
    ρmnG, G_vectors_ρ = compute_overlap_densities(space; n_bands = nbands, Ecut_ratio = 2/3)
    @test G_vectors_ρ == G_vectors
    @test ρmnG .* reshape(sqrt.(kernel_fourier), 1, 1, 1, 1, :) ≈ ΓmnG

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

    # Full grid without cutoff reduction
    ρ_full, G_full = compute_overlap_densities(space; n_bands = nbands, Ecut_ratio = nothing)
    @test length(G_full) == length(scfres.basis.kpoints[1].G_vectors)
    @test size(ρ_full, 5) > size(ρmnG, 5)

    # Test CoulombGramian compression
    cg_alg = CoulombGramian(thresh = 1e-3)
    ΓmnG_cg, _ = compress_coulomb_vertex(ΓmnG, cg_alg)
    val_cg = norm(ΓmnG_cg)
    @test isapprox(val_cg, 2.3182425676526193, rtol = 1e-6)
    @test size(ΓmnG_cg)[1:4] == (nkpt, nbands, nkpt, nbands)
    @test size(ΓmnG_cg, 5) < size(ΓmnG, 5)

    # Test AdaptiveRandomizedSVD compression
    svd_alg = AdaptiveRandomizedSVD(thresh = 1e-3)
    ΓmnG_svd, _ = compress_coulomb_vertex(ΓmnG, svd_alg)
    val_svd = norm(ΓmnG_svd)
    @test isapprox(val_svd, val_cg, rtol = 1e-6)
    @test size(ΓmnG_svd)[1:4] == (nkpt, nbands, nkpt, nbands)
    @test size(ΓmnG_svd, 5) < size(ΓmnG, 5)
end
