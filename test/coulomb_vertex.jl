@testitem "Coulomb Vertex Computation" setup=[TestSystems] tags=[:coulomb_vertex] begin
    using LinearAlgebra

    scfres = TestSystems.setup_water_hf(n_bands_converge=8)
    space = OrbitalSpace(scfres)

    ΓmnG = compute_coulomb_vertex(space; n_bands = scfres.n_bands_converge)

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

    # Test CoulombGramian compression
    cg_alg = CoulombGramian(thresh = 1e-3)
    ΓmnG_cg = compress_coulomb_vertex(ΓmnG, cg_alg)
    val_cg = norm(ΓmnG_cg)
    @test isapprox(val_cg, 2.3182425676526193, rtol = 1e-6)
    @test size(ΓmnG_cg)[1:4] == (nkpt, nbands, nkpt, nbands)
    @test size(ΓmnG_cg, 5) < size(ΓmnG, 5)

    # Test AdaptiveRandomizedSVD compression
    svd_alg = AdaptiveRandomizedSVD(thresh = 1e-3)
    ΓmnG_svd = compress_coulomb_vertex(ΓmnG, svd_alg)
    val_svd = norm(ΓmnG_svd)
    @test isapprox(val_svd, val_cg, rtol = 1e-6)
    @test size(ΓmnG_svd)[1:4] == (nkpt, nbands, nkpt, nbands)
    @test size(ΓmnG_svd, 5) < size(ΓmnG, 5)
end
