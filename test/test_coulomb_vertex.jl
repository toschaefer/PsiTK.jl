@testitem "Coulomb Vertex Computation" setup=[TestSystems] tags=[:coulomb_vertex] begin
    using LinearAlgebra
    
    scfres = TestSystems.setup_water_hf()
    space = OrbitalSpace(scfres)
    
    ΓmnG = compute_coulomb_vertex(space; n_bands=scfres.n_bands_converge)
    
    # Check dimensions
    nkpt = length(scfres.basis.kpoints)
    nbands = scfres.n_bands_converge
    
    # For a gamma point only calculation (nkpt=1), or generally:
    # Dimensions should be (nkpt, nbands, nkpt, nbands, nG_reduced)
    @test size(ΓmnG)[1:4] == (nkpt, nbands, nkpt, nbands)
    @test size(ΓmnG, 5) > 0 # Some G vectors must exist
    
    # Calculate a fingerprint scalar for regression testing
    val = norm(ΓmnG)

    # Check against the reference value
    @test isapprox(val, 2.3235541514083575, rtol=1e-5)
end
