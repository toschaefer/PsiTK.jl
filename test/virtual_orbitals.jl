@testitem "Generate Orbitals" setup=[TestSystems] begin
    using PsiTK
    using DFTK
    using LinearAlgebra

    scfres = TestSystems.setup_h_chain_hf(; Ecut=15)
    basis = scfres.basis
    occ_space = select_orbitals(OrbitalSpace(scfres), OccupiedOrbitals())
    Nfull = length(basis.kpoints[1].G_vectors)
    Nocc = size(occ_space.ψ[1], 2)
    N_virt_target = 12

    # ---------------------------------------------------------
    # 1. Test CanonicalVirtuals (FullDiag)
    # ---------------------------------------------------------
    solver_fd = FullDiagonalizationEigensolver()
    
    # Test N_virt = :all
    target_canon_all = CanonicalVirtuals(:all)
    virt_canon_fd_all = generate_orbitals(target_canon_all, scfres, occ_space, solver_fd)
    
    @test size(virt_canon_fd_all.ψ[1]) == (Nfull, Nfull - Nocc)
    @test virt_canon_fd_all.ψ[1]' * virt_canon_fd_all.ψ[1] ≈ I
    @test norm(occ_space.ψ[1]' * virt_canon_fd_all.ψ[1]) < 1e-6
    @test isapprox(norm(virt_canon_fd_all.ψ[1]), 26.400757564888178, rtol=1e-6)
    @test isapprox(sum(virt_canon_fd_all.eigenvalues[1]), 6318.532948394055, rtol=1e-6)
    # Check Diagonalization Invariant (off-diagonals are zero)
    H_v = virt_canon_fd_all.ψ[1]' * (scfres.ham[1] * virt_canon_fd_all.ψ[1])
    @test norm(H_v - Diagonal(H_v)) < 1e-6
    @test diag(H_v) ≈ virt_canon_fd_all.eigenvalues[1]

    # Test N_virt
    target_canon = CanonicalVirtuals(N_virt_target)
    virt_canon_fd = generate_orbitals(target_canon, scfres, occ_space, solver_fd)
    
    @test size(virt_canon_fd.ψ[1]) == (Nfull, N_virt_target)
    @test virt_canon_fd.ψ[1]' * virt_canon_fd.ψ[1] ≈ I
    @test norm(occ_space.ψ[1]' * virt_canon_fd.ψ[1]) < 1e-6
    @test isapprox(norm(virt_canon_fd.ψ[1]), 3.464101615137754, rtol=1e-6)
    @test isapprox(sum(virt_canon_fd.eigenvalues[1]), 8.215464748722031, rtol=1e-6)
    
    # ---------------------------------------------------------
    # 2. Test CanonicalVirtuals (LOBPCG)
    # ---------------------------------------------------------
    solver_lobpcg = LOBPCGEigensolver(tol=1e-7, maxiter=500)
    virt_canon_lobpcg = generate_orbitals(target_canon, scfres, occ_space, solver_lobpcg)
    
    @test size(virt_canon_lobpcg.ψ[1]) == (Nfull, N_virt_target)
    @test virt_canon_lobpcg.ψ[1]' * virt_canon_lobpcg.ψ[1] ≈ I
    @test norm(occ_space.ψ[1]' * virt_canon_lobpcg.ψ[1]) < 1e-6
    @test isapprox(norm(virt_canon_lobpcg.ψ[1]), 3.464101615137754, rtol=1e-6)
    @test isapprox(sum(virt_canon_lobpcg.eigenvalues[1]), 8.215464748722068, rtol=1e-6)
    # Parity check: LOBPCG and FullDiag should produce the same eigenvalues for CanonicalVirtuals
    @test isapprox(virt_canon_lobpcg.eigenvalues[1], virt_canon_fd.eigenvalues[1], rtol=1e-4)

    # ---------------------------------------------------------
    # 3. Test DensitySpecificVirtuals
    # ---------------------------------------------------------
    target_dsv = DensitySpecificVirtuals(N_virt_target)
    
    # Run once
    virt_dsv_1 = generate_orbitals(target_dsv, scfres, occ_space, solver_lobpcg)
    
    @test size(virt_dsv_1.ψ[1]) == (Nfull, N_virt_target)
    @test virt_dsv_1.is_orthonormal == false
    @test norm(occ_space.ψ[1]' * virt_dsv_1.ψ[1]) < 1e-6
    
    # Canonicalizing the DSVs should restore exact orthonormality and set the flag
    virt_dsv_canon = canonicalize_orbitals(virt_dsv_1, scfres.ham)
    @test virt_dsv_canon.is_orthonormal == true
    @test virt_dsv_canon.ψ[1]' * virt_dsv_canon.ψ[1] ≈ I
end
