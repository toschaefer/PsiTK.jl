@testitem "Generate Orbitals" setup=[TestSystems] tags=[:virtual_orbitals] begin
    using PsiTK
    using DFTK
    using LinearAlgebra

    scfres = TestSystems.setup_h_chain_hf(; Ecut=15)
    basis = scfres.basis
    ham = scfres.ham
    occ_space, _ = split_occupied_virtual(OrbitalSpace(scfres))
    Nfull = length(basis.kpoints[1].G_vectors)
    Nocc = size(occ_space.ψ[1], 2)
    N_virt_target = 12

    # ---------------------------------------------------------
    # 1. CanonicalVirtuals, FullDiagonalization (default for :all)
    # ---------------------------------------------------------
    virt_canon_fd_all = generate_orbitals(CanonicalVirtuals(), occ_space, ham)

    @test size(virt_canon_fd_all.ψ[1]) == (Nfull, Nfull - Nocc)
    @test virt_canon_fd_all.is_orthonormal
    @test virt_canon_fd_all.ψ[1]' * virt_canon_fd_all.ψ[1] ≈ I
    @test norm(occ_space.ψ[1]' * virt_canon_fd_all.ψ[1]) < 1e-6
    @test isapprox(norm(virt_canon_fd_all.ψ[1]), 26.400757564888178, rtol=1e-6)
    @test isapprox(sum(virt_canon_fd_all.eigenvalues[1]), 6318.26601753155, rtol=1e-6)
    # Diagonalization invariant: the Fock operator is diagonal in the virtual space
    H_v = virt_canon_fd_all.ψ[1]' * (ham[1] * virt_canon_fd_all.ψ[1])
    @test norm(H_v - Diagonal(H_v)) < 1e-6
    @test diag(H_v) ≈ virt_canon_fd_all.eigenvalues[1]

    target = CanonicalVirtuals(n_orbitals=N_virt_target)
    virt_canon_fd = generate_orbitals(target, occ_space, ham; solver=FullDiagonalization())

    @test size(virt_canon_fd.ψ[1]) == (Nfull, N_virt_target)
    @test virt_canon_fd.ψ[1]' * virt_canon_fd.ψ[1] ≈ I
    @test norm(occ_space.ψ[1]' * virt_canon_fd.ψ[1]) < 1e-6
    @test isapprox(norm(virt_canon_fd.ψ[1]), 3.464101615137754, rtol=1e-6)
    @test isapprox(sum(virt_canon_fd.eigenvalues[1]), 8.215419405004843, rtol=1e-6)

    # ---------------------------------------------------------
    # 2. CanonicalVirtuals, LOBPCG (default for a fixed number of orbitals)
    # ---------------------------------------------------------
    solver_lobpcg = LOBPCG(tol=1e-7, maxiter=500, callback=identity)
    virt_canon_lobpcg = generate_orbitals(target, occ_space, ham; solver=solver_lobpcg)

    @test size(virt_canon_lobpcg.ψ[1]) == (Nfull, N_virt_target)
    @test virt_canon_lobpcg.ψ[1]' * virt_canon_lobpcg.ψ[1] ≈ I
    @test norm(occ_space.ψ[1]' * virt_canon_lobpcg.ψ[1]) < 1e-6
    @test isapprox(norm(virt_canon_lobpcg.ψ[1]), 3.464101615137754, rtol=1e-6)
    @test isapprox(sum(virt_canon_lobpcg.eigenvalues[1]), 8.21541940500492, rtol=1e-6)
    # LOBPCG and FullDiagonalization solve the same eigenproblem
    λ_fd = virt_canon_fd.eigenvalues[1]
    @test isapprox(virt_canon_lobpcg.eigenvalues[1], λ_fd, rtol=1e-4)

    # ---------------------------------------------------------
    # 3. DensitySpecificVirtuals with both solvers
    # ---------------------------------------------------------
    target = DensitySpecificVirtuals(n_orbitals=N_virt_target)
    virt_dsv = generate_orbitals(target, occ_space, ham; solver=solver_lobpcg)

    @test size(virt_dsv.ψ[1]) == (Nfull, N_virt_target)
    @test virt_dsv.is_orthonormal == false
    @test norm(occ_space.ψ[1]' * virt_dsv.ψ[1]) < 1e-6

    virt_dsv_fd = generate_orbitals(target, occ_space, ham; solver=FullDiagonalization())
    @test size(virt_dsv_fd.ψ[1]) == (Nfull, N_virt_target)
    @test norm(occ_space.ψ[1]' * virt_dsv_fd.ψ[1]) < 1e-6
    @test isapprox(virt_dsv_fd.eigenvalues[1], virt_dsv.eigenvalues[1], rtol=1e-3)

    # Canonicalizing the DSVs restores exact orthonormality and sets the flag
    virt_dsv_canon = canonicalize_orbitals(virt_dsv, ham)
    @test virt_dsv_canon.is_orthonormal == true
    @test virt_dsv_canon.ψ[1]' * virt_dsv_canon.ψ[1] ≈ I

    # ---------------------------------------------------------
    # 4. MaximalExchangeVirtuals: most negative exchange eigenvalues
    # ---------------------------------------------------------
    target = MaximalExchangeVirtuals(n_orbitals=N_virt_target)
    virt_mev = generate_orbitals(target, occ_space, ham; solver=FullDiagonalization())
    @test size(virt_mev.ψ[1]) == (Nfull, N_virt_target)
    @test virt_mev.ψ[1]' * virt_mev.ψ[1] ≈ I
    @test norm(occ_space.ψ[1]' * virt_mev.ψ[1]) < 1e-6
    @test all(virt_mev.eigenvalues[1] .< 0)
    @test issorted(virt_mev.eigenvalues[1])
    virt_mev_lobpcg = generate_orbitals(target, occ_space, ham; solver=solver_lobpcg)
    @test isapprox(virt_mev_lobpcg.eigenvalues[1], virt_mev.eigenvalues[1], rtol=1e-3)
end
