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

    # Canonicalizing a merged space distributes the occupations by aufbau, independent of
    # the merge order, and preserves the occupied subspace of the converged HF (no warning)
    active = @test_nowarn canonicalize_orbitals(
        merge_spaces(occ_space, virt_dsv),
        ham,
    )
    active_rev = canonicalize_orbitals(merge_spaces(virt_dsv, occ_space), ham)
    ψ_active = active.ψ[1]
    f_active = active.occupation[1]
    @test ψ_active' * ψ_active ≈ I
    @test issorted(active.eigenvalues[1])
    @test f_active == sort(vcat(occ_space.occupation[1], virt_dsv.occupation[1]); rev=true)
    @test active_rev.occupation[1] == f_active
    @test active_rev.eigenvalues[1] ≈ active.eigenvalues[1]
    ψ_occ = ψ_active[:, f_active .> 0]
    ψ_occ_hf = occ_space.ψ[1]
    @test norm(ψ_occ - ψ_occ_hf * (ψ_occ_hf' * ψ_occ)) < 1e-4

    # The Fermi level is moved into the gap of the canonical orbitals
    ε_active = active.eigenvalues[1]
    @test maximum(ε_active[f_active .> 0]) < active.εF < minimum(ε_active[f_active .== 0])

    # Splitting does not rely on the occupied orbitals coming first
    occ_split, virt_split = split_occupied_virtual(merge_spaces(virt_dsv, occ_space))
    @test occ_split.ψ[1] == occ_space.ψ[1]
    @test virt_split.ψ[1] == virt_dsv.ψ[1]

    # Mixing HOMO and LUMO by a rotation gives a determinant the Fock operator does not
    # preserve: the aufbau occupation differs from the input one, which is warned about
    X_mixed = copy(ψ_active)
    rotation = [cos(0.1) -sin(0.1); sin(0.1) cos(0.1)]
    X_mixed[:, [Nocc, Nocc + 1]] = X_mixed[:, [Nocc, Nocc + 1]] * rotation
    mixed = OrbitalSpace(
        basis,
        [X_mixed],
        active.eigenvalues,
        active.occupation,
        active.εF,
        true,
    )
    @test_logs (:warn, r"different occupation") canonicalize_orbitals(mixed, ham)

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
