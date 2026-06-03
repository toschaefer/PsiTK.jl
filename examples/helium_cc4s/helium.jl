using DFTK
using PsiTK
using PseudoPotentialData
using LinearAlgebra

function main()
    pd_pbe_family = PseudoFamily("dojo.nc.sr.pbe.v0_5.stringent.upf")

    He = ElementPsp(:He, pd_pbe_family)
    atoms = [He]
    box_length = 10.0
    lattice = [
        [box_length 0.000000 0.00000];
        [0.00000 box_length-0.1 0.00000];
        [0.00000 0.000000 box_length-0.2]
    ]
    positions = [[0.500000, 0.500000, 0.500000]]
    Ecut = 12

    # Run standard DFTK calculations
    model = model_PBE(lattice, atoms, positions)
    basis = PlaneWaveBasis(model; Ecut = Ecut, kgrid = [1, 1, 1])
    println("run PBE")
    scfres_pbe = self_consistent_field(basis; is_converged = ScfConvergenceEnergy(1e-7))

    # use the PBE solution as initial guess for the HF solver
    model = model_HF(lattice, atoms, positions; exx_kernel = Coulomb(ProbeCharge()))
    basis = PlaneWaveBasis(model; Ecut = Ecut, kgrid = [1, 1, 1])
    println("run HF")
    scfres_hf = self_consistent_field(
        basis;
        solver = DFTK.scf_damping_solver(damping = 1.0),
        is_converged = ScfConvergenceEnergy(1e-7),
        tol = 1e-5,
        ρ = scfres_pbe.ρ,
        ψ = scfres_pbe.ψ,
        occupation = scfres_pbe.occupation,
        maxiter = 100,
        diagtolalg = DFTK.AdaptiveDiagtol(; ratio_ρdiff = 5e-4),
        exxalg = AceExx(),
    )

    # Extract the base OrbitalSpace from the HF result
    scf_space = OrbitalSpace(scfres_hf)
    occ_space, _ = split_space_occupied_virtual(scf_space; threshold = 1e-6)

    # Generate a generalized set of orbitals (e.g., 44 DSVs)
    println("Compute DSVs")
    target = DensitySpecificVirtuals(scfres_hf, occ_space; n_orbitals=44)
    dsv_space = generate_orbitals(target, occ_space)
    println(dsv_space.eigenvalues)

    # Merge spaces to create the Active Space using lazy views to save memory
    println("Merge spaces")
    active_space = merge_spaces(occ_space, dsv_space)

    # Canonicalize the active space (diagonalize the Fock Hamiltonian)
    println("Canonicalize Active Space")
    active_space = diagonalize_orbitals(active_space, scfres_hf.ham)

    # Compute the Coulomb Vertex for the Active Space
    println("Compute Coulomb Vertex")
    ΓmnG, G_vectors, kernel_fourier = compute_coulomb_vertex(active_space)
    vertex_alg = CoulombGramian()
    ΓmnF, coulomb_vertex_singular_vectors = compress_coulomb_vertex(ΓmnG, vertex_alg)

    # Dump to the specific correlation solver (Cc4s)
    println("prepare and dump Cc4s files")
    dump_cc4s_files(active_space, ΓmnF, G_vectors, kernel_fourier; coulomb_vertex_singular_vectors = coulomb_vertex_singular_vectors, folder = @__DIR__, force = true)

    println("done")
end

main()
