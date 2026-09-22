using DFTK
using PsiTK
using PseudoPotentialData

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

    # PBE run in DFTK as initial guess for the HF solver
    # (no symmetries, like the HF model, so that both bases share the FFT grid)
    model = model_DFT(lattice, atoms, positions; functionals = PBE(), symmetries = false)
    basis = PlaneWaveBasis(model; Ecut = Ecut, kgrid = [1, 1, 1])
    println("run PBE")
    scfres_pbe = self_consistent_field(basis; is_converged = ScfConvergenceEnergy(1e-7))

    model = model_HF(lattice, atoms, positions; exx_kernel = Coulomb(ProbeCharge()))
    basis = PlaneWaveBasis(model; Ecut = Ecut, kgrid = [1, 1, 1])
    println("run HF")
    scfres_hf = self_consistent_field(
        basis;
        solver = ScfDampingSolver(),
        is_converged = ScfConvergenceEnergy(1e-7),
        ρ = scfres_pbe.ρ,
        ψ = scfres_pbe.ψ,
        occupation = scfres_pbe.occupation,
        maxiter = 100,
        diagtolalg = DFTK.AdaptiveDiagtol(; ratio_ρdiff = 5e-4),
        exxalg = DFTK.AceExx(),
    )

    # Occupied HF orbitals as the starting point
    occ_space = extract_occupied_space(OrbitalSpace(scfres_hf))

    # Generate a compressed virtual space (44 Density Specific Virtuals)
    println("Compute DSVs")
    target = DensitySpecificVirtuals(scfres_hf, occ_space; n_orbitals = 44)
    dsv_space = generate_orbitals(target, occ_space)

    # DSVs are not orthonormal and carry no Fock energies: canonicalize the active space
    # (diagonalize the Fock operator in the merged subspace) before it can be used
    println("Canonicalize Active Space")
    active_space = canonicalize_orbitals(merge_spaces(occ_space, dsv_space), scfres_hf.ham)

    # Compute the Coulomb Vertex for the Active Space
    println("Compute Coulomb Vertex")
    ΓmnG, G_vectors, kernel_fourier = compute_coulomb_vertex(active_space; callback = ShowProgress())
    vertex_alg = CoulombGramian()
    ΓmnF, coulomb_vertex_singular_vectors = compress_coulomb_vertex(ΓmnG, vertex_alg)

    # Dump to the specific correlation solver (Cc4s)
    println("prepare and dump Cc4s files")
    dump_cc4s_files(
        active_space, ΓmnF, G_vectors, kernel_fourier;
        coulomb_vertex_singular_vectors = coulomb_vertex_singular_vectors,
        folder = @__DIR__,
        force = true,
    )

    println("done")
end

main()
