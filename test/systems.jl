@testmodule TestSystems begin
    export setup_water_hf, setup_diamond_hf, setup_he_solid_hf, setup_h_chain_hf

    using DFTK
    using PseudoPotentialData
    using LinearAlgebra

    """
        setup_water_hf(; Ecut=15, box_size=10.0)

    Sets up a simple water molecule in a cubic box and performs a Hartree-Fock
    calculation to return a converged `scfres` object for testing purposes.
    """
    function setup_water_hf(; Ecut = 15, box_size = 10.0, n_bands_converge = nothing)
        pd_pbe_family = PseudoFamily("dojo.nc.sr.pbe.v0_5.stringent.upf")
        O = ElementPsp(:O, pd_pbe_family)
        H = ElementPsp(:H, pd_pbe_family)
        atoms = [O, H, H]

        a = box_size
        lattice = a * I(3)
        kwargs_nbands = isnothing(n_bands_converge) ? (;) : (; nbandsalg=FixedBands(; n_bands_converge))

        # H2O geometry in fractional coordinates
        positions = [
            [1/2, 1/2, 1/2],
            [1/2 + 1.45/a, 1/2 - 1.08/a, 1/2],
            [1/2 - 1.45/a, 1/2 - 1.08/a, 1/2],
        ]

        # --- Run PBE to get a good initial guess
        model_pbe = model_DFT(lattice, atoms, positions; functionals = PBE())
        basis = PlaneWaveBasis(model_pbe; Ecut = Ecut, kgrid = [1, 1, 1])
        scfres_pbe = self_consistent_field(
            basis;
            is_converged = ScfConvergenceEnergy(1e-6),
            kwargs_nbands...
        )

        # --- Run HF using the PBE density and orbitals as a starting point
        model_hf = model_HF(
            lattice,
            atoms,
            positions;
            exx_kernel = DFTK.Coulomb(DFTK.ProbeCharge()),
        )
        basis_hf = PlaneWaveBasis(model_hf; Ecut = Ecut, kgrid = [1, 1, 1])

        scfres_hf = self_consistent_field(
            basis_hf;
            solver = DFTK.scf_damping_solver(damping = 1.0),
            is_converged = ScfConvergenceEnergy(1e-8),
            ρ = scfres_pbe.ρ,
            ψ = scfres_pbe.ψ,
            occupation = scfres_pbe.occupation,
            diagtolalg = DFTK.AdaptiveDiagtol(; ratio_ρdiff = 1e-4),
            exxalg = DFTK.AceExx(),
            kwargs_nbands...
        )

        return scfres_hf
    end

    """
        setup_diamond_hf(; Ecut=20)

    Sets up a conventional diamond cell (8 atoms) at the gamma point only and performs 
    a Hartree-Fock calculation to return a converged `scfres` object.
    """
    function setup_diamond_hf(; Ecut = 20, n_bands_converge = nothing)
        a = 6.74  # Bohr
        lattice = a * I(3)
        pd_pbe_family = PseudoFamily("dojo.nc.sr.pbe.v0_5.stringent.upf")
        C = ElementPsp(:C, pd_pbe_family)
        atoms = fill(C, 8)
        positions = [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
            [0.25, 0.25, 0.25],
            [0.75, 0.75, 0.25],
            [0.75, 0.25, 0.75],
            [0.25, 0.75, 0.75]
        ]
        
        kwargs_nbands = isnothing(n_bands_converge) ? (;) : (; nbandsalg=FixedBands(; n_bands_converge))

        model_pbe = model_DFT(lattice, atoms, positions; functionals = PBE())
        basis = PlaneWaveBasis(model_pbe; Ecut = Ecut, kgrid = [1, 1, 1])
        scfres_pbe = self_consistent_field(
            basis;
            is_converged = ScfConvergenceEnergy(1e-6),
            kwargs_nbands...
        )

        model_hf = model_HF(
            lattice,
            atoms,
            positions;
            exx_kernel = DFTK.Coulomb(DFTK.ProbeCharge()),
        )
        basis_hf = PlaneWaveBasis(model_hf; Ecut = Ecut, kgrid = [1, 1, 1])
        scfres_hf = self_consistent_field(
            basis_hf;
            solver = DFTK.scf_damping_solver(damping = 1.0),
            is_converged = ScfConvergenceEnergy(1e-8),
            ρ = scfres_pbe.ρ,
            ψ = scfres_pbe.ψ,
            occupation = scfres_pbe.occupation,
            diagtolalg = DFTK.AdaptiveDiagtol(; ratio_ρdiff = 1e-4),
            exxalg = DFTK.AceExx(),
            kwargs_nbands...
        )
        
        return scfres_hf
    end

    """
        setup_he_solid_hf(; Ecut=20)

    Sets up a solid Helium supercell (6 atoms, HCP) at the gamma point only and performs 
    a Hartree-Fock calculation to return a converged `scfres` object.
    """
    function setup_he_solid_hf(; Ecut = 20, n_bands_converge = nothing)
        a = 6.75  # Bohr
        c = 11.02 # Bohr
        lattice = [
            3a -a/2 0;
            0 a*sqrt(3)/2 0;
            0 0 c
        ]
        pd_pbe_family = PseudoFamily("dojo.nc.sr.pbe.v0_5.stringent.upf")
        He = ElementPsp(:He, pd_pbe_family)
        atoms = fill(He, 6)
        positions = [
            [0.0, 0.0, 0.0],
            [1/9, 2/3, 1/2],
            [1/3, 0.0, 0.0],
            [4/9, 2/3, 1/2],
            [2/3, 0.0, 0.0],
            [7/9, 2/3, 1/2]
        ]

        kwargs_nbands = isnothing(n_bands_converge) ? (;) : (; nbandsalg=FixedBands(; n_bands_converge))

        model_pbe = model_DFT(lattice, atoms, positions; functionals = PBE())
        basis = PlaneWaveBasis(model_pbe; Ecut = Ecut, kgrid = [1, 1, 1])
        scfres_pbe = self_consistent_field(
            basis;
            is_converged = ScfConvergenceEnergy(1e-6),
            kwargs_nbands...
        )

        model_hf = model_HF(
            lattice,
            atoms,
            positions;
            exx_kernel = DFTK.Coulomb(DFTK.ProbeCharge()),
        )
        basis_hf = PlaneWaveBasis(model_hf; Ecut = Ecut, kgrid = [1, 1, 1])
        scfres_hf = self_consistent_field(
            basis_hf;
            solver = DFTK.scf_damping_solver(damping = 1.0),
            is_converged = ScfConvergenceEnergy(1e-8),
            ρ = scfres_pbe.ρ,
            ψ = scfres_pbe.ψ,
            occupation = scfres_pbe.occupation,
            diagtolalg = DFTK.AdaptiveDiagtol(; ratio_ρdiff = 1e-4),
            exxalg = DFTK.AceExx(),
            kwargs_nbands...
        )
        
        return scfres_hf
    end

    """
        setup_h_chain_hf(; Ecut=15, n_bands_converge=nothing)

    Sets up a Peierls distorted hydrogen chain with 4 atoms in the cell.
    The y and z directions are slightly distorted to break symmetry.
    """
    function setup_h_chain_hf(; Ecut = 15, n_bands_converge = nothing)
        L_x = 6.0
        lattice = [
            L_x 0 0;
            0 1.05*L_x 0;
            0 0 1.10*L_x
        ]
        pd_pbe_family = PseudoFamily("dojo.nc.sr.pbe.v0_5.stringent.upf")
        H = ElementPsp(:H, pd_pbe_family)
        atoms = fill(H, 4)
        
        # Peierls distortion along x
        positions = [
            [0.0, 0.5, 0.5],
            [0.2, 0.5, 0.5],
            [0.5, 0.5, 0.5],
            [0.7, 0.5, 0.5]
        ]

        kwargs_nbands = isnothing(n_bands_converge) ? (;) : (; nbandsalg=FixedBands(; n_bands_converge))

        model_pbe = model_DFT(lattice, atoms, positions; functionals = PBE())
        basis = PlaneWaveBasis(model_pbe; Ecut = Ecut, kgrid = [1, 1, 1])
        scfres_pbe = self_consistent_field(
            basis;
            is_converged = ScfConvergenceEnergy(1e-6),
            kwargs_nbands...
        )

        model_hf = model_HF(
            lattice,
            atoms,
            positions;
            exx_kernel = DFTK.Coulomb(DFTK.ProbeCharge()),
        )
        basis_hf = PlaneWaveBasis(model_hf; Ecut = Ecut, kgrid = [1, 1, 1])
        scfres_hf = self_consistent_field(
            basis_hf;
            solver = DFTK.scf_damping_solver(damping = 1.0),
            is_converged = ScfConvergenceEnergy(1e-8),
            ρ = scfres_pbe.ρ,
            ψ = scfres_pbe.ψ,
            occupation = scfres_pbe.occupation,
            diagtolalg = DFTK.AdaptiveDiagtol(; ratio_ρdiff = 1e-4),
            exxalg = DFTK.AceExx(),
            kwargs_nbands...
        )
        
        return scfres_hf
    end

end
