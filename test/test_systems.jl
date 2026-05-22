@testmodule TestSystems begin
export setup_water_hf

using DFTK
using PseudoPotentialData
using LinearAlgebra

"""
    setup_water_hf(; Ecut=15, box_size=10.0)

Sets up a simple water molecule in a cubic box and performs a Hartree-Fock
calculation to return a converged `scfres` object for testing purposes.
"""
function setup_water_hf(; Ecut=15, box_size=10.0)
    pd_pbe_family = PseudoFamily("dojo.nc.sr.pbe.v0_5.stringent.upf")
    O = ElementPsp(:O, pd_pbe_family)
    H = ElementPsp(:H, pd_pbe_family)
    atoms = [O, H, H]
    
    a = box_size
    lattice = a * I(3)
    n_bands_converge = 8  # corresponds to 4 occupied and 4 virtual bands
    
    # H2O geometry in fractional coordinates
    positions = [
        [1/2, 1/2, 1/2],
        [1/2 + 1.45/a, 1/2 - 1.08/a, 1/2],
        [1/2 - 1.45/a, 1/2 - 1.08/a, 1/2]
    ]

    # --- Run PBE to get a good initial guess
    model_pbe = model_DFT(lattice, atoms, positions; functionals=PBE())
    basis = PlaneWaveBasis(model_pbe; Ecut=Ecut, kgrid=[1, 1, 1])
    scfres_pbe = self_consistent_field(
        basis; 
        is_converged=ScfConvergenceEnergy(1e-5), 
        nbandsalg=FixedBands(; n_bands_converge)
    )
    
    # --- Run HF using the PBE density and orbitals as a starting point
    model_hf = model_HF(
        lattice, 
        atoms, 
        positions; 
        exx_kernel=DFTK.Coulomb(DFTK.ProbeCharge())
    )
    basis_hf = PlaneWaveBasis(model_hf; Ecut=Ecut, kgrid=[1, 1, 1])
    
    scfres_hf = self_consistent_field(
        basis_hf;
        solver=DFTK.scf_damping_solver(damping=1.0),
        is_converged=ScfConvergenceEnergy(1e-7),
        ρ=scfres_pbe.ρ,
        ψ=scfres_pbe.ψ,
        occupation=scfres_pbe.occupation,
        diagtolalg=DFTK.AdaptiveDiagtol(; ratio_ρdiff=1e-4),
        exxalg=DFTK.AceExx(),
        nbandsalg=FixedBands(; n_bands_converge)
    )
    
    return scfres_hf
end

end
