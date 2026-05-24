using Pkg
Pkg.activate(".")
using PsiTK
using DFTK
using LinearAlgebra
using Random

include("test_systems.jl")
using .TestSystems

scfres = setup_water_hf(Ecut=15)
basis = scfres.basis
occ_space = compute_occupied_space(scfres)

# Canonical - FullDiag
target_canon_all = CanonicalVirtuals(:all)
solver_fd = FullDiagonalizationSolver()
virt_canon_fd = generate_orbitals(target_canon_all, scfres, occ_space, solver_fd)

println("FullDiag Canon all N: ", size(virt_canon_fd.ψ[1], 2))
println("FullDiag Canon eigenvalues[1:5]: ", virt_canon_fd.eigenvalues[1][1:5])

target_canon_5 = CanonicalVirtuals(5)
virt_canon_fd_5 = generate_orbitals(target_canon_5, scfres, occ_space, solver_fd)
println("FullDiag Canon 5 eigenvalues: ", virt_canon_fd_5.eigenvalues[1])

# Canonical - LOBPCG
solver_lobpcg = LOBPCGSolver(tol=1e-6, maxiter=500)
virt_canon_lobpcg = generate_orbitals(target_canon_5, scfres, occ_space, solver_lobpcg)
println("LOBPCG Canon 5 eigenvalues: ", virt_canon_lobpcg.eigenvalues[1])

# DSV
target_dsv = DensitySpecificVirtuals(5)
# Run once
virt_dsv_1 = generate_orbitals(target_dsv, scfres, occ_space, solver_lobpcg)
println("DSV 1 norm: ", norm(virt_dsv_1.eigenvalues[1]))
# Run twice to test reproducibility without seed
virt_dsv_2 = generate_orbitals(target_dsv, scfres, occ_space, solver_lobpcg)
println("DSV 2 norm: ", norm(virt_dsv_2.eigenvalues[1]))

