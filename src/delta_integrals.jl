function _ifft_matrix(basis, kpt, ψ_mat)
    Nr = prod(basis.fft_size)
    N = size(ψ_mat, 2)
    ψ_real_flat = zeros(ComplexF64, Nr, N)
    for i in 1:N
        ψ_real_flat[:, i] = vec(DFTK.ifft(basis, kpt, ψ_mat[:, i]))
    end
    return ψ_real_flat
end

@doc raw"""
    compute_delta_integrals(basis, ψ_holes, ::Val{:HH})

Compute the 2-index Delta integrals (Hole-Hole) of the hole orbitals `ψ_holes` (per
k-point, Gamma-only) on the real space grid of `basis`.

# Mathematical Definition
These integrals represent the overlap between two hole orbitals evaluated on the
real-space grid:
```math
\Delta_{ij} = \int d^3r \, \psi_i^*(r) \psi_j(r)
            \approx \sum_{r} \psi_i^*(r) \psi_j(r) \, \Delta V
```
For perfectly orthonormal orbitals and an infinitely dense grid, this evaluates to the
Kronecker delta ``\delta_{ij}``. The grid-based numerical representation is returned.

# Cost
With ``N_\text{occ}`` hole orbitals:
- time: ``O(N_\text{occ} N_r \log N_r + N_r N_\text{occ}^2)``
- memory: ``N_r N_\text{occ}`` complex numbers for the orbitals on the grid
"""
function compute_delta_integrals(basis, ψ_holes, ::Val{:HH})
    # Transform holes to real space and flatten spatial dimensions
    ψ_holes_real_flat = _ifft_matrix(basis, basis.kpoints[1], ψ_holes[1])

    # DeltaIntegrals_ij = sum_r ψ_i^*(r) ψ_j(r) * dvol
    DeltaIntegralsHH = (ψ_holes_real_flat' * ψ_holes_real_flat) .* basis.dvol
    return DeltaIntegralsHH
end

@doc raw"""
    compute_delta_integrals(basis, ψ_particles, ψ_holes, ::Val{:PPHH})

Compute the 4-index Delta integrals (Particle-Particle-Hole-Hole) of the particle orbitals
`ψ_particles` and hole orbitals `ψ_holes` (per k-point, Gamma-only) on the real space grid
of `basis`.

# Mathematical Definition
These integrals represent the 4-orbital pointwise overlap between two particle orbitals
and two hole orbitals:
```math
\Delta_{abij} = \int d^3r \, \psi_a^*(r) \psi_b^*(r) \psi_i(r) \psi_j(r)
              \approx \sum_{r} \psi_a^*(r) \psi_b^*(r) \psi_i(r) \psi_j(r) \, \Delta V
```
These quantities naturally emerge when decomposing two-electron Coulomb integrals using
resolution of the identity or real-space vertex tensors. The returned tensor has the
dimensions `(Nvirt, Nvirt, Nocc, Nocc)`.

# Cost
With ``N_\text{virt}`` particle and ``N_\text{occ}`` hole orbitals:
- time: ``O(N_r N_\text{virt}^2 N_\text{occ}^2)`` for the contraction of the orbital pairs
- memory: ``N_r (N_\text{virt}^2 + N_\text{occ}^2)`` complex numbers for the orbital pairs
  on the grid. This is the practical limit, e.g. 500 virtual orbitals on a ``36^3`` grid
  need about 190 GB.
"""
function compute_delta_integrals(basis, ψ_particles, ψ_holes, ::Val{:PPHH})
    ψ_holes_real_flat = _ifft_matrix(basis, basis.kpoints[1], ψ_holes[1])
    ψ_particles_real_flat = _ifft_matrix(basis, basis.kpoints[1], ψ_particles[1])

    Nr = prod(basis.fft_size)
    Nocc = size(ψ_holes[1], 2)
    Nvirt = size(ψ_particles[1], 2)

    # We want DeltaIntegrals_abij = sum_r ψ_a^*(r) ψ_b^*(r) ψ_i(r) ψ_j(r) * dvol
    # For efficiency with BLAS, we form pairs:
    # V[r, (a,b)] = ψ_a(r) * ψ_b(r)
    # O[r, (i,j)] = ψ_i(r) * ψ_j(r)
    # Then V' * O does the complex conjugate on V!

    V_pairs = zeros(ComplexF64, Nr, Nvirt * Nvirt)
    idx = 1
    for b in 1:Nvirt, a in 1:Nvirt
        @. V_pairs[:, idx] = ψ_particles_real_flat[:, a] * ψ_particles_real_flat[:, b]
        idx += 1
    end

    O_pairs = zeros(ComplexF64, Nr, Nocc * Nocc)
    idx = 1
    for j in 1:Nocc, i in 1:Nocc
        @. O_pairs[:, idx] = ψ_holes_real_flat[:, i] * ψ_holes_real_flat[:, j]
        idx += 1
    end

    # DeltaIntegrals_abij = (V_pairs' * O_pairs) * dvol
    # Size: (Nvirt * Nvirt, Nocc * Nocc)
    DeltaIntegralsPPHH = (V_pairs' * O_pairs) .* basis.dvol

    # Reshape to 4D tensor (Nvirt, Nvirt, Nocc, Nocc)
    return reshape(DeltaIntegralsPPHH, Nvirt, Nvirt, Nocc, Nocc)
end
