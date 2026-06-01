export compute_delta_integrals

function _ifft_matrix(basis, kpt, ψ_mat)
    N_grid = prod(basis.fft_size)
    N_bands = size(ψ_mat, 2)
    ψ_real_flat = zeros(ComplexF64, N_grid, N_bands)
    for i in 1:N_bands
        ψ_real_flat[:, i] = vec(DFTK.ifft(basis, kpt, ψ_mat[:, i]))
    end
    return ψ_real_flat
end

"""
    compute_delta_integrals(holes::OrbitalSpace, ::Val{:HH})

Compute the 2-index Delta integrals (Hole-Hole) on the real space grid.

# Mathematical Definition
These integrals represent the overlap between two hole orbitals evaluated on the real-space grid:
```math
\\Delta_{ij} = \\int d^3r \\, \\psi_i^*(r) \\psi_j(r) \\approx \\sum_{r} \\psi_i^*(r) \\psi_j(r) \\, \\Delta V
```
For perfectly orthonormal orbitals and an infinitely dense grid, this evaluates to the Kronecker delta ``\\delta_{ij}``. The grid-based numerical representation is returned.
"""
function compute_delta_integrals(holes::OrbitalSpace, ::Val{:HH})
    basis = holes.basis
    
    # Transform holes to real space and flatten spatial dimensions
    ψ_holes_real_flat = _ifft_matrix(basis, basis.kpoints[1], holes.ψ[1])
    
    # DeltaIntegrals_ij = sum_r ψ_i^*(r) ψ_j(r) * dvol
    DeltaIntegralsHH = (ψ_holes_real_flat' * ψ_holes_real_flat) .* basis.dvol
    return DeltaIntegralsHH
end

"""
    compute_delta_integrals(particles::OrbitalSpace, holes::OrbitalSpace, ::Val{:PPHH})

Compute the 4-index Delta integrals (Particle-Particle-Hole-Hole) on the real space grid.

# Mathematical Definition
These integrals represent the 4-orbital pointwise overlap between two particle orbitals and two hole orbitals:
```math
\\Delta_{abij} = \\int d^3r \\, \\psi_a^*(r) \\psi_b^*(r) \\psi_i(r) \\psi_j(r) \\approx \\sum_{r} \\psi_a^*(r) \\psi_b^*(r) \\psi_i(r) \\psi_j(r) \\, \\Delta V
```
These quantities naturally emerge when decomposing two-electron Coulomb integrals using resolution of the identity or real-space vertex tensors. 
The returned tensor has the dimensions ``(N_{\\text{particles}}, N_{\\text{particles}}, N_{\\text{holes}}, N_{\\text{holes}})``.
"""
function compute_delta_integrals(particles::OrbitalSpace, holes::OrbitalSpace, ::Val{:PPHH})
    basis = holes.basis
    
    ψ_holes_real_flat = _ifft_matrix(basis, basis.kpoints[1], holes.ψ[1])
    ψ_particles_real_flat = _ifft_matrix(basis, basis.kpoints[1], particles.ψ[1])
    
    N_grid = prod(basis.fft_size)
    N_holes = size(holes.ψ[1], 2)
    N_particles = size(particles.ψ[1], 2)
    
    # We want DeltaIntegrals_abij = sum_r ψ_a^*(r) ψ_b^*(r) ψ_i(r) ψ_j(r) * dvol
    # For efficiency with BLAS, we form pairs:
    # V[r, (a,b)] = ψ_a(r) * ψ_b(r)
    # O[r, (i,j)] = ψ_i(r) * ψ_j(r)
    # Then V' * O does the complex conjugate on V!
    
    V_pairs = zeros(ComplexF64, N_grid, N_particles * N_particles)
    idx = 1
    for b = 1:N_particles, a = 1:N_particles
        @. V_pairs[:, idx] = ψ_particles_real_flat[:, a] * ψ_particles_real_flat[:, b]
        idx += 1
    end
    
    O_pairs = zeros(ComplexF64, N_grid, N_holes * N_holes)
    idx = 1
    for j = 1:N_holes, i = 1:N_holes
        @. O_pairs[:, idx] = ψ_holes_real_flat[:, i] * ψ_holes_real_flat[:, j]
        idx += 1
    end
    
    # DeltaIntegrals_abij = (V_pairs' * O_pairs) * dvol
    # Size: (N_particles * N_particles, N_holes * N_holes)
    DeltaIntegralsPPHH = (V_pairs' * O_pairs) .* basis.dvol
    
    # Reshape to 4D tensor (N_particles, N_particles, N_holes, N_holes)
    return reshape(DeltaIntegralsPPHH, N_particles, N_particles, N_holes, N_holes)
end
