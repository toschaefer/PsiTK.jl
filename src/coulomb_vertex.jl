@doc raw"""
    compute_overlap_densities(basis, ψ_bra, ψ_ket; Ecut_ratio=1.0, callback=identity)
    compute_overlap_densities(basis, ψ; kwargs...)
    compute_overlap_densities(bra_space::OrbitalSpace, ket_space::OrbitalSpace; kwargs...)
    compute_overlap_densities(space::OrbitalSpace; kwargs...)

Compute the overlap densities in reciprocal space
```math
ρ_{mn \bm G} = \frac{1}{\sqrt{Ω}} \int_Ω \; \psi_{m}(\bm r)^∗ \psi_{n}(\bm r)
               \; e^{-i\bm r \bm G}  \; d^3 r
```
i.e. the coefficients of ``\psi_m^∗ \psi_n`` in the orthonormal plane waves
``e^{i\bm r \bm G}/\sqrt{Ω}`` (DFTK's `fft` convention), for all orbitals in `ψ_bra` and
`ψ_ket`. To restrict the orbitals, pass a subset (e.g. via [`select_orbitals`](@ref)).

# Arguments
- `basis`: the `PlaneWaveBasis` of the orbitals
- `ψ_bra`: the bra orbitals per k-point (e.g. occupied orbitals)
- `ψ_ket`: the ket orbitals per k-point (e.g. virtual orbitals)
- `ψ`: orbitals used as both bra and ket. Only the upper triangle ``m ≤ n`` is computed,
  the rest follows from ``ρ_{nm,-\bm G} = ρ_{mn \bm G}^∗``. This shortcut is taken
  whenever `ψ_bra === ψ_ket`.
- `bra_space`, `ket_space`, `space`: [`OrbitalSpace`](@ref)s, a shorthand for passing
  their `basis` and `ψ`
- `Ecut_ratio`: ratio of the plane-wave cutoff (in energy) for the densities relative to
  the orbital cutoff `basis.Ecut` (default: 1.0). Values up to `supersampling^2` of the
  basis (4 for DFTK's default) are allowed since the FFT grid holds products of orbitals
  exactly up to there; `Ecut_ratio=4` therefore yields the exact overlap densities.
- `callback`: called after each orbital pair with `(; step, total_steps)`,
  e.g. `callback=ShowProgress()` for a progress bar (default: no output)

# Returns
A tuple `(ρmnG, G_vectors)`:
- `ρmnG`: the overlap densities as a tensor of shape `(nk, n_bra, nk, n_ket, nG)`.
- `G_vectors`: the corresponding plane-wave vectors.
"""
function compute_overlap_densities(
    basis,
    ψ_bra,
    ψ_ket;
    Ecut_ratio=1.0,
    callback=identity,
)
    all(kpt -> iszero(kpt.coordinate), basis.kpoints) ||
        error("Overlap densities are only implemented for Gamma-point calculations.")
    G_indices = _G_indices_within_cutoff(basis, Ecut_ratio)
    ρmnG = _compute_overlap_densities(basis, ψ_bra, ψ_ket; G_indices, callback)
    return ρmnG, G_vectors(basis)[G_indices]
end
function compute_overlap_densities(basis, ψ; kwargs...)
    return compute_overlap_densities(basis, ψ, ψ; kwargs...)
end
function compute_overlap_densities(
    bra_space::OrbitalSpace,
    ket_space::OrbitalSpace;
    kwargs...,
)
    return compute_overlap_densities(bra_space.basis, bra_space.ψ, ket_space.ψ; kwargs...)
end
function compute_overlap_densities(space::OrbitalSpace; kwargs...)
    return compute_overlap_densities(space.basis, space.ψ; kwargs...)
end

# Linear indices into the full FFT cube G_vectors(basis) of the G vectors with
# |G|^2/2 <= Ecut * Ecut_ratio. The cube covers -cld(N-1,2):fld(N-1,2) per axis; requiring
# the sphere to fit into the smaller side (N-1)÷2 keeps the selection closed under G -> -G.
function _G_indices_within_cutoff(basis, Ecut_ratio)
    Ecut_reduced = basis.Ecut * Ecut_ratio
    # The sphere of radius Gmax extends to |n_i| <= Gmax |a_i| / 2π in integer coordinates
    Gmax = sqrt(2 * Ecut_reduced)
    for (i, a_i) in enumerate(eachcol(basis.model.lattice))
        Gmax * norm(a_i) / 2π <= (basis.fft_size[i] - 1) ÷ 2 ||
            error("Ecut_ratio=$Ecut_ratio exceeds the FFT grid of the basis " *
                  "(at most supersampling^2, i.e. 4 for DFTK's default).")
    end
    recip_lattice = basis.model.recip_lattice
    # vec: linear indices into the cube (findall on the 3D array would give CartesianIndex)
    return findall(vec(G_vectors(basis))) do G
        sum(abs2, recip_lattice * G) / 2 <= Ecut_reduced
    end
end

@doc raw"""
    compute_coulomb_vertex(
        basis,
        ψ_bra,
        ψ_ket;
        interaction_kernel=DFTK.ProbeCharge(DFTK.BareCoulomb()),
        Ecut_ratio=2/3,
        callback=identity
    )
    compute_coulomb_vertex(basis, ψ; kwargs...)
    compute_coulomb_vertex(bra_space::OrbitalSpace, ket_space::OrbitalSpace; kwargs...)
    compute_coulomb_vertex(space::OrbitalSpace; kwargs...)

Compute the Coulomb vertex
```math
Γ_{mn \bm G} = \sqrt{v(\bm G)} \; ρ_{mn \bm G}
```
where $ρ_{mn \bm G}$ are the overlap densities (see [`compute_overlap_densities`](@ref))
and $v(\bm G)$ is the interaction kernel, e.g. the Coulomb potential
```math
v(\bm G) = \frac{4π}{\bm G^2}
```
for all orbitals in `ψ_bra` and `ψ_ket`. To restrict the orbitals, pass a subset (e.g. via
[`select_orbitals`](@ref)).

# Arguments
- `basis`: the `PlaneWaveBasis` of the orbitals
- `ψ_bra`: the bra orbitals per k-point (e.g. occupied orbitals)
- `ψ_ket`: the ket orbitals per k-point (e.g. virtual orbitals)
- `ψ`: orbitals used as both bra and ket, exploiting the symmetry
  ``Γ_{nm,-\bm G} = Γ_{mn \bm G}^∗`` (see [`compute_overlap_densities`](@ref))
- `bra_space`, `ket_space`, `space`: [`OrbitalSpace`](@ref)s, a shorthand for passing
  their `basis` and `ψ`
- `interaction_kernel`: the DFTK `InteractionKernel` (default: bare Coulomb with the
  probe-charge singularity treatment, `ProbeCharge(BareCoulomb())`)
- `Ecut_ratio`: cutoff ratio for the vertex (default: 2/3), see
  [`compute_overlap_densities`](@ref)
- `callback`: called after each orbital pair with `(; step, total_steps)`,
  e.g. `callback=ShowProgress()` for a progress bar (default: no output)

# Returns
A [`DensityFitting`](@ref) holding the uncompressed vertex `Γ` (shape
`(nk, n_bra, nk, n_ket, nG)`), its `G_vectors` and the `kernel_fourier`.
"""
function compute_coulomb_vertex(
    basis,
    ψ_bra,
    ψ_ket;
    interaction_kernel=DFTK.ProbeCharge(DFTK.BareCoulomb()),
    Ecut_ratio=2/3,
    callback=identity,
)
    ρmnG, G_vectors = compute_overlap_densities(basis, ψ_bra, ψ_ket; Ecut_ratio, callback)

    # Kernel on the full FFT cube for momentum transfer q = 0 (Gamma-only)
    G_indices = _G_indices_within_cutoff(basis, Ecut_ratio)
    q = zero(DFTK.Vec3{eltype(basis)})
    kernel_cube = DFTK.eval_kernel_fourier(interaction_kernel, basis, q)
    kernel_fourier = kernel_cube[G_indices]

    # Γ = √v ⊙ ρ along the G axis; ρmnG is not needed anymore, so scale in place
    ΓmnG = ρmnG
    ΓmnG .*= reshape(sqrt.(kernel_fourier), 1, 1, 1, 1, :)

    return DensityFitting(ΓmnG, G_vectors, kernel_fourier, nothing)
end
function compute_coulomb_vertex(basis, ψ; kwargs...)
    return compute_coulomb_vertex(basis, ψ, ψ; kwargs...)
end
function compute_coulomb_vertex(
    bra_space::OrbitalSpace,
    ket_space::OrbitalSpace;
    kwargs...,
)
    return compute_coulomb_vertex(bra_space.basis, bra_space.ψ, ket_space.ψ; kwargs...)
end
function compute_coulomb_vertex(space::OrbitalSpace; kwargs...)
    return compute_coulomb_vertex(space.basis, space.ψ; kwargs...)
end

# This function initially based on code of the experimental "cc4s" branch in DFTK
# written by Michael Herbst
function _compute_overlap_densities(
    basis,
    ψ_bra::AbstractVector{<:AbstractArray{T}},
    ψ_ket::AbstractVector{<:AbstractArray{T}};
    G_indices=eachindex(G_vectors(basis)),
    callback=identity,
) where {T}
    n_kpt = length(basis.kpoints)
    n_bands_bra = size(ψ_bra[1], 2)
    n_bands_ket = size(ψ_ket[1], 2)

    # === Create index to map each stored G to -G on the full FFT cube ===
    Gs = G_vectors(basis)
    G_to_idx = Dict(Gs[i] => i for i in eachindex(Gs))
    idx_minus_G = [G_to_idx[-Gs[i]] for i in G_indices]

    # allocate overlap densities
    ρmnG = zeros(complex(T), n_kpt, n_bands_bra, n_kpt, n_bands_ket, length(G_indices))

    is_symmetric = (ψ_bra === ψ_ket)
    if is_symmetric
        # only upper triangle of ρmnG
        total_steps = (n_bands_bra * (n_bands_bra + 1) ÷ 2) * n_kpt^2
    else
        total_steps = n_bands_bra * n_bands_ket * n_kpt^2
    end
    step = 0

    # TODO:
    # Idea is to make some outer loop over the m-slices
    # for m_slice in m_slices
    #     precalculate ψmk_real for all m in this slice
    # end

    # === Calculate overlap densities ρmnG ===
    @views for (ikn, kptn) in enumerate(basis.kpoints), n in 1:n_bands_ket
        # Prepare ψnk(r)
        ψnk_real = ifft(basis, kptn, ψ_ket[ikn][:, n])

        for (ikm, kptm) in enumerate(basis.kpoints)
            for m in 1:n_bands_bra
                # Compute upper triangle only (m <= n) if spaces are symmetric
                # The lower triangle is filled via Hermitian conjugation below.
                if is_symmetric && m > n
                    continue
                end

                # Prepare ψmk(r)
                # TODO: pre-calculate some of them (not all, the virtual space can be large)
                ψmk_real = ifft(basis, kptm, ψ_bra[ikm][:, m])

                # Overlap density ρ_mn(r) = ψm*(r)ψn(r), FFT'd on the full cube
                overlap_density = fft(basis, conj.(ψmk_real) .* ψnk_real)

                # store entry of the overlap densities
                ρmnG[ikm, m, ikn, n, :] .= overlap_density[G_indices]

                # Fill lower triangle via ρmn(-G) = conjg(ρnmG)
                if is_symmetric && m != n
                    ρmnG[ikn, n, ikm, m, :] .= conj.(overlap_density[idx_minus_G])
                end

                step += 1
                callback((; step, total_steps))
            end
        end
    end
    return ρmnG
end

@doc raw"""
    DensityFitting

Density-fitting (resolution-of-identity) factorization of the electron repulsion integrals
```math
(pr|qs) \propto \sum_F Γ_{rpF}^∗ \, Γ_{qsF}
```
i.e. the Coulomb vertex $Γ$ together with the auxiliary basis $F$ it is expressed in.

# Fields
- `Γ`: the Coulomb vertex tensor of shape `(nk, n_bands, nk, n_bands, NF)`. The auxiliary
  index runs over the plane waves `G_vectors` (uncompressed) or over the compressed index
- `G_vectors`: the plane waves of the uncompressed vertex
- `kernel_fourier`: the interaction kernel evaluated at `G_vectors`
- `singular_vectors`: the transformation of shape `(NG, NF)` from the plane waves to the
  compressed auxiliary index, or `nothing` if uncompressed

See [`compute_coulomb_vertex`](@ref) and [`compress_coulomb_vertex`](@ref).
"""
struct DensityFitting{
    TΓ<:AbstractArray,
    TG<:AbstractVector,
    TV<:AbstractVector,
    TU<:Union{Nothing,AbstractMatrix},
}
    Γ::TΓ
    G_vectors::TG
    kernel_fourier::TV
    singular_vectors::TU
end

@doc raw"""
    compress_coulomb_vertex(fitting::DensityFitting, strategy)

Compress the Coulomb vertex along its auxiliary axis into a smaller auxiliary index $F$,
```math
Γ_{mn F} = \sum_{G} Γ_{mn G} \, U_{G F}
```
where the columns of the transformation $U$ span the dominant subspace of the Coulomb
Gramian $\Gamma^\dagger \Gamma$. How $U$ is determined depends on `strategy`:
- [`CoulombGramian`](@ref): exact diagonalization of the Gramian
- [`AdaptiveRandomizedSVD`](@ref): randomized range finder followed by diagonalization

# Arguments
- `fitting`: the [`DensityFitting`](@ref) from [`compute_coulomb_vertex`](@ref); compressing
  an already compressed fitting accumulates the transformations
- `strategy`: the compression strategy, carrying its own threshold

# Returns
A [`DensityFitting`](@ref) with the compressed `Γ` of shape `(nk, n_bands, nk, n_bands, NF)`
and the accumulated `singular_vectors` of shape `(NG, NF)`.
"""
function compress_coulomb_vertex(fitting::DensityFitting, strategy)
    ΓmnF, U = _compress_coulomb_vertex(fitting.Γ, strategy)
    U_prev = fitting.singular_vectors
    singular_vectors = isnothing(U_prev) ? U : U_prev * U
    return DensityFitting(ΓmnF, fitting.G_vectors, fitting.kernel_fourier, singular_vectors)
end

"""
    _compress_coulomb_vertex(Γ, strategy) -> (Γ_F, U)

Extension point for the compression strategies of [`compress_coulomb_vertex`](@ref):
compress the vertex `Γ` of shape `(nk, n_bands, nk, n_bands, NG)` along its auxiliary axis
to `Γ_F = Γ U` of shape `(nk, n_bands, nk, n_bands, NF)`, with the transformation `U` of
shape `(NG, NF)`.
"""
function _compress_coulomb_vertex end

@doc raw"""
    CoulombGramian(; thresh=1e-6)

Strategy for [`compress_coulomb_vertex`](@ref) through the largest eigenvalues of the
Coulomb Gramian
```math
H = - \Gamma^\dagger \Gamma = U \Lambda U^\dagger
```
The compressed $\Gamma$ is then obtained via $\Gamma_\text{compressed} = \Gamma U$,
where the columns of $U$ are restricted such that $|\lambda| >$ `thresh`.
"""
Base.@kwdef struct CoulombGramian
    thresh::Float64 = 1e-6
end
function _compress_coulomb_vertex(
    ΓmnG::AbstractArray{T,5},
    strategy::CoulombGramian,
) where {T}
    thresh = strategy.thresh
    Γmat = reshape(ΓmnG, prod(size(ΓmnG)[1:4]), size(ΓmnG, 5))
    Npp, NG = size(Γmat)

    H = -Hermitian(Γmat' * Γmat)         # Gramian in full PW basis
    λ, U = eigen(H)                      # diagonalize
    NF = findlast(s -> abs(s) > thresh, λ)  # truncate based on thresh
    if isnothing(NF)
        return ΓmnG, LinearAlgebra.I(size(Γmat, 2))
    else
        ΓmnF = Γmat * U[:, 1:NF]           # rotate
        return reshape(ΓmnF, size(ΓmnG)[1:4]..., NF), U[:, 1:NF]
    end
end

@doc raw"""
    AdaptiveRandomizedSVD(; thresh=1e-6, n_test_vectors=10)

Strategy for [`compress_coulomb_vertex`](@ref) via an adaptive randomized SVD.

The algorithm approximates the range of the row space of $\Gamma$ (the orbital indices
are considered as superindex) through a thin basis Q, such that
```math
\Gamma \approx \Gamma Q Q^\dagger
```
where $\Gamma$ is a $N_{pp} \times N_G$ and $Q$ a $N_G \times N_F$ matrix.
This is done through a stochastic Q and a diagonalization of
```math
H = -\tilde \Gamma^\dagger \tilde \Gamma = U \Lambda U^\dagger
```
where $\tilde \Gamma = \Gamma Q$.
The compressed $\Gamma$ is then obtained via $\Gamma_\text{compressed} = \tilde \Gamma U$,
the effective transformation matrix being $Q U$.

The dimension $N_F$ is found by a preceding adaptive range finder.
This finder iteratively increases the columns of Q (i.e. $N_F$) in steps of
$2\sqrt{N_{pp}}$ and stops when the error for each of `n_test_vectors` stochastic test
vectors $\omega_i$
```math
\varepsilon_i =  \Vert (1 - QQ^\dagger)\Gamma^\dagger \omega_i \Vert
```
is smaller than $\sqrt{\text{thresh}}/2$. With $r$ test vectors this estimator bounds the
true projection error with probability $1 - 10^{-r}$
[Halko, Martinsson, Tropp, SIAM Rev. **53**, 217 (2011), Lemma 4.1]; a single test vector
would stop the finder too early in a small fraction of runs.

TODO: The entire algorithm could be improved by techniques proposed in the following paper:
Fast and accurate randomized algorithms for low-rank tensor decompositions, L. Ma,
E. Solomonik (https://proceedings.neurips.cc/paper_files/paper/2021/hash/cbef46321026d8404bc3216d4774c8a9-Abstract.html)
"""
Base.@kwdef struct AdaptiveRandomizedSVD
    thresh::Float64 = 1e-6
    n_test_vectors::Int = 10
end
function _compress_coulomb_vertex(
    ΓmnG::AbstractArray{T,5},
    strategy::AdaptiveRandomizedSVD,
) where {T}
    thresh = strategy.thresh
    Γmat = reshape(ΓmnG, prod(size(ΓmnG)[1:4]), size(ΓmnG, 5))
    Npp, NG = size(Γmat)

    # === Adaptive Range Finder for NF ===

    # Store blocks of the stochastic guess basis
    Q_blocks = Matrix{T}[]

    # Step size for increasing the basis = 2*(√Npp)
    column_block_size = round(Int, 2 * Npp^0.5)

    # Stochastic test vectors for error estimation
    Ω_test = randn(T, Npp, strategy.n_test_vectors)

    # target error a little smaller than √thresh
    target_error = sqrt(thresh) / 2

    # set current error initially larger than stop criterion
    current_error = 2 * target_error

    # Residuals of the projected test vectors, deflated block by block below
    rem_test = Γmat' * Ω_test

    current_cols = 0

    # Iterate until convergence
    while current_error > target_error && current_cols < NG
        current_block_size = min(column_block_size, NG - current_cols)
        Ω = randn(T, Npp, current_block_size) # Draw a new random block
        Y_block = Γmat' * Ω                   # Project Γ onto Ω

        # Orthogonalize Y_block against existing Q: we do iterated Gram-Schmidt
        # to preserve orthogonality (assuming "twice is enough" rule)
        if !isempty(Q_blocks)
            # first pass
            for Qb in Q_blocks
                coeffs1 = Qb' * Y_block
                Y_block .-= Qb * coeffs1
            end
            # second pass
            for Qb in Q_blocks
                coeffs2 = Qb' * Y_block
                Y_block .-= Qb * coeffs2
            end
        end

        Q_block = Matrix(qr(Y_block).Q) # Orthonormalize block itself (QR)
        push!(Q_blocks, Q_block)        # Update stochastic basis blocks
        current_cols += current_block_size

        # Update current_error incrementally: worst residual over the test vectors
        rem_test .-= Q_block * (Q_block' * rem_test)
        current_error = maximum(norm, eachcol(rem_test))
    end

    # Combine blocks to form the full basis
    Q = reduce(hcat, Q_blocks)

    # === Compression Step ===
    Γ_proj = Γmat * Q                       # Project Γ onto Q
    H = -Hermitian(Γ_proj' * Γ_proj)        # Gramian in Q basis
    λ, U = eigen(H)                         # diagonalize
    NF = findlast(s -> abs(s) > thresh, λ)  # truncate based on thresh
    if isnothing(NF)
        return ΓmnG, LinearAlgebra.I(size(Γmat, 2))
    else
        coulomb_vertex_singular_vectors = Q * U[:, 1:NF]
        ΓmnF = Γ_proj * U[:, 1:NF]           # rotate
        return reshape(ΓmnF, size(ΓmnG)[1:4]..., NF), coulomb_vertex_singular_vectors
    end
end
