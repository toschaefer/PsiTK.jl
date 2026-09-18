@doc raw"""
    compute_overlap_densities(
        bra_space::OrbitalSpace,
        ket_space::OrbitalSpace;
        n_bands_bra=size(bra_space.ψ[1], 2),
        n_bands_ket=size(ket_space.ψ[1], 2),
        Ecut_ratio=1.0,
        callback=identity
    )
    compute_overlap_densities(space::OrbitalSpace; n_bands=size(space.ψ[1], 2), kwargs...)

Compute the overlap densities in reciprocal space
```math
ρ_{mn \bm G} = \int_Ω \; \psi_{m}(\bm r)^∗ \psi_{n}(\bm r)  \; e^{-i\bm r \bm G}  \; d^3 r
```

# Arguments
- `bra_space`: the bra orbital space (e.g. occupied space)
- `ket_space`: the ket orbital space (e.g. virtual space)
- `space`: a single orbital space used as both bra and ket. Only the upper triangle
  ``m ≤ n`` is computed, the rest follows from ``ρ_{nm,-\bm G} = ρ_{mn \bm G}^∗``.
  This shortcut is taken whenever `bra_space === ket_space`.
- `n_bands_bra`: number of bands to be considered from bra_space
- `n_bands_ket`: number of bands to be considered from ket_space
- `n_bands`: number of bands to be considered from `space` (bra and ket alike)
- `Ecut_ratio`: ratio of the plane-wave cutoff (in energy) for the densities relative to
  the orbital cutoff `basis.Ecut` (default: 1.0). Values up to `supersampling^2` of the
  basis (4 for DFTK's default) are allowed since the FFT grid holds products of orbitals
  exactly up to there; `Ecut_ratio=4` therefore yields the exact overlap densities.
- `callback`: called after each orbital pair with `(; step, total_steps)`,
  e.g. `callback=ShowProgress()` for a progress bar (default: no output)

# Returns
A tuple `(ρmnG, G_vectors)`:
- `ρmnG`: the overlap densities as a tensor of shape `(nk, n_bands_bra, nk, n_bands_ket, nG)`.
- `G_vectors`: the corresponding plane-wave vectors.
"""
function compute_overlap_densities(
    bra_space::OrbitalSpace,
    ket_space::OrbitalSpace;
    n_bands_bra = size(bra_space.ψ[1], 2),
    n_bands_ket = size(ket_space.ψ[1], 2),
    Ecut_ratio = 1.0,
    callback = identity,
)
    basis = bra_space.basis
    all(kpt -> iszero(kpt.coordinate), basis.kpoints) ||
        error("Overlap densities are only implemented for Gamma-point calculations.")
    G_indices = _G_indices_within_cutoff(basis, Ecut_ratio)
    ρmnG = _compute_overlap_densities(
        basis,
        bra_space.ψ,
        ket_space.ψ;
        n_bands_bra,
        n_bands_ket,
        G_indices,
        callback,
    )

    return ρmnG, G_vectors(basis)[G_indices]
end

function compute_overlap_densities(space::OrbitalSpace; n_bands = size(space.ψ[1], 2), kwargs...)
    return compute_overlap_densities(space, space; n_bands_bra=n_bands, n_bands_ket=n_bands, kwargs...)
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
    return findall(G -> sum(abs2, recip_lattice * G) / 2 <= Ecut_reduced, vec(G_vectors(basis)))
end

@doc raw"""
    compute_coulomb_vertex(
        bra_space::OrbitalSpace,
        ket_space::OrbitalSpace;
        interaction_kernel=DFTK.Coulomb(DFTK.ProbeCharge()),
        n_bands_bra=size(bra_space.ψ[1], 2),
        n_bands_ket=size(ket_space.ψ[1], 2),
        Ecut_ratio=2/3,
        callback=identity
    )
    compute_coulomb_vertex(space::OrbitalSpace; n_bands=size(space.ψ[1], 2), kwargs...)

Compute the Coulomb vertex
```math
Γ_{mn \bm G} = \sqrt{v(\bm G)} \; ρ_{mn \bm G}
```
where $ρ_{mn \bm G}$ are the overlap densities (see [`compute_overlap_densities`](@ref))
and $v(\bm G)$ is the interaction kernel, e.g. the Coulomb potential
```math
v(\bm G) = \frac{4π}{\bm G^2}
```

# Arguments
- `bra_space`: the bra orbital space (e.g. occupied space)
- `ket_space`: the ket orbital space (e.g. virtual space)
- `space`: a single orbital space used as both bra and ket, exploiting the symmetry
  ``Γ_{nm,-\bm G} = Γ_{mn \bm G}^∗`` (see [`compute_overlap_densities`](@ref))
- `interaction_kernel`: the DFTK interaction kernel to use (default: Coulomb)
- `n_bands_bra`: number of bands to be considered from bra_space
- `n_bands_ket`: number of bands to be considered from ket_space
- `n_bands`: number of bands to be considered from `space` (bra and ket alike)
- `Ecut_ratio`: cutoff ratio for the vertex (default: 2/3), see [`compute_overlap_densities`](@ref)
- `callback`: called after each orbital pair with `(; step, total_steps)`,
  e.g. `callback=ShowProgress()` for a progress bar (default: no output)

# Returns
A tuple `(ΓmnG, G_vectors, kernel_fourier)`:
- `ΓmnG`: the Coulomb vertex tensor in the uncompressed plane-wave basis.
- `G_vectors`: the corresponding plane-wave vectors.
- `kernel_fourier`: the evaluated interaction kernel at the returned G vectors.
"""
function compute_coulomb_vertex(
    bra_space::OrbitalSpace,
    ket_space::OrbitalSpace;
    interaction_kernel = DFTK.Coulomb(DFTK.ProbeCharge()),
    n_bands_bra = size(bra_space.ψ[1], 2),
    n_bands_ket = size(ket_space.ψ[1], 2),
    Ecut_ratio = 2/3,
    callback = identity,
)
    ρmnG, G_vectors = compute_overlap_densities(
        bra_space,
        ket_space;
        n_bands_bra,
        n_bands_ket,
        Ecut_ratio,
        callback,
    )

    # Kernel on the full FFT cube for momentum transfer q = 0 (Gamma-only)
    basis = bra_space.basis
    G_indices = _G_indices_within_cutoff(basis, Ecut_ratio)
    kernel_cube = DFTK.compute_kernel_fourier(interaction_kernel, basis, basis.kpoints[1])
    kernel_fourier = kernel_cube[G_indices]

    # Γ = √v ⊙ ρ along the G axis; ρmnG is not needed anymore, so scale in place
    ΓmnG = ρmnG
    ΓmnG .*= reshape(sqrt.(kernel_fourier), 1, 1, 1, 1, :)

    return ΓmnG, G_vectors, kernel_fourier
end

function compute_coulomb_vertex(space::OrbitalSpace; n_bands = size(space.ψ[1], 2), kwargs...)
    return compute_coulomb_vertex(space, space; n_bands_bra=n_bands, n_bands_ket=n_bands, kwargs...)
end

# This function initially based on code of the experimental "cc4s" branch in DFTK written by Michael Herbst
function _compute_overlap_densities(
    basis,
    ψ_bra::AbstractVector{<:AbstractArray{T}},
    ψ_ket::AbstractVector{<:AbstractArray{T}};
    n_bands_bra = size(ψ_bra[1], 2),
    n_bands_ket = size(ψ_ket[1], 2),
    G_indices = eachindex(G_vectors(basis)),
    callback = identity,
) where {T}
    n_kpt = length(basis.kpoints)

    # === Create index to map each stored G to -G on the full FFT cube ===
    Gs = G_vectors(basis)
    G_to_idx = Dict(Gs[i] => i for i in eachindex(Gs))
    idx_minus_G = [G_to_idx[-Gs[i]] for i in G_indices]

    # allocate overlap densities
    ρmnG = zeros(complex(T), n_kpt, n_bands_bra, n_kpt, n_bands_ket, length(G_indices))

    is_symmetric = (ψ_bra === ψ_ket)
    if is_symmetric
        total_steps = (n_bands_bra*(n_bands_bra+1)÷2)*n_kpt^2 # only upper triangle of ρmnG
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
    @views for (ikn, kptn) in enumerate(basis.kpoints), n = 1:n_bands_ket
        # Prepare ψnk(r)
        ψnk_real = ifft(basis, kptn, ψ_ket[ikn][:, n])

        for (ikm, kptm) in enumerate(basis.kpoints)
            for m = 1:n_bands_bra
                # Compute upper triangle only (m <= n) if spaces are symmetric
                # The lower triangle is filled via Hermitian conjugation below.
                if is_symmetric && m > n
                    continue
                end

                # Prepare ψmk(r)
                # TODO: pre-calculate some of them (not all because virtual space can be large)
                ψmk_real = ifft(basis, kptm, ψ_bra[ikm][:, m])

                # Calculate overlap density ρ_mn(r) = ψm*(r)ψn(r) and FFT it on the full cube
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
    ρmnG
end





@doc raw"""
    compress_coulomb_vertex(ΓmnG::AbstractArray{T,5}, strategy)

Compress the Coulomb vertex along its plane-wave axis $\bm G$ into a smaller auxiliary
field index $F$,
```math
Γ_{mn F} = \sum_{\bm G} Γ_{mn \bm G} \, U_{\bm G F}
```
where the columns of the transformation $U$ span the dominant subspace of the Coulomb
Gramian $\Gamma^\dagger \Gamma$. How $U$ is determined depends on `strategy`:
- [`CoulombGramian`](@ref): exact diagonalization of the Gramian
- [`AdaptiveRandomizedSVD`](@ref): randomized range finder followed by diagonalization

# Arguments
- `ΓmnG`: the uncompressed Coulomb vertex as returned by [`compute_coulomb_vertex`](@ref)
- `strategy`: the compression strategy, carrying its own threshold

# Returns
A tuple `(ΓmnF, coulomb_vertex_singular_vectors)`:
- `ΓmnF`: the compressed vertex of shape `(nk, n_bands, nk, n_bands, NF)`
- `coulomb_vertex_singular_vectors`: the transformation matrix $U$ of shape `(NG, NF)`
"""
function compress_coulomb_vertex end

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
function compress_coulomb_vertex(
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
This finder iteratively increases the columns of Q (i.e. $N_F$) in steps of $2\sqrt{N_{pp}}$
and stops when the error for each of `n_test_vectors` stochastic test vectors $\omega_i$
```math
\varepsilon_i =  \Vert (1 - QQ^\dagger)\Gamma^\dagger \omega_i \Vert
```
is smaller than $\sqrt{\text{thresh}}/2$. With $r$ test vectors this estimator bounds the
true projection error with probability $1 - 10^{-r}$
[Halko, Martinsson, Tropp, SIAM Rev. **53**, 217 (2011), Lemma 4.1]; a single test vector
would stop the finder too early in a small fraction of runs.
"""
Base.@kwdef struct AdaptiveRandomizedSVD
    thresh::Float64 = 1e-6
    n_test_vectors::Int = 10
end
function compress_coulomb_vertex(
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
    column_block_size = round(Int, 2*Npp^0.5)

    # Stochastic test vectors for error estimation
    Ω_test = randn(T, Npp, strategy.n_test_vectors)

    # target error a little smaller than √thresh
    target_error = sqrt(thresh)/2

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
